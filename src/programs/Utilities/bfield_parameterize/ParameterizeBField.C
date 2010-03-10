

#include "chebyshev_fit.C"

//---------------------
// ParameterizeBField
//---------------------
void ParameterizeBField(void)
{
	gROOT->Reset();
	TColor::CreateColorWheel();

	TFile *f = new TFile("bfield.root");
	TTree *Bfield = (TTree*)gROOT->FindObject("Bfield");

	Parameterize(Bfield, "Bz");
	Parameterize(Bfield, "Bx");
}

//-------------------
// Parameterize
//-------------------
void Parameterize(TTree *Bfield, const char *Bi)
{
	// Bounds and poly-order of fits in cm
	//Double_t z_bounds[] = {-20.0, +25.0, +67.0, +100.0, +130.0, +180.0, +236.0, 0.0};
	//Double_t z_bounds[] = {-20.0, -10.0, +5.0, +30.0, +50.0, +60.0, +80.0, +100.0, +120.0, +140.0, +160.0, +180.0, +200.0, +220.0, +240.0, 0.0};
	//Double_t z_bounds[] = {-50.0, +5.0, +50.0, +80.0, +150.0, +200.0, +250.0, +300.0, +350.0, +400.0, +450.0, +500.0, +550.0, +610.0, 0.0};
	//Double_t z_bounds[] = {-50.8, +63.5, +170.18, +254.0, +330.2, +457.2, +599.44, 0.0};
	Double_t z_bounds[] = {-30.0, -10.0, +10.0, +30.0, +67.0, +100.0, +125.0, +148.0, +175.0, 205.0, +236.0, 0.0};
	int order1 = 9;	// max order of poly for Bi vs. z fit (fit to B-field)
	int order2 = 9;	// max order of poly for Pi vs. r fit (fit to parameters of above)
	Double_t rmin_inches = 1.0;
	Double_t rmax_inches = 26.0;
	UInt_t Nr = (UInt_t)(26.0-1.0)+1;
	UInt_t Nsec = 0;
	while(z_bounds[Nsec+1]!=0.0){
		// Convert z_bounds to inches (historic)
		//z_bounds[Nsec]/=2.54;
		
		Nsec++;
	}

	

	// Open file to write results to for later extraction
	stringstream fname;
	fname<<"BfieldParameters_"<<Bi<<".root";
	TFile *fout = new TFile(fname.str().c_str(),"RECREATE");
	
	// Create a Postscript file that we can use to draw multiple plots to
	// (doesn't seem to work for PDF)
	TCanvas *c1 = new TCanvas("c1");
	stringstream ss;
	ss<<"BfieldParameters_"<<Bi<<".ps";
	string gfname = ss.str();
	c1->Print((gfname+"[").c_str());

	int colors[] = {kBlue, kRed, kGreen, kMagenta, kCyan};
	int Ncolors = 5;

	// Create histograms to hold parameters from fits. The naming scheme is:
	//
	// secS_pN   : holds coefficient of poly order "N" from Bi vs. z fit for section "S"
	//             as function of r. This will have one bin on the x-axis for each value
	//             of r fit.
	//
	// secS_ppN  : holds coefficients of poly fit to the corresponding secS_pN histo.
	//             described above. This will have order2+1 bins.
	//
	TH1D *phist = new TH1D("phist","Parameters", 26, 0.5*2.54, 26.5*2.54);
	phist->SetStats(0);
	phist->SetXTitle("R (cm)");
	phist->SetYTitle("Value of Chebyshev coefficient");
	
	// Create histogram to hold z boundaries
	TH1D *z_bounds_hist = new TH1D("z_bounds_hist", "Z boundaries for sections", Nsec+1, 0.5, Nsec+1.5);
	for(UInt_t sec=0; sec<=Nsec; sec++)z_bounds_hist->Fill(sec+1, z_bounds[sec]*2.54);

	// Create a histo to hold each parameter as a function of r for this section
	TH1D *sec_pars[50][50]; // sec_pars[section][p_order]  (not all elements will be used)
	for(UInt_t sec=1; sec<=Nsec; sec++){
		for(UInt_t i=0; i<=order1; i++){
			char hname[256];
			sprintf(hname, "sec%d_p%d", sec, i);
			sec_pars[sec][i] = (TH1D*)phist->Clone(hname);
			
			stringstream title;
			title<<"Parameter "<<i<<" for section "<<sec<<" "<<Bi;
			sec_pars[sec][i]->SetTitle(title.str().c_str());
		}
	}

	// Loop over values of r, finding the parameters for fitting each
	for(Double_t r_inches=rmin_inches; r_inches<=rmax_inches; r_inches+=1.0){

		// Define sections that will be fit and loop over them
		TF1 *func[50];		// save function pointers for each section so we can plot later 
		for(UInt_t sec=1; sec<=Nsec; sec++){
			if(z_bounds[sec]==0.0)break;

			Double_t r = 2.54*r_inches;
			Double_t zmin = z_bounds[sec-1];
			Double_t zmax = z_bounds[sec];

			// Fit Bi vs. z for this value of r
			stringstream funcname;
			funcname<<"sec"<<sec<<"_r"<<(UInt_t)r_inches;
			func[sec] = FitSection(order1, Bfield, Bi,  zmin, zmax, r, funcname.str().c_str());

			// Get parameters
			for(UInt_t i=0; i<=order1; i++){
				sec_pars[sec][i]->SetBinContent(sec_pars[sec][i]->FindBin(r) ,func[sec]->GetParameter(i+2));
			}
		}
		
		// Plot results
		stringstream cut;
		cut<<"(abs(r-"<<r<<")<"<<0.1<<")"; // select a single r bin
		stringstream varexp;
		varexp<<Bi<<":z";
		Bfield->SetMarkerStyle(8);
		Bfield->SetMarkerSize(0.7);
		Bfield->Draw(varexp.str().c_str(), cut.str().c_str());
		
		// Overlay fit results for all sections
		for(UInt_t sec=1; sec<=Nsec; sec++){
			func[sec]->SetLineColor(colors[(sec-1)%Ncolors]);
			
			func[sec]->Draw("same");
		}

		// Draw plot to output file
		c1->Update();
		c1->Print(gfname.c_str());
	}

	// Now we want to fit the parameter histograms and save the
	// parameters describing each. 
	TH1D *pphist = new TH1D("pphist","Parameterization of Parameters", order2+1, -0.5, (double)order2 + 0.5);
	pphist->SetStats(0);
	pphist->SetXTitle("Order of Chebyshev coefficient");
	pphist->SetYTitle("Value of Chebyshev coefficient");

	// Loop over sections
	for(UInt_t sec=1; sec<=Nsec; sec++){

		// Loop over 1st level parameters
		for(UInt_t i=0; i<=order1; i++){

			TF1* f1 = chebyshev_Fit(sec_pars[sec][i], order2);

			char hname[256];
			sprintf(hname, "sec%d_pp%d", sec, i);
			TH1D *sec_ppars = (TH1D*)pphist->Clone(hname);
			for(int j=0; j<=order2; j++)sec_ppars->SetBinContent(j, f1->GetParameter(j+2));
			
			// Draw plot to output file
			f1->SetLineColor(colors[(sec-1)%Ncolors]);
			sec_pars[sec][i]->SetMarkerStyle(8);
			sec_pars[sec][i]->SetMarkerSize(0.7);
			
			sec_pars[sec][i]->Draw("P");
			f1->Draw("same");
			c1->Print(gfname.c_str());
		}
	} 

	// Close out graphics file
	c1->Clear();
	c1->Print((gfname+"]").c_str());
	
	delete phist;
	delete pphist;
	
	fout->Write();
	delete fout;
}

//---------------------
// FitSection
//---------------------
TF1* FitSection(int order, TTree *Bfield, const char *Bi, Double_t z_min_inches, Double_t z_max_inches, Double_t r, const char *fname)
{
	cout<<"Fitting "<<fname<<" at r="<<r<<endl;

	// The following should be defined based on the grid spacing/location
	// of points in the Bfield tree. 
	Double_t delta_z = 2.54; // in cm
	Double_t delta_r = 2.54; // in cm

	// Determine proper histogram limits 
	Double_t z_min = z_min_inches*delta_z - delta_z/2.0;
	Double_t z_max = z_max_inches*delta_z + delta_z/2.0;
	Double_t epsilon = delta_z/1000.0; // used to prevent round-off problems
	UInt_t NbinsZ = floor((z_max - z_min + epsilon)/delta_z);

	// Note that here we explicitly set the limits to be -1 to +1 so that the
	// histogram conforms to the requirements of the chebyshev_FindBestFunction
	// function.
	TH1D *Bi_vs_znorm = new TH1D("Bi_vs_znorm", "Bi as function of modified z coordinate", NbinsZ, z_min, z_max);

	// Make varexp that properly converts from the lab z-coordinates,
	// to the "Chebyshev compatible" coordinates where z_min=-1 and
	// z_max=+1.
	stringstream cut;
	cut<<Bi<<"*(abs(r-"<<r<<")<"<<delta_r/10.0<<")"; // select a single r bin
	Bfield->Project("Bi_vs_znorm", "z", cut.str().c_str());

	TF1 *func = chebyshev_Fit(Bi_vs_znorm, order);
	TF1 *myfunc = (TF1*)func->Clone(fname);

	delete Bi_vs_znorm;
	
	return myfunc;
}


