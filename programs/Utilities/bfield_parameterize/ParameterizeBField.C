

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
	TCanvas *c1 = new TCanvas("c1");
	
	// Open file to write results to for later extraction
	stringstream fname;
	fname<<"BfieldParameters_"<<Bi<<".root";
	TFile *fout = new TFile(fname.str().c_str(),"RECREATE");
	
	// Make histograms to hold parameters as a function of r (in cm).
	TH1D *sec1_pars[10];
	TH1D *sec2_pars[10];
	TH1D *sec3_pars[10];
	TH1D *phist = new TH1D("phist","Parameters", 26, 0.5*2.54, 26.5*2.54);
	phist->SetStats(0);
	phist->SetXTitle("R (cm)");
	phist->SetYTitle("Value of Chebyshev coefficient");
	for(UInt_t i=0; i<10; i++){
		char hname[256];
		sprintf(hname, "sec1_p%d", i);
		sec1_pars[i] = (TH1D*)phist->Clone(hname);
		sprintf(hname, "sec2_p%d", i);
		sec2_pars[i] = (TH1D*)phist->Clone(hname);
		sprintf(hname, "sec3_p%d", i);
		sec3_pars[i] = (TH1D*)phist->Clone(hname);
	}
	
	// Create a Postscript file that we can use to draw multiple plots to
	stringstream ss;
	ss<<"BfieldParameters_"<<Bi<<".ps";
	string gfname = ss.str();
	c1->Print((gfname+"[").c_str());
	
	// Loop over values of r, finding the parameters for fitting each
	for(Double_t r_inches=1.0; r_inches<=26.0; r_inches+=1.0){

		// A single fit does not seem to reproduce the Bx field well
		// so it is done piece-wise with sections roughly covering
		// the CDC(-20 to 67), the FDC(67 to 138) and downstream to
		// the FCAL(138 to 236). The 236 limit seems necessary since
		// the current maps drops off there and the zeros cause bad
		// flucuations in the Chebyshev fit.
		Double_t r = 2.54*r_inches;
		TF1* f1 = FitSection(Bfield, Bi,  -20.0,  +67.0, r, "sec1");
		TF1* f2 = FitSection(Bfield, Bi,  +67.0, +138.0, r, "sec2");
		TF1* f3 = FitSection(Bfield, Bi, +138.0, +236.0, r, "sec3");
		
		// Plot results
		stringstream cut;
		cut<<"(abs(r-"<<r<<")<"<<0.1<<")"; // select a single r bin
		stringstream varexp;
		varexp<<Bi<<":z";
		Bfield->SetMarkerStyle(8);
		Bfield->SetMarkerSize(0.7);
		Bfield->Draw(varexp.str().c_str(), cut.str().c_str());
		
		f1->SetLineColor(kBlue);
		f2->SetLineColor(kRed);
		f3->SetLineColor(kGreen);
		
		f1->Draw("same");
		f2->Draw("same");
		f3->Draw("same");
		c1->Update();
		
		// Draw plot to output file
		c1->Print(gfname.c_str());
		
		// Get parameters
		for(UInt_t i=0; i<10; i++){
			sec1_pars[i]->SetBinContent(sec1_pars[i]->FindBin(r) ,f1->GetParameter(i+2));
			sec2_pars[i]->SetBinContent(sec2_pars[i]->FindBin(r) ,f2->GetParameter(i+2));
			sec3_pars[i]->SetBinContent(sec3_pars[i]->FindBin(r) ,f3->GetParameter(i+2));
		}
	}
	
	// Now we want to fit the parameter histograms and save the
	// parameters describing each. These will be saved in histograms
	// in the output ROOT file.
	UInt_t order = 9;
	TH1D *pphist = new TH1D("pphist","Parameterization of Parameters", order+1, -0.5, (double)order + 0.5);
	pphist->SetStats(0);
	pphist->SetXTitle("Order of Chebyshev coefficient");
	pphist->SetYTitle("Value of Chebyshev coefficient");
	for(UInt_t i=0; i<10; i++){

		char hname[256];
		sprintf(hname, "sec1_pp%d", i);
		TH1D *sec1_pars_pars = (TH1D*)pphist->Clone(hname);
		TF1* f1 = chebyshev_Fit(sec1_pars[i], order);	
		for(int j=0; j<=order; j++)sec1_pars_pars->SetBinContent(j, f1->GetParameter(j+2));

		sprintf(hname, "sec2_pp%d", i);
		TH1D *sec2_pars_pars = (TH1D*)pphist->Clone(hname);
		TF1* f2 = chebyshev_Fit(sec2_pars[i], order);	
		for(int j=0; j<=order; j++)sec2_pars_pars->SetBinContent(j, f2->GetParameter(j+2));

		sprintf(hname, "sec3_pp%d", i);
		TH1D *sec3_pars_pars = (TH1D*)pphist->Clone(hname);
		TF1* f3 = chebyshev_Fit(sec3_pars[i], order);	
		for(int j=0; j<=order; j++)sec3_pars_pars->SetBinContent(j, f3->GetParameter(j+2));
		
		// Draw plots of each to output file
		f1->SetLineColor(kBlue);
		f2->SetLineColor(kRed);
		f3->SetLineColor(kGreen);

		sec1_pars[i]->SetMarkerStyle(8);
		sec2_pars[i]->SetMarkerStyle(8);
		sec3_pars[i]->SetMarkerStyle(8);
		sec1_pars[i]->SetMarkerSize(0.7);
		sec2_pars[i]->SetMarkerSize(0.7);
		sec3_pars[i]->SetMarkerSize(0.7);
		
		sec1_pars[i]->Draw("P");
		f1->Draw("same");
		c1->Print(gfname.c_str());

		sec2_pars[i]->Draw("P");
		f2->Draw("same");
		c1->Print(gfname.c_str());

		sec3_pars[i]->Draw("P");
		f3->Draw("same");
		c1->Print(gfname.c_str());
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
TF1* FitSection(TTree *Bfield, const char *Bi, Double_t z_min_inches, Double_t z_max_inches, Double_t r, const char *fname)
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

	TF1 *func = chebyshev_Fit(Bi_vs_znorm, 9);
	TF1 *myfunc = (TF1*)func->Clone(fname);

	delete Bi_vs_znorm;
	
	return myfunc;
}


