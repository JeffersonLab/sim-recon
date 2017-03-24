// hnamepath: /highlevel/TwoGammaMass
// hnamepath: /highlevel/PiPlusPiMinus
// hnamepath: /highlevel/KPlusKMinus
// hnamepath: /highlevel/PiPlusPiMinusPiZero
// hnamepath: /highlevel/L1bits_gtp
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

{

// This is a trick to get ROOT to use a function for a TF1 without
// defining it in global namespace. This is needed since ROOTSpy
// requires all macros to be nameless. The start of the actual
// macro starts after the class definition.
class FitWrapper{
	public:
	
		// The following provided by Eugene C. (see e-mails to David L.
		// on 3/1/2017 amd 3/2/2017).

		//....................................................
		// mass_bg_1
		//
		// Function for the background in mass plots 
		// a*(x-m0)^b*exp(-c*x)
		//....................................................
		static Double_t mass_bg_1(Double_t *xptr, Double_t *p)
		{
			Double_t x=xptr[0];
			Double_t bg=0;
			if(x-p[1] >0.) bg=p[0]*pow((x-p[1]),p[2])*exp(-p[3]*x);
			return bg;
		}

		//....................................................
		// peak_gs
		//
		// Gaussian - the first parameter is the integral
		// p[0]/sqrt(2pi)/p[2]*exp(-pow((x[0]-p[2])/p[3],2)/2.)
		//....................................................
		static Double_t peak_gs(Double_t *xptr, Double_t *p)
		{
			Double_t x=xptr[0];
			Double_t gs;
			gs=0.;
			if(abs(p[2])>1.E-20){
				gs=p[0]/sqrt(2.*3.1416)/p[2]*exp(-pow((x-p[1])/p[2],2)/2.);
			}
			return gs;
		}

		//....................................................
		// fun_peak_bg_1
		//
		// Mass fit: peaks and background
		//....................................................
		static Double_t fun_peak_bg_1(Double_t *xptr, Double_t *p)
		{
			Double_t x=xptr[0];
			Int_t n_peaks=Int_t(p[0]+0.1);

			Double_t fun = mass_bg_1(xptr,&p[1]);
			for(int i=0;i<n_peaks;i++){
				Int_t iflg=Int_t(p[5+4*i]+0.5); // <0 - reject the points close to the peak, =0 - ignore, =1 - Gauss 
				if(p[5+4*i]<0.) iflg=Int_t(p[5+4*i]-0.5);
				if(iflg == 1){
					fun+=peak_gs( xptr,&p[5+4*i+1]); // add a Gauss peak
				} else if (iflg < 0){
					// Reject the point if in a range of +/- 3*sigma
					if(x>p[5+4*i+2]-p[5+4*i+3]*3. && x<p[5+4*i+2]+p[5+4*i+3]*3.){
						TF1::RejectPoint();
						return 0.;
					}
				}
			}      

			return fun;
		}


		//-----------------------------------
		// FitPeaksWithBackgr
		//-----------------------------------
		static Double_t FitPeaksWithBackgr(
			TH1* h1,       // histogram pointer
			Double_t mass_thresh  = 0.139*2.+0.135,  // mass threshold of the reaction 
			Double_t bg_pos       = 1.2,             // some point of pure background (+/- 5 bins)
			TString type_peaks    = "G",             // types of the peaks:  G -Gaussian, also limits the number of peaks
			                                         //   There are places for 4 peaks, ="GGG" means that 3 Gaussian peak should be used
			Double_t peak_pos_1   = 0.78,            // 1-st peak position (mass)
			Double_t peak_width_1 = 0.03,            // 1-st peak width,  =0 - peak ignored  
			Double_t peak_pos_2   = 0.0,             // 2-nd peak position (mass)
			Double_t peak_width_2 = 0.0,             // 2-nd peak width,  =0 - peak ignored
			Double_t peak_pos_3   = 0.0,             // 2-nd peak position (mass)
			Double_t peak_width_3 = 0.0,             // 2-nd peak width,  =0 - peak ignored
			Double_t peak_pos_4   = 0.0,             // 2-nd peak position (mass)
			Double_t peak_width_4 = 0.0)             // 2-nd peak width,  =0 - peak ignored
		{

			if(h1 == NULL){
//				cout << "FitPeaksWithBackgr:  No histogram name " <<hname << endl;
				return 0;
			}
			TString hname = TString(h1->GetName());
			Int_t nbins=h1->GetNbinsX();
			if(nbins <= 0){
			  cout << "FitPeaksWithBackgr:  Histogram name " <<hname << " N bins "<<nbins << endl;
				return 0;
			}
			Double_t xmin = h1->GetXaxis()->GetXmin();
			Double_t xmax = h1->GetXaxis()->GetXmax();
			Double_t xbin = (xmax-xmin)/nbins;

			gPad->SetLeftMargin(0.12);
			h1->GetYaxis()->SetTitleOffset(0.75);
			h1->GetYaxis()->SetTitle(Form("Combinations / %5.1f MeV",xbin*1000.));
			h1->Draw("");

			if(nbins<20){
	   			  cout << "FitPeaksWithBackgr: No fit - too few bins " << nbins << endl;
	   			  return 0;
			}

			if(h1->GetEntries()<500){
			  cout << "FitPeaksWithBackgr: No fit - too few entries " << h1->GetEntries() << endl;
				return 0;
			}

			// Specify the peaks
			Double_t set_peaks[10][2];   // [i_peak,i_par] = pos,width
			Int_t cod_peaks[10];         // [i_peak] = a code for the peak, =1 - Gauss
			set_peaks[0][0]=peak_pos_1;
			set_peaks[0][1]=peak_width_1;
			set_peaks[1][0]=peak_pos_2;
			set_peaks[1][1]=peak_width_2;
			set_peaks[2][0]=peak_pos_3;
			set_peaks[2][1]=peak_width_3;
			set_peaks[3][0]=peak_pos_4;
			set_peaks[3][1]=peak_width_4;
			Int_t n_peaks=0;            // The number of peaks defined, limited by the length of the string type_peaks 
			for(int i=0; i<int(type_peaks.Length()) && i<4 ; i++){
			  //        for(int i=0; i<1; i++){
			  if(set_peaks[i][1]>0.){
				 set_peaks[n_peaks][0]=set_peaks[i][0]; // compress holes if any 
				 set_peaks[n_peaks][1]=set_peaks[i][1]; 
				 cod_peaks[n_peaks]=1; // Defult - Gauss  
				 if(char(type_peaks(i))=='G' || char(type_peaks(i))=='g') cod_peaks[n_peaks]=1; // room for other function  
				 n_peaks++;
			  }
			}

			// cout << "n_peaks "<<n_peaks<<" "<<type_peaks.Length()<<endl;
			// Range for fitting
			//	Double_t xminfit = mass_thresh>xmin ? mass_thresh:xmin; 
			Double_t xminfit = mass_thresh;
			if(xminfit<xmin) xminfit=xmin;
			Double_t xmaxfit = xmax;

			// Check the position of the specified background
			if( bg_pos - 5.*xbin < xminfit || bg_pos + 5.*xbin > xmaxfit){
			  cout << "FitPeaksWithBackgr: No fit - the specified bg position "<< bg_pos << 
				 " is outside of the fit range, or too close to the edge "<<endl;
			  return 0;
			}
			for(int i=0;i<n_peaks;i++){
			  if( bg_pos + 5.*xbin > set_peaks[i][0]-3.*set_peaks[i][1] &&
         			  bg_pos - 5.*xbin < set_peaks[i][0]+3.*set_peaks[i][1] ){
				 cout << "FitPeaksWithBackgr: No fit - the specified bg position "<< bg_pos << 
	   			" is too close to peak "<< i <<endl;
				 return 0;
			  }
			}


			// Define fit function
			Int_t n_par = 1+4+4*n_peaks; // number of parameters controlling the function
			TString fun_name = hname + TString("_fun");
			TF1 *ftf = (TF1 *)gROOT->FindObject(fun_name);
			if(ftf!=NULL){
			  if(ftf->GetNpar() != n_par){
				 cout << " Function exists with wrong parameter number " << ftf->GetNpar() << endl;
				 ftf->Delete();
				 ftf=NULL;
			  }
			}
			if(ftf == NULL) ftf = new TF1(fun_name, fun_peak_bg_1, 0., 0., n_par);
			ftf->SetRange(xminfit, xmaxfit);
			ftf->SetParName( 0, "Function code .");   // n_peaks+10*BG_type
			ftf->SetParName( 1, "Bkgnd Amp     .");
			ftf->SetParName( 2, "Bkgnd threshold");
			ftf->SetParName( 3, "Bkgnd power   .");
			ftf->SetParName( 4, "Bkgnd exponent.");
			for(int i=0;i<n_peaks;i++){
			  TString nam_p="Peak  ";
			  if(char(type_peaks(i))=='G' || char(type_peaks(i))=='g') nam_p="Gauss ";
			  nam_p.Append(Form("%d",i+1));
			  ftf->SetParName( 5+i*4+0,nam_p+TString(" code  ."));
			  ftf->SetParName( 5+i*4+1,nam_p+TString(" integr."));
			  ftf->SetParName( 5+i*4+2,nam_p+TString(" mean  ."));
			  ftf->SetParName( 5+i*4+3,nam_p+TString(" width ."));
   			  }

			// Set starting parameters. We initially fix the peak parameters
			// so we can fit the background first.

			ftf->FixParameter(0, Double_t(n_peaks));
			ftf->SetParameter(1, 1.0);
			ftf->FixParameter(2, mass_thresh);
			ftf->SetParameter(3, 2.0);
			ftf->SetParameter(4, 4.0);

			for(int i=0;i<n_peaks;i++){  // Set the defined peaks
			  ftf->FixParameter(5+i*4+0,-2.); // at first fit the BG and reject the area around the peaks
			  ftf->FixParameter(5+i*4+1, 0.); // integral
			  ftf->FixParameter(5+i*4+2, set_peaks[i][0]); // position
			  ftf->FixParameter(5+i*4+3, set_peaks[i][1]); // width
			}

			// Find the backround in 10 bins around the set point
			Double_t x1=bg_pos-5.*xbin;
			Double_t x2=bg_pos+5.*xbin;
			Int_t ix1=h1->GetXaxis()->FindBin(x1);
			Int_t ix2=h1->GetXaxis()->FindBin(x2);
			Double_t bg=h1->Integral(ix1,ix2)/(ix2-ix1+1); 
			Double_t fun=ftf->Eval((x1+x2)/2.);
			Double_t norm = 1.;
			if(fun>0.) norm=bg/fun;
			ftf->SetParameter(1, norm); // scale background function to match histo at the set point

   			  ftf->SetParLimits(3,0.3,4.);
   			  ftf->SetParLimits(4,0.5,8.);

			Int_t fit_status = h1->Fit(ftf, "0", "", xminfit, xmaxfit); // Fit the BG, the peaks excluded
   			  if(fit_status != 0 ){
			  cout << "FitPeaksWithBackgr: Fit of backgound failed, code "<< fit_status <<endl;
			  return 0;
			}

			//Add the peaks one by one to the fit, starting with the LAST peak in the list
			Int_t n_ok_fit=0;        // number of peaks fitted successfully 
			for(int i=n_peaks-1;i>=0;i--){ 
			  // Find the size of the peak - integrate the peak and subtract the BG
			  x1=set_peaks[i][0]-set_peaks[i][1]*3.; 
			  x2=set_peaks[i][0]+set_peaks[i][1]*3.; 
			  ix1=h1->GetXaxis()->FindBin(x1);
			  ix2=h1->GetXaxis()->FindBin(x2);
			  Double_t peak_tot=h1->Integral(ix1,ix2); // peak total

			  ftf->FixParameter(5+i*4,0.); // include the BG in the peak area, but ignore the peak
			  Double_t peak_bg=ftf->Integral(x1,x2)/xbin; // BG

			  ftf->SetParameter(5+i*4+1,(peak_tot-peak_bg)*xbin); 
			  // Release peak parameters
			  ftf->FixParameter(5+i*4,Double_t(cod_peaks[i])); // Gaussian peak
			  ftf->ReleaseParameter(5+i*4+1);
			  ftf->ReleaseParameter(5+i*4+2);
			  ftf->ReleaseParameter(5+i*4+3);
			  Double_t wd=ftf->GetParameter(5+i*4+3);
			  ftf->SetParLimits(5+i*4+3,wd*0.3,wd*3.);

			  fit_status = h1->Fit(ftf, "0", "", xminfit, xmaxfit); // Fit the BG + peak
			  if(fit_status != 0 ){
				 cout << "FitPeaksWithBackgr: Fit of peak+backgound failed, code "<< fit_status <<endl;
				 return 0;
			  }
			  n_ok_fit++;

			}

   			  ftf->Draw("same");

   			  // Copy parameters into new function for plotting background
			TString fun_name_bg = fun_name + TString("_bg");
   			  TF1 *fbg = (TF1 *)gROOT->FindObject(fun_name_bg);
   			  if(fbg!=NULL){
			  if(fbg->GetNpar() != n_par){
				 cout << "FitPeaksWithBackgr: Function exists with wrong parameter number " << fbg->GetNpar() << endl;
				 fbg->Delete();
				 fbg=NULL;
			  }
   			  }
   			  if(fbg == NULL) fbg = new TF1("fun_bg", fun_peak_bg_1, 0., 0., n_par); // For some reason Clone doesn't work right here!
   			  fbg->SetParameters(ftf->GetParameters());
			for(int i=0;i<n_peaks;i++){fbg->FixParameter(5+i*4,0.);} // Turn off the peaks 
   			  fbg->SetLineStyle(2);
   			  fbg->SetRange(xminfit, xmaxfit);
   			  fbg->Draw("same");
			//  canv1->Update();

			if(n_peaks != n_ok_fit){
			  cout << "FitPeaksWithBackgr: Fit of " << n_peaks << " failed: successes = " << n_ok_fit << endl;
			  return 0;
			}

			//	Double_t par_peaks[10][3];   // Peak parameters integr,pos,width
			cout << " Peak   Integral                 Center                   Width " << endl;
			Double_t peak_integral=0.;  //  Integral of the 1-st peak
			for(int i=0;i<n_peaks;i++){
			  Double_t par[3], epar[3];
			  for(int j=0;j<3;j++){
				 par[j]=ftf->GetParameter(5+i*4+1+j);
				 epar[j]=ftf->GetParError(5+i*4+1+j);
			  }
			  printf(" %d   %8.1f +/- %7.1f   ",i+1,par[0]/xbin,epar[0]/xbin);
			  printf("  %8.4f +/- %7.4f   ",par[1],epar[1]);
			  printf("  %8.4f +/- %7.4f ",par[2],epar[2]);
			  cout << endl;
			  if(i==0){
				 peak_integral=par[0]/xbin;
			  }
			}

			// Draw line at nominal peak position
			if(n_peaks>0){	  
			  Double_t ymax = 1.05*h1->GetMaximum();
			  Double_t ymin = 0.;
			  TLine lin;
			  lin.SetLineColor(kMagenta);
			  lin.SetLineWidth(1);
			  lin.DrawLine(set_peaks[0][0], ymin, set_peaks[0][0], ymax);

			  char str[256];
			  sprintf(str, "%d MeV", (int)(1000*set_peaks[0][0]));

			  TLatex latex;
			  latex.SetTextAngle(90.0);
			  latex.SetTextSize(0.030);
			  latex.SetTextAlign(21);
			  latex.SetTextColor(kMagenta);
			  latex.DrawLatex(set_peaks[0][0] - (xmax-xmin)*0.02, ymin*0.8+ymax*0.2, str);
			}

			return peak_integral;
		}

		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		// The following methods (gauss_bg1 and FitWithBackground) were used up until
		// 3/3/2017 when they were replaced by the above. (See e-mails from Eugene C.
		// to David L. on 3/1/2017 and 3/2/2017)

		//....................................................
		// gauss_bg1
		//
		// Fit function that allows excluded range to be
		// specified using last 2 parameters. Functional
		// form provided by E. Chudakov.
		//....................................................
		static Double_t gauss_bg1(Double_t *xptr, Double_t *p)
		{
			Double_t x = xptr[0];
			Double_t xexcl_min = p[7];
			Double_t xexcl_max = p[8];
			if( x>xexcl_min && x<xexcl_max ){
				TF1::RejectPoint();
				return 0;
			}

			Double_t signal = p[0]*TMath::Gaus(x, p[1],p[2]);
			Double_t bkgnd  = p[3]*pow((x-p[4]),p[5])*exp(-x*p[6]);

			return signal + bkgnd;
		}

		//-----------------------------------
		// FitWithBackground
		//-----------------------------------
		static Double_t FitWithBackground(
			TH1* h1,
			Double_t peak_pos,
			Double_t peak_width,
			Double_t mass_thresh,
			Double_t xmaxfit=0.0)
		{

			// If too few events then just plot histogram and return
			if(h1->GetEntries()<100){
				h1->Draw();
				return 0.0;
			}

			// Make unique names for signal and background functions
			// for each histogram fit so they can be displayed
			// simultaneously.
			char ftfname[256];
			char fbgname[256];
			sprintf(ftfname, "f%s_signal", h1->GetName());
			sprintf(fbgname, "f%s_bkgrnd", h1->GetName());

			// Define fit function
			TF1 *ftf = (TF1 *)gROOT->FindObject(ftfname);
			if(!ftf){
				ftf = new TF1(ftfname, gauss_bg1, 0.0, 0.0, 9);
				ftf->SetParName(0, "Gauss Amp");
				ftf->SetParName(1, "Gauss mean");
				ftf->SetParName(2, "Gauss sigma");
				ftf->SetParName(3, "Bkgnd Amp");
				ftf->SetParName(4, "Bkgnd offset");
				ftf->SetParName(5, "Bkgnd exponent");
				ftf->SetParName(6, "Bkgnd expo-rate");
				ftf->SetParName(7, "xmin excluded region");
				ftf->SetParName(8, "xmax excluded region");
			}

			// Threshold parameter is either mass thresh or histogram low edge
			Double_t xmin = h1->GetXaxis()->GetXmin();

			// Set starting parameters. We initially fix the peak parameters
			// so we can fit the background first.
			ftf->FixParameter(0, 0.5*h1->GetBinContent(h1->FindBin(peak_pos)));
			ftf->FixParameter(1, peak_pos);
			ftf->FixParameter(2, peak_width);
			ftf->SetParameter(3, 1.0);
			ftf->SetParameter(4, mass_thresh>xmin ? mass_thresh:xmin);
			ftf->SetParameter(5, 2.0);
			ftf->SetParameter(6, 4.0);

			// Limits of signal region to exclude from initial fit
			Double_t xexcl_1 = peak_pos - 3.0*peak_width;
			Double_t xexcl_2 = peak_pos + 3.0*peak_width;
			ftf->FixParameter(7, xexcl_1);  // Set excluded region min
			ftf->FixParameter(8, xexcl_2);  // Set excluded region max

			// Find limits of initial background fit and do it
			Double_t xminfit = mass_thresh;
			if(xmaxfit==0.0) xmaxfit = h1->GetXaxis()->GetXmax();
			Int_t minbin_int = h1->FindBin(xexcl_2);
			Int_t maxbin_int = minbin_int + 5; // Integrate 6 bins
			Double_t xmin_int = h1->GetBinLowEdge(minbin_int);   // low edge of integration region
			Double_t xmax_int = h1->GetBinLowEdge(maxbin_int+1); // high edge of integration region
			Double_t norm = h1->GetBinContent(minbin_int, maxbin_int)/h1->GetBinWidth(minbin_int)/ftf->Integral(xmin_int, xmax_int);
			ftf->SetParameter(3, norm); // scale background function to match histo integral near edge of excluded region
			h1->Fit(ftf, "0", "", xminfit, xmaxfit);

			// Release peak parameters and fit to full range
			ftf->ReleaseParameter(0);
			ftf->ReleaseParameter(1);
			ftf->ReleaseParameter(2);
			ftf->ReleaseParameter(9);
			ftf->FixParameter(7, -1.0E6);  // disable excluded region
			ftf->FixParameter(8, -1.0E6);  // disable excluded region
			h1->Fit(ftf, "", "", xminfit, xmaxfit);

			// Copy parameters into new function for plotting background
			TF1 *fbg = (TF1 *)gROOT->FindObject(fbgname);
			if(!fbg) fbg = new TF1(fbgname, gauss_bg1, xminfit, xmaxfit, ftf->GetNpar()); // For some reason Clone doesn't work right here!
			fbg->SetParameters(ftf->GetParameters());
			fbg->SetParameter(0, 0.0); // zero out peak	
			fbg->SetLineStyle(2);
			fbg->SetLineColor(kMagenta);
			fbg->Draw("same");

			// Draw line at nominal peak position
			double max = 1.05*h1->GetMaximum();
			TLine lin;
			lin.SetLineColor(kMagenta);
			lin.SetLineWidth(1);
			lin.DrawLine(peak_pos, 0.0, peak_pos, max);

			char str[256];
			sprintf(str, "%d MeV", (int)(1000*peak_pos));

			TLatex latex;
			latex.SetTextAngle(90.0);
			latex.SetTextSize(0.035);
			latex.SetTextAlign(21);
			latex.SetTextColor(kMagenta);
			latex.DrawLatex(peak_pos - 0.005, max/2.0, str);

			// Get number of signal particless
			Double_t I = ftf->Integral(xminfit, xmaxfit) - fbg->Integral(xminfit, xmaxfit);
			I /= h1->GetBinWidth(1);

			return I;
		}

};

	//------------------------- Macro starts here ------------------------

	vector<bool> trig(6, true); // triggers to include 

	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH1* TwoGammaMass        = (TH1*)gDirectory->Get("TwoGammaMass");
	TH1* PiPlusPiMinus       = (TH1*)gDirectory->Get("PiPlusPiMinus");
	TH1* KPlusKMinus         = (TH1*)gDirectory->Get("KPlusKMinus");
	TH1* PiPlusPiMinusPiZero = (TH1*)gDirectory->Get("PiPlusPiMinusPiZero");
	TH1* L1bits_gtp          = (TH1*)gDirectory->Get("L1bits_gtp");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("PID", "PID", 1200, 600); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	double Ntrig_tot = 0.0;
	if(L1bits_gtp){
		for(int itrig=1; itrig<=6; itrig++){
			if(trig[itrig-1]) Ntrig_tot += (double)L1bits_gtp->GetBinContent(itrig);
		}
	}	


	TLatex latex;

	//----------- Pi0 --------------
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TwoGammaMass != NULL)
	{

		TwoGammaMass->GetXaxis()->SetTitleSize(0.05);
		TwoGammaMass->GetYaxis()->SetTitleSize(0.05);
		TwoGammaMass->GetXaxis()->SetLabelSize(0.05);
		TwoGammaMass->GetYaxis()->SetLabelSize(0.035);
		TwoGammaMass->SetStats(0);

		// Fit to pi0 peak
		TF1 *fun = (TF1*)gDirectory->FindObjectAny("fun_pi0_fit");
		if(!fun)fun = new TF1("fun_pi0_fit", "gaus(0) + pol2(3)");

		// Fit once with fixed parameters to force finding of polynomial params
		fun->FixParameter(0, TwoGammaMass->GetBinContent(TwoGammaMass->FindBin(0.134))*0.1);
		fun->FixParameter(1, 0.134);
		fun->FixParameter(2, 0.01);
		fun->SetParameter(3, 0.0);
		fun->SetParameter(4, 0.0);
		fun->SetParameter(5, 0.0);

		// Region of interest for fit
		double lo = 0.080;
		double hi = 0.190;

		// Fit and Draw
		TwoGammaMass->Fit(fun, "", "", lo, hi);

		// Release gaussian parameters and fit again
		fun->ReleaseParameter(0);
		fun->ReleaseParameter(1);
		fun->ReleaseParameter(2);

		// Fit and Draw again (histogram and function)
		TwoGammaMass->Fit(fun, "", "", lo, hi);

		// Second function for drawing background
		TF1 *fun2 = (TF1*)gDirectory->FindObjectAny("fun_pi0_fit2");
		if(!fun2) fun2 = new TF1("fun_pi0_fit", "pol2(0)" , lo, hi);
		double pars[10];
		fun->GetParameters(pars);
		fun2->SetParameters(&pars[3]);
		fun2->SetLineColor(kMagenta);
		fun2->SetLineStyle(2);
		fun2->Draw("same");
		
		double max = 1.05*TwoGammaMass->GetMaximum();
		TLine lin;
		lin.SetLineColor(kMagenta);
		lin.SetLineWidth(1);
		lin.DrawLine(0.135, 0.0, 0.135, max);
		
		TLatex latex;
		latex.SetTextAngle(90.0);
		latex.SetTextSize(0.035);
		latex.SetTextAlign(21);
		latex.SetTextColor(kMagenta);
		latex.DrawLatex(0.131, max/2.0, "135 MeV");

		// Get number of pi0's
		double I = fun->GetParameter(0)*fun->GetParameter(2)*sqrt(TMath::TwoPi());
		I /= TwoGammaMass->GetBinWidth(1);
		char str[256];
		sprintf(str, "num. #pi^{o} : %g", I);

		latex.SetTextColor(kBlack);
		latex.SetTextAngle(0.0);
		latex.SetTextAlign(11);
		latex.SetTextSize(0.075);
		latex.DrawLatex(0.175, max*3.0/4.0, str);
		
		// Print rate per trigger
		if(Ntrig_tot>0.0){
			sprintf(str, "%3.1f per 1k triggers", I/Ntrig_tot*1000.0);
			latex.SetTextSize(0.06);
			latex.DrawLatex(0.3, max*0.65, str);
		}
		
		// Print number of L1 triggers
		latex.SetTextSize(0.05);
		latex.SetTextAlign(12);
		if(L1bits_gtp){
			sprintf(str, "trig bit 1 (FCAL/BCAL): %g", (double)L1bits_gtp->GetBinContent(1));
			latex.DrawLatex(0.4, max*0.5, str);
			sprintf(str, "trig bit 3 (BCAL): %g", (double)L1bits_gtp->GetBinContent(3));
			latex.DrawLatex(0.4, max*0.4, str);
			sprintf(str, "trig bit 4 (PS): %g", (double)L1bits_gtp->GetBinContent(4));
			latex.DrawLatex(0.4, max*0.3, str);
		}	
	}


	//----------- Phi --------------
	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(KPlusKMinus != NULL)
	{
		KPlusKMinus->GetXaxis()->SetTitleSize(0.05);
		KPlusKMinus->GetYaxis()->SetTitleSize(0.05);
		KPlusKMinus->GetXaxis()->SetLabelSize(0.05);
		KPlusKMinus->GetYaxis()->SetLabelSize(0.035);
		KPlusKMinus->SetStats(0);
		KPlusKMinus->GetXaxis()->SetRangeUser(0.8, 2.0);
		
		Double_t I = FitWrapper::FitPeaksWithBackgr(KPlusKMinus, 0.494*2, 1.8, "GG", 1.02, 0.01, 1.22, 0.05);

		if(I>0.0){
			char str[256];
			sprintf(str, "num. #phi : %g", I);

			double max = 1.05*KPlusKMinus->GetMaximum();
			latex.SetTextColor(kBlack);
			latex.SetTextAngle(0.0);
			latex.SetTextAlign(11);
			latex.SetTextSize(0.075);
			latex.DrawLatex(1.4, max*3.0/4.0, str);

			// Print rate per trigger
			if(Ntrig_tot>0.0){
				sprintf(str, "%3.3f per 1k triggers", I/Ntrig_tot*1000.0);
				latex.SetTextSize(0.06);
				latex.DrawLatex(1.4, max*0.65, str);
			}
		}

#if 0
		// Only do fit if there are at least 25 entries in the
		// bin at 1020MeV
		Int_t Npeak = KPlusKMinus->GetBinContent(KPlusKMinus->FindBin(1.020));
		double max = 1.05*KPlusKMinus->GetMaximum();
		if(Npeak < 25){
			KPlusKMinus->Draw();
		}else{
		
			// Fit to phi peak
			TF1 *fun = (TF1*)gDirectory->FindObjectAny("fun_phi_fit");
			if(!fun)fun = new TF1("fun_phi_fit", "[0]*TMath::Voigt(x-[1], [2], [3]) + pol2(4)");

			// Fit once with fixed parameters to force finding of polynomial params
			fun->FixParameter(0, Npeak*0.5);
			fun->FixParameter(1, 1.020);
			fun->FixParameter(2, 0.2);
			fun->FixParameter(3, 0.1);
			fun->SetParameter(4, 0.0);
			fun->FixParameter(5, 0.0);
			fun->SetParameter(6, 0.0);
			fun->SetParameter(7, 0.0);
			//fun->SetParameter(8, 0.0);

			// Region of interest for fit
			double lo = 0.98;
			double hi = 1.07;

			// Fit and Draw
			KPlusKMinus->Fit(fun, "", "", lo, hi);

			// Release Voigt parameters and fit again
			fun->ReleaseParameter(0);
			fun->ReleaseParameter(1);
			fun->ReleaseParameter(2);
			fun->ReleaseParameter(3);

			// Fit and Draw again (histogram and function)
			KPlusKMinus->Fit(fun, "", "", lo, hi);

			// Second function for drawing background
			TF1 *fun2 = (TF1*)gDirectory->FindObjectAny("fun_phi_fit2");
			if(!fun2) fun2 = new TF1("fun_phi_fit2", "pol3(0)" , lo, hi);
			double pars[10];
			fun->GetParameters(pars);
			fun2->SetParameters(&pars[4]);
			fun2->SetLineColor(kMagenta);
			fun2->SetLineStyle(2);
			fun2->Draw("same");

			// Get number of rho's
			double I = fun->Integral(lo, hi) - fun2->Integral(lo,hi);
			I /= TwoGammaMass->GetBinWidth(1);
			char str[256];
			sprintf(str, "num. #phi : %g", I);

			latex.SetTextColor(kBlack);
			latex.SetTextAngle(0.0);
			latex.SetTextAlign(11);
			latex.SetTextSize(0.075);
			latex.DrawLatex(0.81, max*0.93, str);

			// Print rate per trigger
			if(Ntrig_tot>0.0){
				sprintf(str, "%3.3f per 1k triggers", I/Ntrig_tot*1000.0);
				latex.SetTextSize(0.06);
				latex.DrawLatex(0.81, max*0.85, str);
			}
		}

		TLine lin;
		lin.SetLineColor(kMagenta);
		lin.SetLineWidth(1);
		lin.DrawLine(1.020, 0.0, 1.020, max);
		
		latex.SetTextAngle(90.0);
		latex.SetTextSize(0.035);
		latex.SetTextAlign(21);
		latex.SetTextColor(kMagenta);
		latex.DrawLatex(1.015, max/2.0, "1020 MeV");
#endif
	}

	//----------- Rho --------------
	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(PiPlusPiMinus != NULL)
	{
		PiPlusPiMinus->GetXaxis()->SetTitleSize(0.05);
		PiPlusPiMinus->GetYaxis()->SetTitleSize(0.05);
		PiPlusPiMinus->GetXaxis()->SetLabelSize(0.05);
		PiPlusPiMinus->GetYaxis()->SetLabelSize(0.035);
		PiPlusPiMinus->SetStats(0);

		//Double_t I = FitWrapper::FitWithBackground(PiPlusPiMinus, 0.770, 0.1, 0.3, 1.6);
		Double_t I = FitWrapper::FitPeaksWithBackgr(PiPlusPiMinus, 0.139*2, 1.8, "G", 0.770, 0.085);
		
		if(I>0.0){
			char str[256];
			sprintf(str, "num. #rho : %g", I);

			double max = 1.05*PiPlusPiMinus->GetMaximum();
			latex.SetTextColor(kBlack);
			latex.SetTextAngle(0.0);
			latex.SetTextAlign(11);
			latex.SetTextSize(0.075);
			latex.DrawLatex(1.005, max*3.0/4.0, str);

			// Print rate per trigger
			if(Ntrig_tot>0.0){
				sprintf(str, "%3.3f per 1k triggers", I/Ntrig_tot*1000.0);
				latex.SetTextSize(0.06);
				latex.DrawLatex(1.010, max*0.65, str);
			}
		}
	}

	//----------- Omega --------------
	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(PiPlusPiMinusPiZero != NULL)
	{
		PiPlusPiMinusPiZero->GetXaxis()->SetTitleSize(0.05);
		PiPlusPiMinusPiZero->GetYaxis()->SetTitleSize(0.05);
		PiPlusPiMinusPiZero->GetXaxis()->SetLabelSize(0.05);
		PiPlusPiMinusPiZero->GetYaxis()->SetLabelSize(0.035);
		PiPlusPiMinusPiZero->SetStats(0);
	
		//Double_t I = FitWrapper::FitWithBackground(PiPlusPiMinusPiZero, 0.782, 0.03, 0.42, 1.6);
		Double_t I = FitWrapper::FitPeaksWithBackgr(PiPlusPiMinusPiZero, 0.139*2+0.135, 1.3, "G", 0.782, 0.009);
		
		if(I>0.0){
			char str[256];
			sprintf(str, "num. #omega : %g", I);

			double max = 1.05*PiPlusPiMinusPiZero->GetMaximum();
			latex.SetTextColor(kBlack);
			latex.SetTextAngle(0.0);
			latex.SetTextAlign(11);
			latex.SetTextSize(0.075);
			latex.DrawLatex(1.005, max*0.8, str);

			// Print rate per trigger
			if(Ntrig_tot>0.0){
				sprintf(str, "%3.3f per 1k triggers", I/Ntrig_tot*1000.0);
				latex.SetTextSize(0.06);
				latex.DrawLatex(1.010, max*0.65, str);
			}
		}

	}
}
