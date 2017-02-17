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
		KPlusKMinus->GetXaxis()->SetRangeUser(0.8, 1.4);

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

		Double_t I = FitWrapper::FitWithBackground(PiPlusPiMinus, 0.770, 0.1, 0.3, 1.6);
		
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
	
		Double_t I = FitWrapper::FitWithBackground(PiPlusPiMinusPiZero, 0.782, 0.03, 0.42, 1.6);
		
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
