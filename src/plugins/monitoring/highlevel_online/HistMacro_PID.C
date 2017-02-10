// hnamepath: /highlevel/TwoGammaMass
// hnamepath: /highlevel/PiPlusPiMinus
// hnamepath: /highlevel/PiPlusPiMinusPiZero
// hnamepath: /highlevel/BetaVsP
// hnamepath: /highlevel/L1bits_gtp
//
// e-mail: davidl@jlab.org
// e-mail: pmatt@jlab.org
// e-mail: staylor@jlab.org
// e-mail: sdobbs@jlab.org
//

//TwoGammaMass
//PiPlusPiMinus
//PiPlusPiMinusPiZero
//BetaVsP

{
	vector<bool> trig(6, true); // triggers to include 

	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("highlevel");
	if(!locDirectory)
		return;
	locDirectory->cd();

	TH1* TwoGammaMass        = (TH1*)gDirectory->Get("TwoGammaMass");
	TH1* PiPlusPiMinus       = (TH1*)gDirectory->Get("PiPlusPiMinus");
	TH1* PiPlusPiMinusPiZero = (TH1*)gDirectory->Get("PiPlusPiMinusPiZero");
	TH2* BetaVsP             = (TH2*)gDirectory->Get("BetaVsP");
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

	//Draw
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

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(BetaVsP != NULL)
	{
		BetaVsP->GetXaxis()->SetTitleSize(0.05);
		BetaVsP->GetYaxis()->SetTitleSize(0.045);
		BetaVsP->GetXaxis()->SetLabelSize(0.05);
		BetaVsP->GetYaxis()->SetLabelSize(0.05);
		BetaVsP->SetStats(0);
		BetaVsP->Draw("colz");
		//gPad->SetLogz();
		gPad->Update();
	}

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
		
		// Fit to rho0 peak
		TF1 *fun = (TF1*)gDirectory->FindObjectAny("fun_rho0_fit");
		if(!fun)fun = new TF1("fun_rho0_fit", "[0]*TMath::Voigt(x-[1], [2], [3]) + pol2(4)");

		// Fit once with fixed parameters to force finding of polynomial params
		fun->FixParameter(0, PiPlusPiMinus->GetBinContent(PiPlusPiMinus->FindBin(0.770))*0.5);
		fun->FixParameter(1, 0.770);
		fun->FixParameter(2, 0.2);
		fun->FixParameter(3, 0.1);
		fun->SetParameter(4, 0.0);
		fun->FixParameter(5, 0.0);
		fun->SetParameter(6, 0.0);
		fun->SetParameter(7, 0.0);
		//fun->SetParameter(8, 0.0);

		// Region of interest for fit
		double lo = 0.4;
		double hi = 1.1;

		// Fit and Draw
		PiPlusPiMinus->Fit(fun, "", "", lo, hi);

		// Release Voigt parameters and fit again
		fun->ReleaseParameter(0);
		fun->ReleaseParameter(1);
		fun->ReleaseParameter(2);
		fun->ReleaseParameter(3);

		// Fit and Draw again (histogram and function)
		PiPlusPiMinus->Fit(fun, "", "", lo, hi);

		// Second function for drawing background
		TF1 *fun2 = (TF1*)gDirectory->FindObjectAny("fun_rho0_fit2");
		if(!fun2) fun2 = new TF1("fun_rho0_fit2", "pol3(0)" , lo, hi);
		double pars[10];
		fun->GetParameters(pars);
		fun2->SetParameters(&pars[4]);
		fun2->SetLineColor(kMagenta);
		fun2->SetLineStyle(2);
		fun2->Draw("same");

		double max = 1.05*PiPlusPiMinus->GetMaximum();
		TLine lin;
		lin.SetLineColor(kMagenta);
		lin.SetLineWidth(1);
		lin.DrawLine(0.770, 0.0, 0.770, max);
		
		TLatex latex;
		latex.SetTextAngle(90.0);
		latex.SetTextSize(0.035);
		latex.SetTextAlign(21);
		latex.SetTextColor(kMagenta);
		latex.DrawLatex(0.765, max/2.0, "770 MeV");

		// Get number of rho's
		double I = fun->Integral(lo, hi) - fun2->Integral(lo,hi);
		I /= TwoGammaMass->GetBinWidth(1);
		char str[256];
		sprintf(str, "num. #rho : %g", I);

		latex.SetTextColor(kBlack);
		latex.SetTextAngle(0.0);
		latex.SetTextAlign(11);
		latex.SetTextSize(0.075);
		latex.DrawLatex(1.005, max*3.0/4.0, str);
		
		// Print rate per trigger
		if(Ntrig_tot>0.0){
			sprintf(str, "%3.3f per 1k triggers", I/Ntrig_tot*1000.0);
			latex.SetTextSize(0.06);
			latex.DrawLatex(1.25, max*0.65, str);
		}
	}

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

		// Fit to rho0 peak
		TF1 *fun = (TF1*)gDirectory->FindObjectAny("fun_omega_fit");
		if(!fun)fun = new TF1("fun_omega_fit", "[0]*TMath::Voigt(x-[1], [2], [3]) + pol3(4)");

		// Fit once with fixed parameters to force finding of polynomial params
		fun->FixParameter(0, PiPlusPiMinusPiZero->GetBinContent(PiPlusPiMinusPiZero->FindBin(0.782))*0.5);
		fun->FixParameter(1, 0.782);
		fun->FixParameter(2, 0.2);
		fun->FixParameter(3, 0.1);
		fun->SetParameter(4, 0.0);
		fun->FixParameter(5, 0.0);
		fun->SetParameter(6, 0.0);
		fun->SetParameter(7, 0.0);
		//fun->SetParameter(8, 0.0);

		// Region of interest for fit
		double lo = 0.6;
		double hi = 1.0;

		// Fit and Draw
		PiPlusPiMinusPiZero->Fit(fun, "", "", lo, hi);

		// Release Voigt parameters and fit again
		fun->ReleaseParameter(0);
		fun->ReleaseParameter(1);
		fun->ReleaseParameter(2);
		fun->ReleaseParameter(3);

		// Fit and Draw again (histogram and function)
		PiPlusPiMinusPiZero->Fit(fun, "", "", lo, hi);

		// Second function for drawing background
		TF1 *fun2 = (TF1*)gDirectory->FindObjectAny("fun_omega_fit2");
		if(!fun2) fun2 = new TF1("fun_omega_fit2", "pol3(0)" , lo, hi);
		double pars[20];
		fun->GetParameters(pars);
		fun2->SetParameters(&pars[4]);
		fun2->SetLineColor(kMagenta);
		fun2->SetLineStyle(2);
		fun2->Draw("same");
		
		double max = 1.05*PiPlusPiMinusPiZero->GetMaximum();
		TLine lin;
		lin.SetLineColor(kMagenta);
		lin.SetLineWidth(1);
		lin.DrawLine(0.782, 0.0, 0.782, max);
		
		TLatex latex;
		latex.SetTextAngle(90.0);
		latex.SetTextSize(0.035);
		latex.SetTextAlign(21);
		latex.SetTextColor(kMagenta);
		latex.DrawLatex(0.777, max/2.0, "782 MeV");

		// Get number of omega's
		double I = fun->Integral(lo, hi) - fun2->Integral(lo,hi);
		I /= TwoGammaMass->GetBinWidth(1);
		char str[256];
		sprintf(str, "num. #omega : %g", I);

		latex.SetTextColor(kBlack);
		latex.SetTextAngle(0.0);
		latex.SetTextAlign(11);
		latex.SetTextSize(0.075);
		latex.DrawLatex(1.005, max*0.8, str);
		
		// Print rate per trigger
		if(Ntrig_tot>0.0){
			sprintf(str, "%3.3f per 1k triggers", I/Ntrig_tot*1000.0);
			latex.SetTextSize(0.06);
			latex.DrawLatex(1.25, max*0.72, str);
		}
	}
}
