// hnamepath: /highlevel/TwoGammaMass
// hnamepath: /highlevel/PiPlusPiMinus
// hnamepath: /highlevel/PiPlusPiMinusPiZero
// hnamepath: /highlevel/BetaVsP

//TwoGammaMass
//PiPlusPiMinus
//PiPlusPiMinusPiZero
//BetaVsP

{
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

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("PID", "PID", 1200, 600); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);


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
		TwoGammaMass->Draw();
		
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
		PiPlusPiMinus->Draw();
		
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
		PiPlusPiMinusPiZero->Draw();
		
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
	}
}
