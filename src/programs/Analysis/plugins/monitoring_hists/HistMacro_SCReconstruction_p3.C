// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/SCDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/SCDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/SCDeltaTVsPhi
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/SCDeltaTVsTheta

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Matched_ReconstructedPID/Pi-");
	TH1I* locHist_SCDeltaT_PiMinus = (TH1I*)gDirectory->Get("SCDeltaT"); //pi-
	TH2I* locHist_SCDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("SCDeltaTVsP"); //pi-
	TH2I* locHist_SCDeltaTVsPhi_PiMinus = (TH2I*)gDirectory->Get("SCDeltaTVsPhi"); //pi-
	TH2I* locHist_SCDeltaTVsTheta_PiMinus = (TH2I*)gDirectory->Get("SCDeltaTVsTheta"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("SCReconstruction_p3", "SCReconstruction_p3", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaT_PiMinus != NULL)
	{
		locHist_SCDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaTVsP_PiMinus != NULL)
	{
		locHist_SCDeltaTVsP_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsP_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsP_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsP_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsP_PiMinus->Draw("COLZ");
	}


	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaTVsPhi_PiMinus != NULL)
	{
		locHist_SCDeltaTVsPhi_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsPhi_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsPhi_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsPhi_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsPhi_PiMinus->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaTVsTheta_PiMinus != NULL)
	{
		locHist_SCDeltaTVsTheta_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsTheta_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsTheta_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsTheta_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsTheta_PiMinus->Draw("COLZ");
	}
}

