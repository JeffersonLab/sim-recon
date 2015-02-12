// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaTVsPhi
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaTVsTheta

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Matched_ReconstructedPID/Proton");
	TH1I* locHist_SCDeltaT_Proton = (TH1I*)gDirectory->Get("SCDeltaT"); //proton
	TH2I* locHist_SCDeltaTVsP_Proton = (TH2I*)gDirectory->Get("SCDeltaTVsP"); //proton
	TH2I* locHist_SCDeltaTVsPhi_Proton = (TH2I*)gDirectory->Get("SCDeltaTVsPhi"); //proton
	TH2I* locHist_SCDeltaTVsTheta_Proton = (TH2I*)gDirectory->Get("SCDeltaTVsTheta"); //proton

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("SCReconstruction_p2", "SCReconstruction_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaT_Proton != NULL)
	{
		locHist_SCDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_Proton->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaTVsP_Proton != NULL)
	{
		locHist_SCDeltaTVsP_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsP_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsP_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsP_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsP_Proton->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaTVsPhi_Proton != NULL)
	{
		locHist_SCDeltaTVsPhi_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsPhi_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsPhi_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsPhi_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsPhi_Proton->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaTVsTheta_Proton != NULL)
	{
		locHist_SCDeltaTVsTheta_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsTheta_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_SCDeltaTVsTheta_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsTheta_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaTVsTheta_Proton->Draw("COLZ");
	}
}

