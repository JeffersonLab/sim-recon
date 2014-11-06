// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/FCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/FCALShowerDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/FCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/FCALShowerDeltaTVsP

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Matched_ReconstructedPID/Proton");
	TH1I* locHist_FCALShowerDeltaT_Proton = (TH1I*)gDirectory->Get("FCALShowerDeltaT"); //proton
	TH2I* locHist_FCALShowerDeltaTVsP_Proton = (TH2I*)gDirectory->Get("FCALShowerDeltaTVsP"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_FCALShowerDeltaT_PiMinus = (TH1I*)gDirectory->Get("FCALShowerDeltaT"); //pi-
	TH2I* locHist_FCALShowerDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("FCALShowerDeltaTVsP"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("FCALReconstruction_p2", "FCALReconstruction_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw

	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaT_Proton != NULL)
	{
		locHist_FCALShowerDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_Proton->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaT_PiMinus != NULL)
	{
		locHist_FCALShowerDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaTVsP_Proton != NULL)
	{
		locHist_FCALShowerDeltaTVsP_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_Proton->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaTVsP_PiMinus != NULL)
	{
		locHist_FCALShowerDeltaTVsP_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaTVsP_PiMinus->Draw("COLZ");
	}
}

