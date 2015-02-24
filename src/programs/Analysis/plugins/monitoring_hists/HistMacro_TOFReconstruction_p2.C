// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Proton/TOFDeltaT
// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Proton/TOFDeltaTVsP
// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Pi-/TOFDeltaT
// hnamepath: /Independent/Hist_DetectorMatchParams/ReconstructedPID/Pi-/TOFDeltaTVsP

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatchParams");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("ReconstructedPID/Proton");
	TH1I* locHist_TOFDeltaT_Proton = (TH1I*)gDirectory->Get("TOFDeltaT"); //proton
	TH2I* locHist_TOFDeltaTVsP_Proton = (TH2I*)gDirectory->Get("TOFDeltaTVsP"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_TOFDeltaT_PiMinus = (TH1I*)gDirectory->Get("TOFDeltaT"); //pi-
	TH2I* locHist_TOFDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("TOFDeltaTVsP"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TOFReconstruction_p2", "TOFReconstruction_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFDeltaT_Proton != NULL)
	{
		locHist_TOFDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_Proton->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFDeltaTVsP_Proton != NULL)
	{
		locHist_TOFDeltaTVsP_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaTVsP_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaTVsP_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaTVsP_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaTVsP_Proton->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFDeltaT_PiMinus != NULL)
	{
		locHist_TOFDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFDeltaTVsP_PiMinus != NULL)
	{
		locHist_TOFDeltaTVsP_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaTVsP_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaTVsP_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaTVsP_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaTVsP_PiMinus->Draw("COLZ");
	}
}

