// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/BCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/BCALShowerDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/BCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/BCALShowerDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/PID/q-/BetaVsP_BCAL
// hnamepath: /Independent/Hist_DetectorStudies/PID/q+/BetaVsP_BCAL

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Matched_ReconstructedPID/Proton");
	TH1I* locHist_BCALShowerDeltaT_Proton = (TH1I*)gDirectory->Get("BCALShowerDeltaT"); //proton
	TH2I* locHist_BCALShowerDeltaTVsP_Proton = (TH2I*)gDirectory->Get("BCALShowerDeltaTVsP"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_BCALShowerDeltaT_PiMinus = (TH1I*)gDirectory->Get("BCALShowerDeltaT"); //pi-
	TH2I* locHist_BCALShowerDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("BCALShowerDeltaTVsP"); //pi-

	gDirectory->cd("../../PID/q-");
	TH2I* locHist_BetaVsP_QMinus = (TH2I*)gDirectory->Get("BetaVsP_BCAL"); //q-
	gDirectory->cd("../q+");
	TH2I* locHist_BetaVsP_QPlus = (TH2I*)gDirectory->Get("BetaVsP_BCAL"); //q-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("BCALReconstruction_p2", "BCALReconstruction_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerDeltaT_Proton != NULL)
	{
		locHist_BCALShowerDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_Proton->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerDeltaT_PiMinus != NULL)
	{
		locHist_BCALShowerDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerDeltaTVsP_Proton != NULL)
	{
		locHist_BCALShowerDeltaTVsP_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaTVsP_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaTVsP_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaTVsP_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaTVsP_Proton->Draw("COLZ");
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerDeltaTVsP_PiMinus != NULL)
	{
		locHist_BCALShowerDeltaTVsP_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaTVsP_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaTVsP_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaTVsP_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaTVsP_PiMinus->Draw("COLZ");
	}

	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QPlus != NULL)
	{
		locHist_BetaVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QMinus != NULL)
	{
		locHist_BetaVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->Draw("COLZ");
	}
}

