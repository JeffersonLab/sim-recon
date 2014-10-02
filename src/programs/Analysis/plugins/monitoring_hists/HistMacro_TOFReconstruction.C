// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/TOFPointYVsX
// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/TOFPointEnergy
// hnamepath: /Independent/Hist_DetectorStudies/dEdxVsPByCharge/q-/TOFdEdXVsP
// hnamepath: /Independent/Hist_DetectorStudies/dEdxVsPByCharge/q+/TOFdEdXVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/TOFDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/TOFDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/TOFDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/TOFDeltaTVsP

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Reconstruction");
	TH2I* locHist_TOFPointYVsX = (TH2I*)gDirectory->Get("TOFPointYVsX");
	TH1I* locHist_TOFPointEnergy = (TH1I*)gDirectory->Get("TOFPointEnergy");

	gDirectory->cd("../dEdxVsPByCharge/q-");
	TH2I* locHist_TOFdEdXVsP_QMinus = (TH2I*)gDirectory->Get("TOFdEdXVsP"); //q-
	gDirectory->cd("../q+");
	TH2I* locHist_TOFdEdXVsP_QPlus = (TH2I*)gDirectory->Get("TOFdEdXVsP"); //q-

	gDirectory->cd("../../Matched_ReconstructedPID/Proton");
	TH1I* locHist_TOFDeltaT_Proton = (TH1I*)gDirectory->Get("TOFDeltaT"); //proton
	TH2I* locHist_TOFDeltaTVsP_Proton = (TH2I*)gDirectory->Get("TOFDeltaTVsP"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_TOFDeltaT_PiMinus = (TH1I*)gDirectory->Get("TOFDeltaT"); //pi-
	TH2I* locHist_TOFDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("TOFDeltaTVsP"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TOFReconstruction"); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 3);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFPointYVsX != NULL)
	{
		locHist_TOFPointYVsX->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFPointYVsX->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFPointYVsX->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFPointYVsX->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFPointYVsX->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFPointEnergy != NULL)
	{
		locHist_TOFPointEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFPointEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFPointEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFPointEnergy->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFDeltaT_Proton != NULL)
	{
		locHist_TOFDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_Proton->Draw();
	}

	locCanvas->cd(5);
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

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFdEdXVsP_QPlus != NULL)
	{
		locHist_TOFdEdXVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(7);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFDeltaT_PiMinus != NULL)
	{
		locHist_TOFDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(8);
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

	locCanvas->cd(9);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFdEdXVsP_QMinus != NULL)
	{
		locHist_TOFdEdXVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QMinus->Draw("COLZ");
	}
}

