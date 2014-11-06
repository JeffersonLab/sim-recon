// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/FCALShowerYVsX
// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/FCALShowerEnergy
// hnamepath: /Independent/Hist_DetectorStudies/Not-Matched/FCALNeutralShowerEnergy
// hnamepath: /Independent/Hist_DetectorStudies/Not-Matched/FCALNeutralShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Not-Matched/FCALNeutralShowerDeltaTVsE
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
	gDirectory->cd("Reconstruction");
	TH2I* locHist_FCALShowerYVsX = (TH2I*)gDirectory->Get("FCALShowerYVsX");
	TH1I* locHist_FCALShowerEnergy = (TH1I*)gDirectory->Get("FCALShowerEnergy");

	gDirectory->cd("../Not-Matched");
	TH1I* locHist_FCALNeutralShowerEnergy = (TH1I*)gDirectory->Get("FCALNeutralShowerEnergy"); //photon(?)
	TH1I* locHist_FCALNeutralShowerDeltaT = (TH1I*)gDirectory->Get("FCALNeutralShowerDeltaT"); //photon(?)
	TH2I* locHist_FCALNeutralShowerDeltaTVsE = (TH2I*)gDirectory->Get("FCALNeutralShowerDeltaTVsE"); //photon(?)

	gDirectory->cd("../Matched_ReconstructedPID/Proton");
	TH1I* locHist_FCALShowerDeltaT_Proton = (TH1I*)gDirectory->Get("FCALShowerDeltaT"); //proton
	TH2I* locHist_FCALShowerDeltaTVsP_Proton = (TH2I*)gDirectory->Get("FCALShowerDeltaTVsP"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_FCALShowerDeltaT_PiMinus = (TH1I*)gDirectory->Get("FCALShowerDeltaT"); //pi-
	TH2I* locHist_FCALShowerDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("FCALShowerDeltaTVsP"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("FCALReconstruction", "FCALReconstruction", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 3);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerYVsX != NULL)
	{
		locHist_FCALShowerYVsX->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerYVsX->GetYaxis()->SetTitleSize(0.05);
		locHist_FCALShowerYVsX->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerYVsX->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerYVsX->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerEnergy != NULL)
	{
		locHist_FCALShowerEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerEnergy->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALNeutralShowerEnergy != NULL)
	{
		locHist_FCALNeutralShowerEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALNeutralShowerEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALNeutralShowerEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALNeutralShowerEnergy->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALNeutralShowerDeltaT != NULL)
	{
		locHist_FCALNeutralShowerDeltaT->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALNeutralShowerDeltaT->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALNeutralShowerDeltaT->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALNeutralShowerDeltaT->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaT_Proton != NULL)
	{
		locHist_FCALShowerDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_Proton->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerDeltaT_PiMinus != NULL)
	{
		locHist_FCALShowerDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALShowerDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(7);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALNeutralShowerDeltaTVsE != NULL)
	{
		locHist_FCALNeutralShowerDeltaTVsE->GetXaxis()->SetTitleSize(0.05);
		locHist_FCALNeutralShowerDeltaTVsE->GetYaxis()->SetTitleSize(0.05);
		locHist_FCALNeutralShowerDeltaTVsE->GetXaxis()->SetLabelSize(0.05);
		locHist_FCALNeutralShowerDeltaTVsE->GetYaxis()->SetLabelSize(0.05);
		locHist_FCALNeutralShowerDeltaTVsE->Draw("COLZ");
	}

	locCanvas->cd(8);
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

	locCanvas->cd(9);
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

