// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/BCALShowerPhiVsZ
// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/BCALShowerEnergy
// hnamepath: /Independent/Hist_DetectorStudies/Not-Matched/BCALNeutralShowerEnergy
// hnamepath: /Independent/Hist_DetectorStudies/Not-Matched/BCALNeutralShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Not-Matched/BCALNeutralShowerDeltaTVsE
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/BCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/BCALShowerDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/BCALShowerDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/PiMinus/BCALShowerDeltaTVsP

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Reconstruction");
	TH2I* locHist_BCALShowerPhiVsZ = (TH2I*)gDirectory->Get("BCALShowerPhiVsZ");
	TH1I* locHist_BCALShowerEnergy = (TH1I*)gDirectory->Get("BCALShowerEnergy");

	gDirectory->cd("../Not-Matched");
	TH1I* locHist_BCALNeutralShowerEnergy = (TH1I*)gDirectory->Get("BCALNeutralShowerEnergy"); //photon(?)
	TH1I* locHist_BCALNeutralShowerDeltaT = (TH1I*)gDirectory->Get("BCALNeutralShowerDeltaT"); //photon(?)
	TH2I* locHist_BCALNeutralShowerDeltaTVsE = (TH2I*)gDirectory->Get("BCALNeutralShowerDeltaTVsE"); //photon(?)

	gDirectory->cd("../Matched_ReconstructedPID/Proton");
	TH1I* locHist_BCALShowerDeltaT_Proton = (TH1I*)gDirectory->Get("BCALShowerDeltaT"); //proton
	TH2I* locHist_BCALShowerDeltaTVsP_Proton = (TH2I*)gDirectory->Get("BCALShowerDeltaTVsP"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_BCALShowerDeltaT_PiMinus = (TH1I*)gDirectory->Get("BCALShowerDeltaT"); //pi-
	TH2I* locHist_BCALShowerDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("BCALShowerDeltaTVsP"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("BCALReconstruction"); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 3);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerPhiVsZ != NULL)
	{
		locHist_BCALShowerPhiVsZ->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerPhiVsZ->GetYaxis()->SetTitleSize(0.05);
		locHist_BCALShowerPhiVsZ->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerPhiVsZ->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerPhiVsZ->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerEnergy != NULL)
	{
		locHist_BCALShowerEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerEnergy->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALNeutralShowerEnergy != NULL)
	{
		locHist_BCALNeutralShowerEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALNeutralShowerEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALNeutralShowerEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALNeutralShowerEnergy->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALNeutralShowerDeltaT != NULL)
	{
		locHist_BCALNeutralShowerDeltaT->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALNeutralShowerDeltaT->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALNeutralShowerDeltaT->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALNeutralShowerDeltaT->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerDeltaT_Proton != NULL)
	{
		locHist_BCALShowerDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_Proton->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerDeltaT_PiMinus != NULL)
	{
		locHist_BCALShowerDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALShowerDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALShowerDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(7);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALNeutralShowerDeltaTVsE != NULL)
	{
		locHist_BCALNeutralShowerDeltaTVsE->GetXaxis()->SetTitleSize(0.05);
		locHist_BCALNeutralShowerDeltaTVsE->GetYaxis()->SetTitleSize(0.05);
		locHist_BCALNeutralShowerDeltaTVsE->GetXaxis()->SetLabelSize(0.05);
		locHist_BCALNeutralShowerDeltaTVsE->GetYaxis()->SetLabelSize(0.05);
		locHist_BCALNeutralShowerDeltaTVsE->Draw("COLZ");
	}

	locCanvas->cd(8);
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

	locCanvas->cd(9);
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
}

