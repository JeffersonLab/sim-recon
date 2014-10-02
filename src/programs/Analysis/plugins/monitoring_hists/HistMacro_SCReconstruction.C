// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/SCHitEnergy
// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/SCHitEnergyVsSector
// hnamepath: /Independent/Hist_DetectorStudies/dEdxVsPByCharge/q-/SCdEdXVsP
// hnamepath: /Independent/Hist_DetectorStudies/dEdxVsPByCharge/q+/SCdEdXVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaT
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaTVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaTVsPhi
// hnamepath: /Independent/Hist_DetectorStudies/Matched_ReconstructedPID/Proton/SCDeltaTVsTheta
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
	gDirectory->cd("Reconstruction");
	TH1I* locHist_SCHitEnergy = (TH1I*)gDirectory->Get("SCHitEnergy");
	TH2I* locHist_SCHitEnergyVsSector = (TH2I*)gDirectory->Get("SCHitEnergyVsSector");

	gDirectory->cd("../dEdxVsPByCharge/q-");
	TH2I* locHist_SCdEdXVsP_QMinus = (TH2I*)gDirectory->Get("SCdEdXVsP"); //q-
	gDirectory->cd("../q+");
	TH2I* locHist_SCdEdXVsP_QPlus = (TH2I*)gDirectory->Get("SCdEdXVsP"); //q-

	gDirectory->cd("../../Matched_ReconstructedPID/Proton");
	TH1I* locHist_SCDeltaT_Proton = (TH1I*)gDirectory->Get("SCDeltaT"); //proton
	TH2I* locHist_SCDeltaTVsP_Proton = (TH2I*)gDirectory->Get("SCDeltaTVsP"); //proton
	TH2I* locHist_SCDeltaTVsPhi_Proton = (TH2I*)gDirectory->Get("SCDeltaTVsPhi"); //proton
	TH2I* locHist_SCDeltaTVsTheta_Proton = (TH2I*)gDirectory->Get("SCDeltaTVsTheta"); //proton

	gDirectory->cd("../Pi-");
	TH1I* locHist_SCDeltaT_PiMinus = (TH1I*)gDirectory->Get("SCDeltaT"); //pi-
	TH2I* locHist_SCDeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("SCDeltaTVsP"); //pi-
	TH2I* locHist_SCDeltaTVsPhi_PiMinus = (TH2I*)gDirectory->Get("SCDeltaTVsPhi"); //pi-
	TH2I* locHist_SCDeltaTVsTheta_PiMinus = (TH2I*)gDirectory->Get("SCDeltaTVsTheta"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("SCReconstruction"); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(4, 3);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCHitEnergy != NULL)
	{
		locHist_SCHitEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_SCHitEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_SCHitEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_SCHitEnergy->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCHitEnergyVsSector != NULL)
	{
		locHist_SCHitEnergyVsSector->GetXaxis()->SetTitleSize(0.05);
		locHist_SCHitEnergyVsSector->GetYaxis()->SetTitleSize(0.05);
		locHist_SCHitEnergyVsSector->GetXaxis()->SetLabelSize(0.05);
		locHist_SCHitEnergyVsSector->GetYaxis()->SetLabelSize(0.05);
		locHist_SCHitEnergyVsSector->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCdEdXVsP_QPlus != NULL)
	{
		locHist_SCdEdXVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_SCdEdXVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_SCdEdXVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_SCdEdXVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_SCdEdXVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCdEdXVsP_QMinus != NULL)
	{
		locHist_SCdEdXVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_SCdEdXVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_SCdEdXVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_SCdEdXVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_SCdEdXVsP_QMinus->Draw("COLZ");
	}

	//Proton
	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaT_Proton != NULL)
	{
		locHist_SCDeltaT_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaT_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_Proton->Draw();
	}

	locCanvas->cd(6);
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

	locCanvas->cd(7);
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

	locCanvas->cd(8);
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

	//PiMinus
	locCanvas->cd(9);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCDeltaT_PiMinus != NULL)
	{
		locHist_SCDeltaT_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_SCDeltaT_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_SCDeltaT_PiMinus->Draw();
	}

	locCanvas->cd(10);
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


	locCanvas->cd(11);
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

	locCanvas->cd(12);
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

