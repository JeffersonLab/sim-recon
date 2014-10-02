// hnamepath: /Independent/Hist_DetectedParticleKinematics/PID/BetaVsP_Q+
// hnamepath: /Independent/Hist_DetectedParticleKinematics/PID/BetaVsP_Q-
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Pi+/PVsTheta
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Pi+/PhiVsTheta
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Pi-/PVsTheta
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Pi-/PhiVsTheta
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Proton/PVsTheta
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Proton/PhiVsTheta
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Photon/PVsTheta
// hnamepath: /Independent/Hist_DetectedParticleKinematics/Photon/PhiVsTheta

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectedParticleKinematics");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("PID");
	TH2I* locHist_BetaVsPQPlus = (TH2I*)gDirectory->Get("BetaVsP_Q+");
	TH2I* locHist_BetaVsPQMinus = (TH2I*)gDirectory->Get("BetaVsP_Q-");
	gDirectory->cd("../Pi+");
	TH2I* locHist_PVsTheta_PiPlus = (TH2I*)gDirectory->Get("PVsTheta");
	TH2I* locHist_PhiVsTheta_PiPlus = (TH2I*)gDirectory->Get("PhiVsTheta");
	gDirectory->cd("../Pi-");
	TH2I* locHist_PVsTheta_PiMinus = (TH2I*)gDirectory->Get("PVsTheta");
	TH2I* locHist_PhiVsTheta_PiMinus = (TH2I*)gDirectory->Get("PhiVsTheta");
	gDirectory->cd("../Proton");
	TH2I* locHist_PVsTheta_Proton = (TH2I*)gDirectory->Get("PVsTheta");
	TH2I* locHist_PhiVsTheta_Proton = (TH2I*)gDirectory->Get("PhiVsTheta");
	gDirectory->cd("../Photon");
	TH2I* locHist_PVsTheta_Photon = (TH2I*)gDirectory->Get("PVsTheta");
	TH2I* locHist_PhiVsTheta_Photon = (TH2I*)gDirectory->Get("PhiVsTheta");


	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Kinematics"); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(4, 3);

	//Draw
	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsPQPlus != NULL)
	{
		locHist_BetaVsPQPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQPlus->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsPQMinus != NULL)
	{
		locHist_BetaVsPQMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQMinus->Draw("COLZ");
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_PiPlus != NULL)
	{
		locHist_PVsTheta_PiPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_PiPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_PiPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_PiPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_PiPlus->Draw("COLZ");
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_PiMinus != NULL)
	{
		locHist_PVsTheta_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_PiMinus->Draw("COLZ");
	}

	locCanvas->cd(7);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_Proton != NULL)
	{
		locHist_PVsTheta_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Proton->Draw("COLZ");
	}

	locCanvas->cd(8);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_Photon != NULL)
	{
		locHist_PVsTheta_Photon->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Photon->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Photon->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Photon->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Photon->Draw("COLZ");
	}

	locCanvas->cd(9);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PhiVsTheta_PiPlus != NULL)
	{
		locHist_PhiVsTheta_PiPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_PiPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_PiPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_PiPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_PiPlus->Draw("COLZ");
	}

	locCanvas->cd(10);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PhiVsTheta_PiMinus != NULL)
	{
		locHist_PhiVsTheta_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_PiMinus->Draw("COLZ");
	}

	locCanvas->cd(11);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PhiVsTheta_Proton != NULL)
	{
		locHist_PhiVsTheta_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_Proton->Draw("COLZ");
	}

	locCanvas->cd(12);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PhiVsTheta_Photon != NULL)
	{
		locHist_PhiVsTheta_Photon->GetXaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_Photon->GetYaxis()->SetTitleSize(0.05);
		locHist_PhiVsTheta_Photon->GetXaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_Photon->GetYaxis()->SetLabelSize(0.05);
		locHist_PhiVsTheta_Photon->Draw("COLZ");
	}
}

