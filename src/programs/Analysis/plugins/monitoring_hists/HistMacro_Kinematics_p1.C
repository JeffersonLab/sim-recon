// hnamepath: /Independent/Hist_DetectedParticleKinematics/PID/BetaVsP_Q+
// hnamepath: /Independent/Hist_DetectedParticleKinematics/PID/BetaVsP_Q-
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
	gDirectory->cd("../Photon");
	TH2I* locHist_PVsTheta_Photon = (TH2I*)gDirectory->Get("PVsTheta");
	TH2I* locHist_PhiVsTheta_Photon = (TH2I*)gDirectory->Get("PhiVsTheta");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Kinematics_p1", "Kinematics_p1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
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

	locCanvas->cd(2);
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

	locCanvas->cd(3);
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

	locCanvas->cd(4);
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

