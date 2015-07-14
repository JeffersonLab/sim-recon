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

	//Beta-vs-p functions
	TF1* locBetaVsPFunc_Proton = new TF1("BetaVsPFunc_Proton", "x/sqrt(x*x + [0]*[0])", 0.0, 12.0);
	locBetaVsPFunc_Proton->SetParameter(0, 0.938272046);
	locBetaVsPFunc_Proton->SetLineWidth(2);
	locBetaVsPFunc_Proton->SetLineColor(kBlack);
	locBetaVsPFunc_Proton->SetNpx(1000);

	TF1* locBetaVsPFunc_Kaon = new TF1("BetaVsPFunc_Kaon", "x/sqrt(x*x + [0]*[0])", 0.0, 12.0);
	locBetaVsPFunc_Kaon->SetParameter(0, 0.493677);
	locBetaVsPFunc_Kaon->SetLineWidth(2);
	locBetaVsPFunc_Kaon->SetLineColor(kBlack);
	locBetaVsPFunc_Kaon->SetNpx(1000);

	TF1* locBetaVsPFunc_Pion = new TF1("BetaVsPFunc_Pion", "x/sqrt(x*x + [0]*[0])", 0.0, 12.0);
	locBetaVsPFunc_Pion->SetParameter(0, 0.13957018);
	locBetaVsPFunc_Pion->SetLineWidth(2);
	locBetaVsPFunc_Pion->SetLineColor(kBlack);
	locBetaVsPFunc_Pion->SetNpx(1000);

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
		locHist_BetaVsPQPlus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_BetaVsPQPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQPlus->GetListOfFunctions()->Add(locBetaVsPFunc_Proton);
		locHist_BetaVsPQPlus->GetListOfFunctions()->Add(locBetaVsPFunc_Kaon);
		locHist_BetaVsPQPlus->GetListOfFunctions()->Add(locBetaVsPFunc_Pion);
		locHist_BetaVsPQPlus->Draw("COLZ");
                gPad->SetLogz();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsPQMinus != NULL)
	{
		locHist_BetaVsPQMinus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_BetaVsPQMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsPQMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsPQMinus->GetListOfFunctions()->Add(locBetaVsPFunc_Pion);
		locHist_BetaVsPQMinus->GetListOfFunctions()->Add(locBetaVsPFunc_Kaon);
		locHist_BetaVsPQMinus->Draw("COLZ");
                gPad->SetLogz();
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

