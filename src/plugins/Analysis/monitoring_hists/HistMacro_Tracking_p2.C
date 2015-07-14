// hnamepath: /Independent/Hist_Reconstruction/Tracking/q+/PVsTheta_Candidates_q+
// hnamepath: /Independent/Hist_Reconstruction/Tracking/q-/PVsTheta_Candidates_q-
// hnamepath: /Independent/Hist_Reconstruction/Tracking/q+/PVsTheta_WireBased_q+
// hnamepath: /Independent/Hist_Reconstruction/Tracking/q-/PVsTheta_WireBased_q-
// hnamepath: /Independent/Hist_Reconstruction/Tracking/q+/PVsTheta_TimeBased_q+
// hnamepath: /Independent/Hist_Reconstruction/Tracking/q-/PVsTheta_TimeBased_q-

{
	//Goto Path
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_Reconstruction");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Tracking/q+");
	TH2I* locHist_PVsTheta_Candidates_QPlus = (TH2I*)gDirectory->Get("PVsTheta_Candidates_q+");
	TH2I* locHist_PVsTheta_WireBased_QPlus = (TH2I*)gDirectory->Get("PVsTheta_WireBased_q+");
	TH2I* locHist_PVsTheta_TimeBased_QPlus = (TH2I*)gDirectory->Get("PVsTheta_TimeBased_q+");

	gDirectory->cd("../q-");
	TH2I* locHist_PVsTheta_Candidates_QMinus = (TH2I*)gDirectory->Get("PVsTheta_Candidates_q-");
	TH2I* locHist_PVsTheta_WireBased_QMinus = (TH2I*)gDirectory->Get("PVsTheta_WireBased_q-");
	TH2I* locHist_PVsTheta_TimeBased_QMinus = (TH2I*)gDirectory->Get("PVsTheta_TimeBased_q-");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Tracking_p2", "Tracking_p2", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_Candidates_QPlus != NULL)
	{
		locHist_PVsTheta_Candidates_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Candidates_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Candidates_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Candidates_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Candidates_QPlus->Draw("COLZ");
		gPad->SetLogz();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_WireBased_QPlus != NULL)
	{
		locHist_PVsTheta_WireBased_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_WireBased_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_WireBased_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_WireBased_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_WireBased_QPlus->Draw("COLZ");
		gPad->SetLogz();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_TimeBased_QPlus != NULL)
	{
		locHist_PVsTheta_TimeBased_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_TimeBased_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_TimeBased_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_TimeBased_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_TimeBased_QPlus->Draw("COLZ");
		gPad->SetLogz();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_Candidates_QMinus != NULL)
	{
		locHist_PVsTheta_Candidates_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Candidates_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_Candidates_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Candidates_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_Candidates_QMinus->Draw("COLZ");
		gPad->SetLogz();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_WireBased_QMinus != NULL)
	{
		locHist_PVsTheta_WireBased_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_WireBased_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_WireBased_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_WireBased_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_WireBased_QMinus->Draw("COLZ");
		gPad->SetLogz();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PVsTheta_TimeBased_QMinus != NULL)
	{
		locHist_PVsTheta_TimeBased_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_TimeBased_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_PVsTheta_TimeBased_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_TimeBased_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_PVsTheta_TimeBased_QMinus->Draw("COLZ");
		gPad->SetLogz();
	}
}

