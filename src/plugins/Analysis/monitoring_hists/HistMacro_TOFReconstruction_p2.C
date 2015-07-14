// hnamepath: /Independent/Hist_DetectorPID/TOF/DeltaTVsP_Proton
// hnamepath: /Independent/Hist_DetectorPID/TOF/DeltaTVsP_Pi+
// hnamepath: /Independent/Hist_DetectorPID/TOF/DeltaTVsP_Pi-

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorPID");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("TOF");
	TH2I* locHist_DeltaTVsP_Proton = (TH2I*)gDirectory->Get("DeltaTVsP_Proton"); //proton
	TH2I* locHist_DeltaTVsP_PiPlus = (TH2I*)gDirectory->Get("DeltaTVsP_Pi+"); //pi+
	TH2I* locHist_DeltaTVsP_PiMinus = (TH2I*)gDirectory->Get("DeltaTVsP_Pi-"); //pi-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TOFReconstruction_p2", "TOFReconstruction_p2", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_DeltaTVsP_PiMinus != NULL)
	{
		locHist_DeltaTVsP_PiMinus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_DeltaTVsP_PiMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_DeltaTVsP_PiMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_DeltaTVsP_PiMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_DeltaTVsP_PiMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_DeltaTVsP_PiMinus->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_DeltaTVsP_PiMinus != NULL)
	{
		TH1I* locHist = (TH1I*)locHist_DeltaTVsP_PiMinus->ProjectionY("DeltaTVsP_PiMinus_1D");
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_DeltaTVsP_PiPlus != NULL)
	{
		locHist_DeltaTVsP_PiPlus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_DeltaTVsP_PiPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_DeltaTVsP_PiPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_DeltaTVsP_PiPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_DeltaTVsP_PiPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_DeltaTVsP_PiPlus->Draw("COLZ");
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_DeltaTVsP_PiPlus != NULL)
	{
		TH1I* locHist = (TH1I*)locHist_DeltaTVsP_PiPlus->ProjectionY("DeltaTVsP_PiPlus_1D");
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_DeltaTVsP_Proton != NULL)
	{
		locHist_DeltaTVsP_Proton->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_DeltaTVsP_Proton->GetXaxis()->SetTitleSize(0.05);
		locHist_DeltaTVsP_Proton->GetXaxis()->SetLabelSize(0.05);
		locHist_DeltaTVsP_Proton->GetYaxis()->SetTitleSize(0.05);
		locHist_DeltaTVsP_Proton->GetYaxis()->SetLabelSize(0.05);
		locHist_DeltaTVsP_Proton->Draw("COLZ");
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_DeltaTVsP_Proton != NULL)
	{
		TH1I* locHist = (TH1I*)locHist_DeltaTVsP_Proton->ProjectionY("DeltaTVsP_Proton_1D");
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}
}

