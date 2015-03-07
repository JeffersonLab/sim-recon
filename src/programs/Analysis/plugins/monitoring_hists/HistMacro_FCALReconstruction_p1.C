// hnamepath: /Independent/Hist_Reconstruction/FCALShowerYVsX
// hnamepath: /Independent/Hist_Reconstruction/FCALShowerEnergy
// hnamepath: /Independent/Hist_Neutrals/FCALNeutralShowerEnergy
// hnamepath: /Independent/Hist_Neutrals/FCALNeutralShowerDeltaT
// hnamepath: /Independent/Hist_Neutrals/FCALNeutralShowerDeltaTVsE

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_Reconstruction");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2I* locHist_FCALShowerYVsX = (TH2I*)gDirectory->Get("FCALShowerYVsX");
	TH1I* locHist_FCALShowerEnergy = (TH1I*)gDirectory->Get("FCALShowerEnergy");

	gDirectory->cd("../Hist_Neutrals");
	TH1I* locHist_FCALNeutralShowerEnergy = (TH1I*)gDirectory->Get("FCALNeutralShowerEnergy"); //photon(?)
	TH1I* locHist_FCALNeutralShowerDeltaT = (TH1I*)gDirectory->Get("FCALNeutralShowerDeltaT"); //photon(?)
	TH2I* locHist_FCALNeutralShowerDeltaTVsE = (TH2I*)gDirectory->Get("FCALNeutralShowerDeltaTVsE"); //photon(?)

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("FCALReconstruction_p1", "FCALReconstruction_p1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerYVsX != NULL)
	{
		locHist_FCALShowerYVsX->GetXaxis()->SetTitleSize(0.04);
		locHist_FCALShowerYVsX->GetYaxis()->SetTitleSize(0.04);
		locHist_FCALShowerYVsX->GetXaxis()->SetLabelSize(0.04);
		locHist_FCALShowerYVsX->GetYaxis()->SetLabelSize(0.04);
		locHist_FCALShowerYVsX->Draw("COLZ");
		gPad->SetLogz();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALShowerEnergy != NULL)
	{
		locHist_FCALShowerEnergy->GetXaxis()->SetTitleSize(0.04);
		locHist_FCALShowerEnergy->GetXaxis()->SetLabelSize(0.04);
		locHist_FCALShowerEnergy->GetYaxis()->SetLabelSize(0.04);
		locHist_FCALShowerEnergy->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALNeutralShowerEnergy != NULL)
	{
		locHist_FCALNeutralShowerEnergy->GetXaxis()->SetTitleSize(0.04);
		locHist_FCALNeutralShowerEnergy->GetXaxis()->SetLabelSize(0.04);
		locHist_FCALNeutralShowerEnergy->GetYaxis()->SetLabelSize(0.04);
		locHist_FCALNeutralShowerEnergy->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALNeutralShowerDeltaT != NULL)
	{
		locHist_FCALNeutralShowerDeltaT->GetXaxis()->SetTitleSize(0.04);
		locHist_FCALNeutralShowerDeltaT->GetXaxis()->SetLabelSize(0.04);
		locHist_FCALNeutralShowerDeltaT->GetYaxis()->SetLabelSize(0.04);
		locHist_FCALNeutralShowerDeltaT->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCALNeutralShowerDeltaTVsE != NULL)
	{
		locHist_FCALNeutralShowerDeltaTVsE->GetXaxis()->SetTitleSize(0.04);
		locHist_FCALNeutralShowerDeltaTVsE->GetYaxis()->SetTitleSize(0.04);
		locHist_FCALNeutralShowerDeltaTVsE->GetXaxis()->SetLabelSize(0.04);
		locHist_FCALNeutralShowerDeltaTVsE->GetYaxis()->SetLabelSize(0.04);
		locHist_FCALNeutralShowerDeltaTVsE->Draw("COLZ");
	}
}

