// hnamepath: /Independent/Hist_Reconstruction/BCALShowerPhiVsZ
// hnamepath: /Independent/Hist_Reconstruction/BCALShowerEnergy
// hnamepath: /Independent/Hist_Neutrals/BCALNeutralShowerEnergy
// hnamepath: /Independent/Hist_Neutrals/BCALNeutralShowerDeltaT
// hnamepath: /Independent/Hist_Neutrals/BCALNeutralShowerDeltaTVsE

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_Reconstruction");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2I* locHist_BCALShowerPhiVsZ = (TH2I*)gDirectory->Get("BCALShowerPhiVsZ");
	TH1I* locHist_BCALShowerEnergy = (TH1I*)gDirectory->Get("BCALShowerEnergy");

	gDirectory->cd("../Hist_Neutrals");
	TH1I* locHist_BCALNeutralShowerEnergy = (TH1I*)gDirectory->Get("BCALNeutralShowerEnergy"); //photon(?)
	TH1I* locHist_BCALNeutralShowerDeltaT = (TH1I*)gDirectory->Get("BCALNeutralShowerDeltaT"); //photon(?)
	TH2I* locHist_BCALNeutralShowerDeltaTVsE = (TH2I*)gDirectory->Get("BCALNeutralShowerDeltaTVsE"); //photon(?)

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("BCALReconstruction_p1", "BCALReconstruction_p1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerPhiVsZ != NULL)
	{
		locHist_BCALShowerPhiVsZ->GetXaxis()->SetTitleSize(0.04);
		locHist_BCALShowerPhiVsZ->GetYaxis()->SetTitleSize(0.04);
		locHist_BCALShowerPhiVsZ->GetXaxis()->SetLabelSize(0.04);
		locHist_BCALShowerPhiVsZ->GetYaxis()->SetLabelSize(0.04);
		locHist_BCALShowerPhiVsZ->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALShowerEnergy != NULL)
	{
		locHist_BCALShowerEnergy->GetXaxis()->SetTitleSize(0.04);
		locHist_BCALShowerEnergy->GetXaxis()->SetLabelSize(0.04);
		locHist_BCALShowerEnergy->GetYaxis()->SetLabelSize(0.04);
		locHist_BCALShowerEnergy->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALNeutralShowerEnergy != NULL)
	{
		locHist_BCALNeutralShowerEnergy->GetXaxis()->SetTitleSize(0.04);
		locHist_BCALNeutralShowerEnergy->GetXaxis()->SetLabelSize(0.04);
		locHist_BCALNeutralShowerEnergy->GetYaxis()->SetLabelSize(0.04);
		locHist_BCALNeutralShowerEnergy->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALNeutralShowerDeltaT != NULL)
	{
		locHist_BCALNeutralShowerDeltaT->GetXaxis()->SetTitleSize(0.04);
		locHist_BCALNeutralShowerDeltaT->GetXaxis()->SetLabelSize(0.04);
		locHist_BCALNeutralShowerDeltaT->GetYaxis()->SetLabelSize(0.04);
		locHist_BCALNeutralShowerDeltaT->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCALNeutralShowerDeltaTVsE != NULL)
	{
		locHist_BCALNeutralShowerDeltaTVsE->GetXaxis()->SetTitleSize(0.04);
		locHist_BCALNeutralShowerDeltaTVsE->GetYaxis()->SetTitleSize(0.04);
		locHist_BCALNeutralShowerDeltaTVsE->GetXaxis()->SetLabelSize(0.04);
		locHist_BCALNeutralShowerDeltaTVsE->GetYaxis()->SetLabelSize(0.04);
		locHist_BCALNeutralShowerDeltaTVsE->Draw("COLZ");
	}
}

