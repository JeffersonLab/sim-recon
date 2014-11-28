// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/SCHitEnergy
// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/SCHitEnergyVsSector
// hnamepath: /Independent/Hist_DetectorStudies/PID/q-/SCdEdXVsP
// hnamepath: /Independent/Hist_DetectorStudies/PID/q+/SCdEdXVsP
// hnamepath: /Independent/Hist_DetectorStudies/PID/q-/BetaVsP_SC
// hnamepath: /Independent/Hist_DetectorStudies/PID/q+/BetaVsP_SC

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

	gDirectory->cd("../PID/q-");
	TH2I* locHist_SCdEdXVsP_QMinus = (TH2I*)gDirectory->Get("SCdEdXVsP"); //q-
	TH2I* locHist_BetaVsP_QMinus = (TH2I*)gDirectory->Get("BetaVsP_SC"); //q-
	gDirectory->cd("../q+");
	TH2I* locHist_SCdEdXVsP_QPlus = (TH2I*)gDirectory->Get("SCdEdXVsP"); //q+
	TH2I* locHist_BetaVsP_QPlus = (TH2I*)gDirectory->Get("BetaVsP_SC"); //q+

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("SCReconstruction_p1", "SCReconstruction_p1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SCHitEnergy != NULL)
	{
		locHist_SCHitEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_SCHitEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_SCHitEnergy->Draw();
	}

	locCanvas->cd(4);
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

	locCanvas->cd(2);
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

	locCanvas->cd(5);
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

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QPlus != NULL)
	{
		locHist_BetaVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BetaVsP_QMinus != NULL)
	{
		locHist_BetaVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_BetaVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_BetaVsP_QMinus->Draw("COLZ");
	}
}

