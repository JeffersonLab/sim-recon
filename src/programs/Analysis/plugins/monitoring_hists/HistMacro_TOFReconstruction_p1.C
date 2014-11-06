// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/TOFPointYVsX
// hnamepath: /Independent/Hist_DetectorStudies/Reconstruction/TOFPointEnergy
// hnamepath: /Independent/Hist_DetectorStudies/dEdxVsPByCharge/q-/TOFdEdXVsP
// hnamepath: /Independent/Hist_DetectorStudies/dEdxVsPByCharge/q+/TOFdEdXVsP

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Reconstruction");
	TH2I* locHist_TOFPointYVsX = (TH2I*)gDirectory->Get("TOFPointYVsX");
	TH1I* locHist_TOFPointEnergy = (TH1I*)gDirectory->Get("TOFPointEnergy");

	gDirectory->cd("../dEdxVsPByCharge/q-");
	TH2I* locHist_TOFdEdXVsP_QMinus = (TH2I*)gDirectory->Get("TOFdEdXVsP"); //q-
	gDirectory->cd("../q+");
	TH2I* locHist_TOFdEdXVsP_QPlus = (TH2I*)gDirectory->Get("TOFdEdXVsP"); //q-

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TOFReconstruction_p1", "TOFReconstruction_p1", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFPointYVsX != NULL)
	{
		locHist_TOFPointYVsX->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFPointYVsX->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFPointYVsX->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFPointYVsX->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFPointYVsX->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFPointEnergy != NULL)
	{
		locHist_TOFPointEnergy->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFPointEnergy->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFPointEnergy->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFPointEnergy->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFdEdXVsP_QPlus != NULL)
	{
		locHist_TOFdEdXVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFdEdXVsP_QMinus != NULL)
	{
		locHist_TOFdEdXVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFdEdXVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFdEdXVsP_QMinus->Draw("COLZ");
	}
}

