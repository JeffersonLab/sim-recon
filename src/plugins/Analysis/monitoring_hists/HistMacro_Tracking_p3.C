// hnamepath: /Independent/Hist_DetectorPID/CDC/dEdxVsP_q-
// hnamepath: /Independent/Hist_DetectorPID/CDC/dEdxVsP_q+
// hnamepath: /Independent/Hist_DetectorPID/FDC/dEdxVsP_q-
// hnamepath: /Independent/Hist_DetectorPID/FDC/dEdxVsP_q+

{
	//Goto Path
	TDirectory *locInitDirectory = gDirectory;
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorPID");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("CDC");
	TH2I* locHist_CDCdEdXVsP_QMinus = (TH2I*)gDirectory->Get("dEdXVsP_q-"); //q-
	TH2I* locHist_CDCdEdXVsP_QPlus = (TH2I*)gDirectory->Get("dEdXVsP_q+"); //q+
	gDirectory->cd("../FDC");
	TH2I* locHist_FDCdEdXVsP_QMinus = (TH2I*)gDirectory->Get("dEdXVsP_q-"); //q-
	TH2I* locHist_FDCdEdXVsP_QPlus = (TH2I*)gDirectory->Get("dEdXVsP_q+"); //q+

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Tracking_p3", "Tracking_p3", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_CDCdEdXVsP_QPlus != NULL)
	{
		locHist_CDCdEdXVsP_QPlus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_CDCdEdXVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_CDCdEdXVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_CDCdEdXVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_CDCdEdXVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_CDCdEdXVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_CDCdEdXVsP_QMinus != NULL)
	{
		locHist_CDCdEdXVsP_QMinus->GetXaxis()->SetRangeUser(0.0, 6.0);
		locHist_CDCdEdXVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_CDCdEdXVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_CDCdEdXVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_CDCdEdXVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_CDCdEdXVsP_QMinus->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FDCdEdXVsP_QPlus != NULL)
	{
		locHist_FDCdEdXVsP_QPlus->GetXaxis()->SetTitleSize(0.05);
		locHist_FDCdEdXVsP_QPlus->GetYaxis()->SetTitleSize(0.05);
		locHist_FDCdEdXVsP_QPlus->GetXaxis()->SetLabelSize(0.05);
		locHist_FDCdEdXVsP_QPlus->GetYaxis()->SetLabelSize(0.05);
		locHist_FDCdEdXVsP_QPlus->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FDCdEdXVsP_QMinus != NULL)
	{
		locHist_FDCdEdXVsP_QMinus->GetXaxis()->SetTitleSize(0.05);
		locHist_FDCdEdXVsP_QMinus->GetYaxis()->SetTitleSize(0.05);
		locHist_FDCdEdXVsP_QMinus->GetXaxis()->SetLabelSize(0.05);
		locHist_FDCdEdXVsP_QMinus->GetYaxis()->SetLabelSize(0.05);
		locHist_FDCdEdXVsP_QMinus->Draw("COLZ");
	}
}

