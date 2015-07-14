// hnamepath: /RF/DeltaT_RF_TAGH/FDCRF_TaggerDeltaT
// hnamepath: /RF/DeltaT_RF_TAGH/TOFRF_TaggerDeltaT
// hnamepath: /RF/DeltaT_RF_TAGH/TAGHRF_TaggerDeltaT
// hnamepath: /RF/DeltaT_RF_TAGH/PSCRF_TaggerDeltaT

{
	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("DeltaT_RF_TAGH");
	TH1I* locHist_FDCRF_TaggerDeltaT = (TH1I*)gDirectory->Get("FDCRF_TaggerDeltaT");
	TH1I* locHist_TOFRF_TaggerDeltaT = (TH1I*)gDirectory->Get("TOFRF_TaggerDeltaT");
	TH1I* locHist_TAGHRF_TaggerDeltaT = (TH1I*)gDirectory->Get("TAGHRF_TaggerDeltaT");
	TH1I* locHist_PSCRF_TaggerDeltaT = (TH1I*)gDirectory->Get("PSCRF_TaggerDeltaT");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_p3", "RF_p3", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FDCRF_TaggerDeltaT != NULL)
	{
		locHist_FDCRF_TaggerDeltaT->GetXaxis()->SetTitleSize(0.05);
		locHist_FDCRF_TaggerDeltaT->GetYaxis()->SetTitleSize(0.05);
		locHist_FDCRF_TaggerDeltaT->GetXaxis()->SetLabelSize(0.05);
		locHist_FDCRF_TaggerDeltaT->GetYaxis()->SetLabelSize(0.05);
		locHist_FDCRF_TaggerDeltaT->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFRF_TaggerDeltaT != NULL)
	{
		locHist_TOFRF_TaggerDeltaT->GetXaxis()->SetTitleSize(0.05);
		locHist_TOFRF_TaggerDeltaT->GetYaxis()->SetTitleSize(0.05);
		locHist_TOFRF_TaggerDeltaT->GetXaxis()->SetLabelSize(0.05);
		locHist_TOFRF_TaggerDeltaT->GetYaxis()->SetLabelSize(0.05);
		locHist_TOFRF_TaggerDeltaT->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TAGHRF_TaggerDeltaT != NULL)
	{
		locHist_TAGHRF_TaggerDeltaT->GetXaxis()->SetTitleSize(0.05);
		locHist_TAGHRF_TaggerDeltaT->GetYaxis()->SetTitleSize(0.05);
		locHist_TAGHRF_TaggerDeltaT->GetXaxis()->SetLabelSize(0.05);
		locHist_TAGHRF_TaggerDeltaT->GetYaxis()->SetLabelSize(0.05);
		locHist_TAGHRF_TaggerDeltaT->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PSCRF_TaggerDeltaT != NULL)
	{
		locHist_PSCRF_TaggerDeltaT->GetXaxis()->SetTitleSize(0.05);
		locHist_PSCRF_TaggerDeltaT->GetYaxis()->SetTitleSize(0.05);
		locHist_PSCRF_TaggerDeltaT->GetXaxis()->SetLabelSize(0.05);
		locHist_PSCRF_TaggerDeltaT->GetYaxis()->SetLabelSize(0.05);
		locHist_PSCRF_TaggerDeltaT->Draw();
	}
}
