// hnamepath: /RF/DeltaT_RF_FirstTime/FDCRF_FirstTimeDeltaT
// hnamepath: /RF/DeltaT_RF_FirstTime/TOFRF_FirstTimeDeltaT
// hnamepath: /RF/DeltaT_RF_FirstTime/TAGHRF_FirstTimeDeltaT
// hnamepath: /RF/DeltaT_RF_FirstTime/PSCRF_FirstTimeDeltaT

int RFMacro_TDCConversion(void)
{
	gDirectory->cd("/"); //return to file base directory

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("DeltaT_RF_FirstTime");
	TH1I* locHist_FDCRF_FirstTimeDeltaT = (TH1I*)gDirectory->Get("FDCRF_FirstTimeDeltaT");
	TH1I* locHist_TOFRF_FirstTimeDeltaT = (TH1I*)gDirectory->Get("TOFRF_FirstTimeDeltaT");
	TH1I* locHist_TAGHRF_FirstTimeDeltaT = (TH1I*)gDirectory->Get("TAGHRF_FirstTimeDeltaT");
	TH1I* locHist_PSCRF_FirstTimeDeltaT = (TH1I*)gDirectory->Get("PSCRF_FirstTimeDeltaT");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_TDCConversion", "RF_TDCConversion", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFRF_FirstTimeDeltaT != NULL)
	{
		TH1I* locHist = locHist_TOFRF_FirstTimeDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TAGHRF_FirstTimeDeltaT != NULL)
	{
		TH1I* locHist = locHist_TAGHRF_FirstTimeDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PSCRF_FirstTimeDeltaT != NULL)
	{
		TH1I* locHist = locHist_PSCRF_FirstTimeDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FDCRF_FirstTimeDeltaT != NULL)
	{
		TH1I* locHist = locHist_FDCRF_FirstTimeDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

    return 1;
}
