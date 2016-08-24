// hnamepath: /RF/DeltaT_RF_TAGH/FDCRF_TaggerDeltaT
// hnamepath: /RF/DeltaT_RF_TAGH/TOFRF_TaggerDeltaT
// hnamepath: /RF/DeltaT_RF_TAGH/TAGHRF_TaggerDeltaT
// hnamepath: /RF/DeltaT_RF_TAGH/PSCRF_TaggerDeltaT

int RFMacro_TaggerComparison(void)
{
	gStyle->SetOptStat(1111);
	gDirectory->cd("/"); //return to file base directory

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
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
		locCanvas = new TCanvas("RF_TaggerDeltaT", "RF_TaggerDeltaT", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFRF_TaggerDeltaT != NULL)
	{
		TH1I* locHist = locHist_TOFRF_TaggerDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TAGHRF_TaggerDeltaT != NULL)
	{
		TH1I* locHist = locHist_TAGHRF_TaggerDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PSCRF_TaggerDeltaT != NULL)
	{
		TH1I* locHist = locHist_PSCRF_TaggerDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FDCRF_TaggerDeltaT != NULL)
	{
		TH1I* locHist = locHist_FDCRF_TaggerDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

    return 1;
}
