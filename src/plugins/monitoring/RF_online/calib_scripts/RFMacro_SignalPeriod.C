// hnamepath: /RF/RF_SignalPeriod/FDCRF_SignalPeriod
// hnamepath: /RF/RF_SignalPeriod/TOFRF_SignalPeriod
// hnamepath: /RF/RF_SignalPeriod/TAGHRF_SignalPeriod
// hnamepath: /RF/RF_SignalPeriod/PSCRF_SignalPeriod

int RFMacro_SignalPeriod(void)
{
	gStyle->SetOptStat(1111);
	gDirectory->cd("/"); //return to file base directory

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("RF_SignalPeriod");
	TH1I* locHist_FDCRF_SignalPeriod = (TH1I*)gDirectory->Get("FDCRF_SignalPeriod");
	TH1I* locHist_TOFRF_SignalPeriod = (TH1I*)gDirectory->Get("TOFRF_SignalPeriod");
	TH1I* locHist_TAGHRF_SignalPeriod = (TH1I*)gDirectory->Get("TAGHRF_SignalPeriod");
	TH1I* locHist_PSCRF_SignalPeriod = (TH1I*)gDirectory->Get("PSCRF_SignalPeriod");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_SignalPeriod", "RF_SignalPeriod", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFRF_SignalPeriod != NULL)
	{
		TH1I* locHist = locHist_TOFRF_SignalPeriod;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TAGHRF_SignalPeriod != NULL)
	{
		TH1I* locHist = locHist_TAGHRF_SignalPeriod;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PSCRF_SignalPeriod != NULL)
	{
		TH1I* locHist = locHist_PSCRF_SignalPeriod;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FDCRF_SignalPeriod != NULL)
	{
		TH1I* locHist = locHist_FDCRF_SignalPeriod;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

    return 1;
}
