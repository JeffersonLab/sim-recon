// hnamepath: /RF/RF_BeamBunchPeriod/RFBeamBunchPeriod

int RFMacro_BeamBunchPeriod(void)
{
	gDirectory->cd("/"); //return to file base directory

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF Bunch Period Histograms
	gDirectory->cd("RF_BeamBunchPeriod");
	TH1I* locHist_RFBeamBunchPeriod = (TH1I*)gDirectory->Get("RFBeamBunchPeriod");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_BeamBunchPeriod", "RF_BeamBunchPeriod", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();

	//Draw
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFBeamBunchPeriod != NULL)
	{
		TH1I* locHist = locHist_RFBeamBunchPeriod;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

    return 1;
}
