// hnamepath: /HLDetectorTiming/TAGM/TAGMHit Matched time
// hnamepath: /HLDetectorTiming/TAGM/TAGMHit TDC_ADC Difference
// hnamepath: /HLDetectorTiming/TAGH/TAGHHit Matched time
// hnamepath: /HLDetectorTiming/TAGH/TAGHHit TDC_ADC Difference

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
    TH1I* TAGH_Timing        = (TH1I*)gDirectory->Get("TAGH/TAGHHit Matched time");
    TH1I* TAGM_Timing        = (TH1I*)gDirectory->Get("TAGM/TAGMHit Matched time");
	TH2I* TAGH_TDCADC_Timing = (TH2I*)gDirectory->Get("TAGH/TAGHHit TDC_ADC Difference");
	TH2I* TAGM_TDCADC_Timing = (TH2I*)gDirectory->Get("TAGM/TAGMHit TDC_ADC Difference");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TaggerTiming", "TaggerTiming", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGH_Timing != NULL)
	{
        TAGH_Timing->Draw();
        TAGH_Timing->SetFillColor(kGray);
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGH_TDCADC_Timing != NULL)
	{
		TAGH_TDCADC_Timing->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGM_Timing != NULL)
	{
        TAGM_Timing->Draw();
        TAGM_Timing->SetFillColor(kGray);
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGM_TDCADC_Timing != NULL)
	{
		TAGM_TDCADC_Timing->Draw("colz");
	}
}

