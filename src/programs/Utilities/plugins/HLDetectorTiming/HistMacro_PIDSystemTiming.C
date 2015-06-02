// hnamepath: /HLDetectorTiming/SC/SCHit Matched time
// hnamepath: /HLDetectorTiming/SC/SCHit TDC_ADC Difference
// hnamepath: /HLDetectorTiming/TOF/TOFHit Matched time
// hnamepath: /HLDetectorTiming/TOF/TOFHit TDC_ADC Difference

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
    TH1I* SC_Timing          = (TH1I*)gDirectory->Get("SC/SCHit Matched time");
    TH1I* TOF_Timing         = (TH1I*)gDirectory->Get("TOF/TOFHit Matched time");
	TH2I* SC_TDCADC_Timing   = (TH2I*)gDirectory->Get("SC/SCHit TDC_ADC Difference");
	TH2I* TOF_TDCADC_Timing  = (TH2I*)gDirectory->Get("TOF/TOFHit TDC_ADC Difference");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("PIDSystemTiming", "PIDSystemTiming", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(SC_Timing != NULL)
	{
        SC_Timing->Draw();
        SC_Timing->SetFillColor(kGray);
	}
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched SC TDC/ADC Times");
        text->Draw();
    }

    locCanvas->cd(2);
    gPad->SetTicks();
    gPad->SetGrid();
    if(SC_TDCADC_Timing != NULL)
    {
        SC_TDCADC_Timing->Draw("colz");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched SC TDC/ADC Times");
        text->Draw();
    }

    locCanvas->cd(3);
    gPad->SetTicks();
    gPad->SetGrid();
    if(TOF_Timing != NULL)
    {
        TOF_Timing->Draw();
        TOF_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched TOF TDC/ADC Times");
        text->Draw();
    }

    locCanvas->cd(4);
    gPad->SetTicks();
    gPad->SetGrid();
    if(TOF_TDCADC_Timing != NULL)
    {
        TOF_TDCADC_Timing->Draw("colz");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched TOF TDC/ADC Times");
        text->Draw();
    }
}

