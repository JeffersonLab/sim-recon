// hnamepath: /HLDetectorTiming/CDC/CDCHit time
// hnamepath: /HLDetectorTiming/FDC/FDCHit Cathode time
// hnamepath: /HLDetectorTiming/FDC/FDCHit Wire time
// hnamepath: /HLDetectorTiming/TRACKING/Earliest CDC Time Minus Matched SC Time
// hnamepath: /HLDetectorTiming/TRACKING/Earliest Flight-time Corrected CDC Time
// hnamepath: /HLDetectorTiming/TRACKING/Earliest Flight-time Corrected FDC Time

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
    TH1I* CDC_Timing     = (TH1I*)gDirectory->Get("CDC/CDCHit time");
	TH1I* CDC_SC_Timing  = (TH1I*)gDirectory->Get("TRACKING/Earliest CDC Time Minus Matched SC Time");
	TH1I* CDC_Earliest_Time  = (TH1I*)gDirectory->Get("TRACKING/Earliest Flight-time Corrected CDC Time");
    TH1I* FDC_Strip_Timing     = (TH1I*)gDirectory->Get("FDC/FDCHit Cathode time");
    TH1I* FDC_Wire_Timing     = (TH1I*)gDirectory->Get("FDC/FDCHit Wire time");
	TH1I* FDC_Earliest_Time  = (TH1I*)gDirectory->Get("TRACKING/Earliest Flight-time Corrected FDC Time");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TrackMatchedTiming", "TrackMatchedTiming", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(CDC_Timing != NULL)
	{
        CDC_Timing->Draw();
        CDC_Timing->SetFillColor(kGray);
	}
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No CDC Times in data");
        text->Draw();
    }

    locCanvas->cd(2);
    gPad->SetTicks();
    gPad->SetGrid();
    if(CDC_Earliest_Time != NULL)
    {
        CDC_Earliest_Time->Draw();
        CDC_Earliest_Time->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No CDC tracks matced to BCAL/SC with reasonable FOM");
        text->Draw();
    }
    
    locCanvas->cd(3);
    gPad->SetTicks();
    gPad->SetGrid();
    if(CDC_SC_Timing != NULL)
    {
        CDC_SC_Timing->Draw();
        CDC_SC_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No CDC tracks matced to SC with reasonable FOM");
        text->Draw();
    }

    locCanvas->cd(4);
    gPad->SetTicks();
    gPad->SetGrid();
    if(FDC_Strip_Timing != NULL)
    {
        FDC_Strip_Timing->Draw();
        FDC_Strip_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No FDC Cathode Hit times");
        text->Draw();
    }

    locCanvas->cd(5);
    gPad->SetTicks();
    gPad->SetGrid();
    if(FDC_Wire_Timing != NULL)
    {
        FDC_Wire_Timing->Draw();
        FDC_Wire_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No FDC Wire Hit times");
        text->Draw();
    }

    locCanvas->cd(6);
    gPad->SetTicks();
    gPad->SetGrid();
    if(FDC_Earliest_Time != NULL)
    {
        FDC_Earliest_Time->Draw();
        FDC_Earliest_Time->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No FDC tracks matced to SC/TOF with reasonable FOM");
        text->Draw();
    }
}

