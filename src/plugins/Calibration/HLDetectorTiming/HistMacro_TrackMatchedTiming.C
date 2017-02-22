// hnamepath: /HLDetectorTiming/CDC/CDCHit time
// hnamepath: /HLDetectorTiming/FDC/FDCHit time
// hnamepath: /HLDetectorTiming/TRACKING/Earliest CDC Time Minus Matched SC Time
// hnamepath: /HLDetectorTiming/TRACKING/FCAL - SC Target Time
// hnamepath: /HLDetectorTiming/TRACKING/TOF - SC Target Time
// hnamepath: /HLDetectorTiming/TRACKING/BCAL - SC Target Time

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
    TH1I* CDC_Timing     = (TH1I*)gDirectory->Get("CDC/CDCHit time");
    TH1I* FDC_Timing     = (TH1I*)gDirectory->Get("FDC/FDCHit Cathode time");
	TH1I* CDC_SC_Timing  = (TH1I*)gDirectory->Get("TRACKING/Earliest CDC Time Minus Matched SC Time");
	TH1I* FCAL_SC_Timing = (TH1I*)gDirectory->Get("TRACKING/FCAL - SC Target Time");
    TH1I* TOF_SC_Timing  = (TH1I*)gDirectory->Get("TRACKING/TOF - SC Target Time");
    TH1I* BCAL_SC_Timing = (TH1I*)gDirectory->Get("TRACKING/BCAL - SC Target Time");

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

    locCanvas->cd(3);
    gPad->SetTicks();
    gPad->SetGrid();
    if(FCAL_SC_Timing != NULL)
    {
        FCAL_SC_Timing->Draw();
        FCAL_SC_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matches to SC and FCAL with reasonable FOM");
        text->Draw();
    }

    locCanvas->cd(4);
    gPad->SetTicks();
    gPad->SetGrid();
    if(FDC_Timing != NULL)
    {
        FDC_Timing->Draw();
        FDC_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No FDC Cathode Hit times");
        text->Draw();
    }

    locCanvas->cd(5);
    gPad->SetTicks();
    gPad->SetGrid();
    if(TOF_SC_Timing != NULL)
    {
        TOF_SC_Timing->Draw();
        TOF_SC_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matches to SC and TOF with reasonable FOM");
        text->Draw();
    }

    locCanvas->cd(6);
    gPad->SetTicks();
    gPad->SetGrid();
    if(BCAL_SC_Timing != NULL)
    {
        BCAL_SC_Timing->Draw();
        BCAL_SC_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matches to SC and BCAL with reasonable FOM");
        text->Draw();
    }
}

