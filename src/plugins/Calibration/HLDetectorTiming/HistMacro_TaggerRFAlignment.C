// hnamepath: /HLDetectorTiming/TRACKING/TAGM - RFBunch Time
// hnamepath: /HLDetectorTiming/TRACKING/TAGH - RFBunch Time
// hnamepath: /HLDetectorTiming/TRACKING/Tagger - RFBunch Time
// hnamepath: /HLDetectorTiming/TRACKING/Tagger - RFBunch 1D Time

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
    TH1I* Tagger_RF_1D = (TH1I*)gDirectory->Get("TRACKING/Tagger - RFBunch 1D Time");
    TH2I* Tagger_RF_2D = (TH2I*)gDirectory->Get("TRACKING/Tagger - RFBunch Time");
    TH2I* TAGH_RF_2D   = (TH2I*)gDirectory->Get("TRACKING/TAGH - RFBunch Time");
    TH2I* TAGM_RF_2D   = (TH2I*)gDirectory->Get("TRACKING/TAGM - RFBunch Time");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TaggerRFAlignment", "TaggerRFAlignment", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGH_RF_2D != NULL)
	{
        TAGH_RF_2D->Draw("colz");
	}
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched TAGH/RFBunch Times");
        text->Draw();
    }

    locCanvas->cd(2);
    gPad->SetTicks();
    gPad->SetGrid();
    if(Tagger_RF_2D != NULL)
    {
        Tagger_RF_2D->Draw("colz");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched Tagger/RFBunch Times");
        text->Draw();
    }

    locCanvas->cd(3);
    gPad->SetTicks();
    gPad->SetGrid();
    if(TAGM_RF_2D != NULL)
    {
        TAGM_RF_2D->Draw("colz");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched TAGM/RFBunch Times");
        text->Draw();
    }

    locCanvas->cd(4);
    gPad->SetTicks();
    gPad->SetGrid();
    if(Tagger_RF_1D != NULL)
    {
        Tagger_RF_1D->Draw();
        Tagger_RF_1D->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched Tagger/RFBunch Times");
        text->Draw();
    }
}

