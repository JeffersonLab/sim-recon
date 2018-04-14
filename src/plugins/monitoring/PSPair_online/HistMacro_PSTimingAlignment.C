// hnamepath: /HLDetectorTiming/TAGM/TAGMHit Matched time
// hnamepath: /HLDetectorTiming/TAGM/TAGMHit TDC_ADC Difference
// hnamepath: /HLDetectorTiming/TAGH/TAGHHit Matched time
// hnamepath: /HLDetectorTiming/TAGH/TAGHHit TDC_ADC Difference

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("PSPair_online");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH1I* TAGH_PStdiffediff  = (TH1I*)gDirectory->Get("PSPair/PSC_PS_TAGH/PSTAGH_tdiffVsEdiff");
	TH1I* TAGH_PStimeVsE     = (TH1I*)gDirectory->Get("PSPair/PSC_PS_TAGH/PSTAGH_timeVsE");
	TH1I* TAGM_PStdiffediff  = (TH1I*)gDirectory->Get("PSPair/PSC_PS_TAGM/PSTAGM_tdiffVsEdiff");
	TH1I* TAGM_PStimeVsE     = (TH1I*)gDirectory->Get("PSPair/PSC_PS_TAGM/PSTAGM_timeVsE");

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
	//gPad->SetGrid();
	if(TAGH_PStdiffediff != NULL)
	{
	  TAGH_PStdiffediff->GetYaxis()->SetRangeUser(-25.,25.);
	    TAGH_PStdiffediff->Draw("COLZ");
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No Matched PS/TAGH Times");
	  text->Draw();
	}

    locCanvas->cd(2);
    gPad->SetTicks();
    gPad->SetGrid();
    if(TAGH_PStimeVsE != NULL)
    {
      TAGH_PStimeVsE->GetYaxis()->SetRangeUser(-20,20.);
        TAGH_PStimeVsE->Draw("COLZ");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched PS/TAGH Times");
        text->Draw();
    }

    locCanvas->cd(3);
    gPad->SetTicks();
    gPad->SetGrid();
    if(TAGM_PStdiffediff != NULL)
    {
      TAGM_PStdiffediff->GetYaxis()->SetRangeUser(-25.,25.);
        TAGM_PStdiffediff->Draw("COLZ");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched PS/TAGM Times");
        text->Draw();
    }

    locCanvas->cd(4);
    gPad->SetTicks();
    gPad->SetGrid();
    if(TAGM_PStimeVsE != NULL)
    {
	  TAGM_PStimeVsE->GetYaxis()->SetRangeUser(-20.,20.);
        TAGM_PStimeVsE->Draw("COLZ");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched PS/TAGM Times");
        text->Draw();
    }
}

