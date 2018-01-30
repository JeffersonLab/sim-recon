// hnamepath: /HLDetectorTiming/FCAL/FCALHit time
// hnamepath: /HLDetectorTiming/BCAL/BCALHit TDC time
// hnamepath: /HLDetectorTiming/BCAL/BCALHit Upstream Per Channel TDC-ADC Hit Time
// hnamepath: /HLDetectorTiming/BCAL/BCALHit Downstream Per Channel TDC-ADC Hit Time

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
    TH1I* FCAL_Timing          = (TH1I*)gDirectory->Get("FCAL/FCALHit time");
    TH1I* BCAL_Timing          = (TH1I*)gDirectory->Get("BCAL/BCALHit TDC time");
    TH2I* BCAL_U_TDCADC_Timing = (TH2I*)gDirectory->Get("BCAL/BCALHit Upstream Per Channel TDC-ADC Hit Time");
    TH2I* BCAL_D_TDCADC_Timing = (TH2I*)gDirectory->Get("BCAL/BCALHit Downstream Per Channel TDC-ADC Hit Time");

    //Get/Make Canvas
    TCanvas *locCanvas = NULL;
    if(TVirtualPad::Pad() == NULL)
        locCanvas = new TCanvas("CalorimeterTiming", "CalorimeterTiming", 1200, 800); //for testing
    else
        locCanvas = gPad->GetCanvas();
    locCanvas->Divide(2, 2);

    //Draw
    locCanvas->cd(1);
    gPad->SetTicks();
    gPad->SetGrid();
    if(FCAL_Timing != NULL)
    {
        FCAL_Timing->Draw();
        FCAL_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No FCAL Hit Times");
        text->Draw();
    }

    locCanvas->cd(2);
    gPad->SetTicks();
    gPad->SetGrid();
    if(BCAL_U_TDCADC_Timing != NULL)
    {
        BCAL_U_TDCADC_Timing->Draw("colz");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched Upstream BCAL TDC/ADC Times");
        text->Draw();
    }

    locCanvas->cd(3);
    gPad->SetTicks();
    gPad->SetGrid();
    if(BCAL_Timing != NULL)
    {
        BCAL_Timing->Draw();
        BCAL_Timing->SetFillColor(kGray);
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No BCAL TDC Times");
        text->Draw();
    }

    locCanvas->cd(4);
    gPad->SetTicks();
    gPad->SetGrid();
    if(BCAL_D_TDCADC_Timing != NULL)
    {
        BCAL_D_TDCADC_Timing->Draw("colz");
    }
    else{
        TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
        text->AddText("No Matched BCAL Downstream TDC/ADC Times");
        text->Draw();
    }
}

