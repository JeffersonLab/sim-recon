// hnamepath: /HLDetectorTiming/TRACKING/TAGM - SC Target Time
// hnamepath: /HLDetectorTiming/TRACKING/TAGH - SC Target Time
// hnamepath: /HLDetectorTiming/TRACKING/Tagger - SC Target Time
// hnamepath: /HLDetectorTiming/TRACKING/Tagger - SC 1D Target Time

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
    TH1I* Tagger_SC_1D = (TH1I*)gDirectory->Get("TRACKING/Tagger - SC 1D Target Time");
    TH2I* Tagger_SC_2D = (TH2I*)gDirectory->Get("TRACKING/Tagger - SC Target Time");
    TH2I* TAGH_SC_2D   = (TH2I*)gDirectory->Get("TRACKING/TAGH - SC Target Time");
    TH2I* TAGM_SC_2D   = (TH2I*)gDirectory->Get("TRACKING/TAGM - SC Target Time");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TaggerSCAlignment", "TaggerSCAlignment", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGH_SC_2D != NULL)
	{
        TAGH_SC_2D->Draw("colz");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(Tagger_SC_2D != NULL)
	{
		Tagger_SC_2D->Draw("colz");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGM_SC_2D != NULL)
	{
        TAGM_SC_2D->Draw("colz");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(Tagger_SC_1D != NULL)
	{
		Tagger_SC_1D->Draw();
        Tagger_SC_1D->SetFillColor(kGray);
	}
}

