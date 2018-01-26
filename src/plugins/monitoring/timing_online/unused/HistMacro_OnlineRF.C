// hnamepath: /HLDetectorTiming/TRACKING/TAGM - RFBunch 1D Time
// hnamepath: /HLDetectorTiming/TRACKING/Tagger - RFBunch 1D Time
// hnamepath: /HLDetectorTiming/TRACKING/SC - RF Time
// hnamepath: /HLDetectorTiming/TRACKING/FCAL - RF Time
// hnamepath: /HLDetectorTiming/TRACKING/TOF - RF Time
// hnamepath: /HLDetectorTiming/TRACKING/BCAL - RF Time

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH1I* TAGM_RF_Timing     = (TH1I*)gDirectory->Get("TRACKING/TAGM - RFBunch 1D Time");
	TH1I* Tagger_RF_Timing     = (TH1I*)gDirectory->Get("TRACKING/Tagger - RFBunch 1D Time");
	TH1I* SC_RF_Timing  = (TH1I*)gDirectory->Get("TRACKING/SC - RF Time");
	TH1I* FCAL_RF_Timing = (TH1I*)gDirectory->Get("TRACKING/FCAL - RF Time");
	TH1I* TOF_RF_Timing  = (TH1I*)gDirectory->Get("TRACKING/TOF - RF Time");
	TH1I* BCAL_RF_Timing = (TH1I*)gDirectory->Get("TRACKING/BCAL - RF Time");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("OnlineRF", "OnlineRF", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(SC_RF_Timing != NULL)
	{
	    SC_RF_Timing->Draw();
	    SC_RF_Timing->SetFillColor(kGray);
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No CDC tracks matched to SC with reasonable FOM");
	  text->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TOF_RF_Timing != NULL)
	{
	  TOF_RF_Timing->Draw();
	  TOF_RF_Timing->SetFillColor(kGray);
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No track matches to TOF with reasonable FOM");
	  text->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(BCAL_RF_Timing != NULL)
	{
	  BCAL_RF_Timing->Draw();
	  BCAL_RF_Timing->SetFillColor(kGray);
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No track matches to BCAL with reasonable FOM");
	  text->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(FCAL_RF_Timing != NULL)
	{
	  FCAL_RF_Timing->Draw();
	  FCAL_RF_Timing->SetFillColor(kGray);
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No track matches to FCAL with reasonable FOM");
	  text->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(Tagger_RF_Timing != NULL)
	{
	  Tagger_RF_Timing->Draw();
	  Tagger_RF_Timing->SetFillColor(kGray);
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No Matches to SC and TOF with reasonable FOM");
	  text->Draw();
	}
	
	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGM_RF_Timing != NULL)
	{
	  TAGM_RF_Timing->Draw();
	  TAGM_RF_Timing->SetFillColor(kGray);
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No Matches to SC and BCAL with reasonable FOM");
	  text->Draw();
	}
}

