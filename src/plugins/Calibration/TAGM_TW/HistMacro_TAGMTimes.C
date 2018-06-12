// hnamepath:  /TAGM_TW/adc_rf_all
// hnamepath:  /TAGM_TW/tdc_adc_all
// hnamepath:  /TAGM_TW/t_adc_all
{
	//Goto Path
        TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("TAGM_TW");
        //TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("st_Tresolution");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2I* TAGM_ADC_RF_Times = (TH2I*)gDirectory->Get("adc_rf_all");
	TH2I* TAGM_TDC_ADC_Times = (TH2I*)gDirectory->Get("tdc_adc_all");
	TH2I* TAGM_T_ADC_Times = (TH2I*)gDirectory->Get("t_adc_all");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("TAGM_TW", "TAGM_TW", 1200, 600); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3,1);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGM_ADC_RF_Times != NULL)
	{	  
	    TAGM_ADC_RF_Times->SetStats(0);
	    //TAGM_ADC_RF_Times->SetTitle("Corrected Time vs. Z  -  Sector 3");
	    //TAGM_ADC_RF_Times->GetYaxis()->SetTitle("Corrected Time - RF Time (ns)");
	    TAGM_ADC_RF_Times->Draw("COLZ");
	    locCanvas->Update();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No TAGM Hits!");
	  text->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGM_TDC_ADC_Times != NULL)
	{	  
	    TAGM_TDC_ADC_Times->SetStats(0);
	    //TAGM_TDC_ADC_Times->SetTitle("Corrected Time vs. Z  -  Sector 3");
	    //TAGM_TDC_ADC_Times->GetYaxis()->SetTitle("Corrected Time - RF Time (ns)");
	    TAGM_TDC_ADC_Times->Draw("COLZ");
	    locCanvas->Update();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No TAGM Hits!");
	  text->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(TAGM_T_ADC_Times != NULL)
	{	  
	    TAGM_T_ADC_Times->SetStats(0);
	    //TAGM_T_ADC_Times->SetTitle("Corrected Time vs. Z  -  Sector 3");
	    //TAGM_T_ADC_Times->GetYaxis()->SetTitle("Corrected Time - RF Time (ns)");
	    TAGM_T_ADC_Times->Draw("COLZ");
	    locCanvas->Update();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No TAGM Hits!");
	  text->Draw();
	}



}

