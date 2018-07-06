// hnamepath:  /ST_Tresolution/h2_CorrectedTime_z_3
// hnamepath:  /ST_Tresolution/h2_CorrectedTime_z_10
// hnamepath:  /ST_Tresolution/h2_CorrectedTime_z_18
// hnamepath:  /ST_Tresolution/h2_CorrectedTime_z_25
{
	//Goto Path
        TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("ST_Tresolution");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	TH2I* ST_Paddle3_Resolution = (TH2I*)gDirectory->Get("h2_CorrectedTime_z_3");
	TH2I* ST_Paddle10_Resolution = (TH2I*)gDirectory->Get("h2_CorrectedTime_z_10");
	TH2I* ST_Paddle18_Resolution = (TH2I*)gDirectory->Get("h2_CorrectedTime_z_18");
	TH2I* ST_Paddle25_Resolution = (TH2I*)gDirectory->Get("h2_CorrectedTime_z_25");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("STResolution", "STResolution", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(ST_Paddle3_Resolution != NULL)
	{	  
	    ST_Paddle3_Resolution->Rebin2D(5,1);
	    ST_Paddle3_Resolution->SetStats(0);
	    ST_Paddle3_Resolution->SetTitle("Corrected Time vs. Z  -  Sector 3");
	    ST_Paddle3_Resolution->GetYaxis()->SetRangeUser(-2.,2.);
	    ST_Paddle3_Resolution->GetYaxis()->SetTitle("Corrected Time - RF Time (ns)");
	    ST_Paddle3_Resolution->Draw("COLZ");
	    locCanvas->Update();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No hits in SC Sector 3!");
	  text->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(ST_Paddle10_Resolution != NULL)
	{
	    ST_Paddle10_Resolution->Rebin2D(5,1);
	    ST_Paddle10_Resolution->SetStats(0);
	    ST_Paddle10_Resolution->SetTitle("Corrected Time vs. Z  -  Sector 10");
	    ST_Paddle10_Resolution->GetYaxis()->SetRangeUser(-2.,2.);
	    ST_Paddle10_Resolution->GetYaxis()->SetTitle("Corrected Time - RF Time (ns)");
	    ST_Paddle10_Resolution->Draw("COLZ");
	    locCanvas->Update();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No hits in SC Sector 10!");
	  text->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(ST_Paddle18_Resolution != NULL)
	{
	    ST_Paddle18_Resolution->Rebin2D(5,1);
	    ST_Paddle18_Resolution->SetStats(0);
	    ST_Paddle18_Resolution->SetTitle("Corrected Time vs. Z  -  Sector 18");
	    ST_Paddle18_Resolution->GetYaxis()->SetRangeUser(-2.,2.);
	    ST_Paddle18_Resolution->GetYaxis()->SetTitle("Corrected Time - RF Time (ns)");
	    ST_Paddle18_Resolution->Draw("COLZ");
	    locCanvas->Update();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No hits in SC Sector 18!");
	  text->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(ST_Paddle25_Resolution != NULL)
	{
	    ST_Paddle25_Resolution->Rebin2D(5,1);
	    ST_Paddle25_Resolution->SetStats(0);
	    ST_Paddle25_Resolution->SetTitle("Corrected Time vs. Z  -  Sector 25");
	    ST_Paddle25_Resolution->GetYaxis()->SetRangeUser(-2.,2.);
	    ST_Paddle25_Resolution->GetYaxis()->SetTitle("Corrected Time - RF Time (ns)");
	    ST_Paddle25_Resolution->Draw("COLZ");
	    locCanvas->Update();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No hits in SC Sector 25!");
	  text->Draw();
	}

}

