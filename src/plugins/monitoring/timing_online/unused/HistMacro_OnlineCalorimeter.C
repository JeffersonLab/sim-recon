// hnamepath: /HLDetectorTiming/BCAL/BCALHit ADC time
// hnamepath: /HLDetectorTiming/BCAL/BCALHit TDC time
// hnamepath: /HLDetectorTiming/BCAL/BCALHit Downstream Per Channel TDC-ADC Hit Time
// hnamepath: /HLDetectorTiming/BCAL/BCALHit Upstream Per Channel TDC-ADC Hit Time
// hnamepath: /HLDetectorTiming/FCAL/FCALHit time
// hnamepath: /HLDetectorTiming/FCAL/FCALHit Local Time
// hnamepath: /HLDetectorTiming/FCAL/FCALHit Occupancy

{
	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("HLDetectorTiming");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Setpoints
	double nominalFCALTime = 17.;
	double nominalBCALADCTime = 17.;
	double nominalBCALTDCTime = 17.;

	//Get Histograms
	TH1I* BCAL_ADC_Timing = (TH1I*)gDirectory->Get("BCAL/BCALHit ADC time");
	TH1I* BCAL_TDC_Timing = (TH1I*)gDirectory->Get("BCAL/BCALHit TDC time");
	TH1I* BCAL_TDC_ADC_DS_Timing = (TH1I*)gDirectory->Get("BCAL/BCALHit Downstream Per Channel TDC-ADC Hit Time");
	TH1I* BCAL_TDC_ADC_US_Timing = (TH1I*)gDirectory->Get("BCAL/BCALHit Upstream Per Channel TDC-ADC Hit Time");
	TH1I* FCAL_ADC_Timing = (TH1I*)gDirectory->Get("FCAL/FCALHit time");
	TH1I* FCAL_Local_Timing = (TH1I*)gDirectory->Get("FCAL/FCALHit Local Time");
	TH1I* FCAL_Occupancy = (TH1I*)gDirectory->Get("FCAL/FCALHit Occupancy");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("OnlineCalorimeter", "OnlineCalorimeter", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(FCAL_ADC_Timing != NULL)
	{
	    FCAL_ADC_Timing->GetXaxis()->SetRangeUser(-200,200);
	    FCAL_ADC_Timing->Draw();
	    FCAL_ADC_Timing->SetFillColor(kGray);
	    locCanvas->Update();

	    TLine *ln = new TLine(nominalFCALTime, gPad->GetUymin(), nominalFCALTime, gPad->GetUymax());
	    ln->SetLineColor(2);
	    ln->Draw();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No FCAL ADC hits!");
	  text->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(BCAL_ADC_Timing != NULL)
	{
	    BCAL_ADC_Timing->GetXaxis()->SetRangeUser(-200,200);
	    BCAL_ADC_Timing->Draw();
	    BCAL_ADC_Timing->SetFillColor(kGray);
	    locCanvas->Update();

	    TLine *ln = new TLine(nominalFCALTime, gPad->GetUymin(), nominalBCALADCTime, gPad->GetUymax());
	    ln->SetLineColor(2);
	    ln->Draw();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No BCAL ADC hits!");
	  text->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(BCAL_TDC_Timing != NULL)
	{
	    BCAL_TDC_Timing->GetXaxis()->SetRangeUser(-200,200);
	    BCAL_TDC_Timing->Draw();
	    BCAL_TDC_Timing->SetFillColor(kGray);
	    locCanvas->Update();

	    TLine *ln = new TLine(nominalFCALTime, gPad->GetUymin(), nominalBCALTDCTime, gPad->GetUymax());
	    ln->SetLineColor(2);
	    ln->Draw();
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No BCAL TDC hits!");
	  text->Draw();
	}


	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(FCAL_Local_Timing != NULL && FCAL_Occupancy != NULL)
	{
	    TH2F* FCAL_Avg_Timing = (TH2F*)FCAL_Local_Timing->Clone( "FCALHit Local Time" );
	    FCAL_Avg_Timing->Divide( FCAL_Occupancy );

	    for( int x = 1; x <= FCAL_Avg_Timing->GetNbinsX(); ++x ){
	      for( int y = 1; y <= FCAL_Avg_Timing->GetNbinsY(); ++y ){

		if( FCAL_Occupancy->GetBinContent( x, y ) == 0 ){

		  // set this off scale low so unused blocks
		  // appear white in the plot
		  FCAL_Avg_Timing->SetBinContent( x, y, -1E3 );
		}
	      }
	    }

	    FCAL_Avg_Timing->SetMinimum(-10 );
	    FCAL_Avg_Timing->SetMaximum( 10 );
	    FCAL_Avg_Timing->SetStats( 0 );
	    FCAL_Avg_Timing->Draw( "colz" );
	  
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No FCAL ADC hits!");
	  text->Draw();
	}


	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(BCAL_TDC_ADC_DS_Timing != NULL)
	{
	    BCAL_TDC_ADC_DS_Timing->Draw("COLZ");
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No BCAL ADC and TDC hits!");
	  text->Draw();
	}


	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(BCAL_TDC_ADC_US_Timing != NULL)
	{
	    BCAL_TDC_ADC_US_Timing->Draw("COLZ");
	}
	else{
	  TPaveText *text = new TPaveText(0.1, 0.4, 0.9, 0.6);
	  text->AddText("No BCAL ADC and TDC hits!");
	  text->Draw();
	}


}

