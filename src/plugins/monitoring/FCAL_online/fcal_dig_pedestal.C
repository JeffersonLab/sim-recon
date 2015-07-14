// hnamepath:  /fcal/digOcc2D
// hnamepath:  /fcal/digPed2D
// hnamepath:  /fcal/digPedSq2D
// hnamepath:  /fcal/digPed
// hnamepath:  /fcal/digPedChan

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("fcal");
  if(dir) dir->cd();

  TH2F* digOcc2D = (TH2F*)gDirectory->FindObjectAny("digOcc2D");
  TH2F* digPed2D = (TH2F*)gDirectory->FindObjectAny("digPed2D");
  TH2F* digPedSq2D = (TH2F*)gDirectory->FindObjectAny("digPedSq2D");
  TH1I* digPed = (TH1I*)gDirectory->FindObjectAny("digPed");
  TProfile* digPedChan = (TProfile*)gDirectory->FindObjectAny("digPedChan");
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL Monitor", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide( 2, 2 );
  
  if( digPed ){

    digPed->SetStats( 0 );
    digPed->SetFillColor( kRed );
    c1->cd( 1 );
    digPed->Draw();
  }
  
  if( digPedChan ){

    digPedChan->SetStats( 0 );
    c1->cd( 2 );
    digPedChan->Draw();
  }

  if( digPed2D && digOcc2D && digPed ){

    TH2F* digPed2DAvg = (TH2F*)digPed2D->Clone( "digPed2DAvg" );
    digPed2DAvg->Divide( digOcc2D );
    double avgPed = digPed->GetMean();
    digPed2DAvg->SetTitle( "FCAL Pedestal - Average Pedestal" );

    for( int x = 1; x <= digPed2DAvg->GetNbinsX(); ++x ){
      for( int y = 1; y <= digPed2DAvg->GetNbinsY(); ++y ){

	digPed2DAvg->SetBinContent( x, y, digPed2DAvg->GetBinContent( x, y ) - avgPed );
      }
    }

    digPed2DAvg->SetStats( 0 );
    digPed2DAvg->SetMinimum( -0.2*avgPed );
    digPed2DAvg->SetMaximum(  0.2*avgPed );
    c1->cd( 3 );
    digPed2DAvg->Draw( "colz" );
  }

  if( digPed2D && digPedSq2D && digOcc2D && digPed ){

    TH2F* digPed2DAvg = (TH2F*)digPed2D->Clone( "digPed2DAvg" );
    TH2F* digPed2DRMS = (TH2F*)digPed2D->Clone( "digPed2DRMS" );
    TH2F* digPedSq2DAvg = (TH2F*)digPedSq2D->Clone( "digPedSq2DAvg" );
    digPed2DAvg->Divide( digOcc2D );
    digPedSq2DAvg->Divide( digOcc2D );

    digPed2DRMS->SetTitle( "FCAL Pedestal RMS [ADC Counts]" );

    for( int x = 1; x <= digPed2DAvg->GetNbinsX(); ++x ){
      for( int y = 1; y <= digPed2DAvg->GetNbinsY(); ++y ){

	double var = digPedSq2DAvg->GetBinContent( x, y );
	var -= ( digPed2DAvg->GetBinContent( x, y ) * 
		 digPed2DAvg->GetBinContent( x, y ) );
	
	if( digOcc2D->GetBinContent( x, y ) != 0 ){

	  digPed2DRMS->SetBinContent( x, y, TMath::Sqrt( var ) );
	}
	else{
	  
	  // set these below minimum so that they appear
	  // white in the histogram
	  digPed2DRMS->SetBinContent( x, y, -1 );
	}
      }
    }

    digPed2DRMS->SetStats( 0 );
    digPed2DRMS->SetMaximum( 10 );
    digPed2DRMS->SetMinimum( 0 );
    c1->cd( 4 );
    digPed2DRMS->Draw( "colz" );
  }
}
