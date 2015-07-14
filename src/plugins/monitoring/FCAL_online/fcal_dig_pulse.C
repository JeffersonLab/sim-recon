// hnamepath:  /fcal/digOcc2D
// hnamepath:  /fcal/digPeakV2D
// hnamepath:  /fcal/digPeakV
// hnamepath:  /fcal/digN

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("fcal");
  if(dir) dir->cd();

  TH2F* digOcc2D   = (TH2F*)gDirectory->FindObjectAny("digOcc2D");
  TH1I* digN       = (TH1I*)gDirectory->FindObjectAny("digN");
  TH2F* digPeakV2D = (TH2F*)gDirectory->FindObjectAny("digPeakV2D");
  TH1I* digPeakV   = (TH1I*)gDirectory->FindObjectAny("digPeakV");

  double nEvents = ( digN ? digN->GetEntries() : 0 );
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL Monitor", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide( 2, 2 );

  if( digN ){

    digN->SetStats( 0 );
    digN->SetFillColor( kRed );
    c1->cd( 1 );
    digN->Draw();
  }
  
  if( digOcc2D ){

    TH2F* digOcc2DAvg = (TH2F*)digOcc2D->Clone( "digOcc2DAvg" );

    digOcc2DAvg->SetTitle( "FCAL Pulse Occupancy per Event" );

    for( int x = 1; x <= digOcc2DAvg->GetNbinsX(); ++x ){
      for( int y = 1; y <= digOcc2DAvg->GetNbinsY(); ++y ){

	digOcc2DAvg->SetBinContent( x, y, digOcc2DAvg->GetBinContent( x, y ) / nEvents );
      }
    }

    digOcc2DAvg->SetStats( 0 );
    c1->cd( 2 );
    digOcc2DAvg->Draw( "colz" );
  }

  if( digPeakV ){

    digPeakV->SetStats( 0 );
    digPeakV->SetFillColor( kRed );
    c1->cd( 3 );
    digPeakV->Draw();
  }

  if( digPeakV2D && digOcc2D ){

    TH2F* digPeakV2DAvg = (TH2F*)digPeakV2D->Clone( "digPeakV2DAvg" );
    digPeakV2DAvg->Divide( digOcc2D );
    digPeakV2DAvg->SetMinimum( 0 );
    digPeakV2DAvg->SetMaximum( 2 );
    digPeakV2DAvg->SetStats( 0 );
    c1->cd( 4 );
    digPeakV2DAvg->Draw( "colz" );
  }

}
