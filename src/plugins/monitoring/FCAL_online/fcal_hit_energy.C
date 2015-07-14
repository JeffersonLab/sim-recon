// hnamepath:  /fcal/hitOcc2D
// hnamepath:  /fcal/hitN
// hnamepath:  /fcal/hitETot
// hnamepath:  /fcal/hitE2D

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("fcal");
  if(dir) dir->cd();

  TH2F* hitOcc2D = (TH2F*)gDirectory->FindObjectAny("hitOcc2D");
  TH1I* hitN     = (TH1I*)gDirectory->FindObjectAny("hitN");
  TH1I* hitETot  = (TH1I*)gDirectory->FindObjectAny("hitETot");
  TH2F* hitE2D   = (TH2F*)gDirectory->FindObjectAny("hitE2D");

  double nEvents = ( hitETot ? hitETot->GetEntries() : 0 );
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL Monitor", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide( 2, 2 );

  if( hitN ){

    hitN->SetStats( 0 );
    hitN->SetFillColor( kYellow );
    c1->cd( 1 );
    hitN->Draw();
  }
  
  if( hitETot ){

    hitETot->SetStats( 0 );
    hitETot->SetFillColor( kYellow );
    c1->cd( 2 );
    hitETot->Draw();
  }

  if( hitOcc2D ){

    TH2F* hitOcc2DAvg = (TH2F*)hitOcc2D->Clone( "hitOcc2DAvg" );

    hitOcc2DAvg->SetTitle( "FCAL Hit Occupancy per Event" );

    for( int x = 1; x <= hitOcc2DAvg->GetNbinsX(); ++x ){
      for( int y = 1; y <= hitOcc2DAvg->GetNbinsY(); ++y ){

	hitOcc2DAvg->SetBinContent( x, y, hitOcc2DAvg->GetBinContent( x, y ) / nEvents );
      }
    }

    hitOcc2DAvg->SetStats( 0 );
    c1->cd( 3 );
    hitOcc2DAvg->Draw( "colz" );
  }

  if( hitE2D && hitOcc2D ){

    TH2F* hitE2DAvg = (TH2F*)hitE2D->Clone( "hitE2DAvg" );
    hitE2DAvg->Divide( hitOcc2D );
    hitE2DAvg->SetStats( 0 );
    hitE2DAvg->SetMinimum( 10 );
    c1->cd( 4 );
    hitE2DAvg->Draw( "colz" );
  }
}
