// hnamepath:  /fcal/clus2GMass
// hnamepath:  /fcal/clusPhi
// hnamepath:  /fcal/clusXYHigh
// hnamepath:  /fcal/clusXYLow

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("fcal");
  if(dir) dir->cd();

  TH1I* clus2GMass = (TH1I*)gDirectory->FindObjectAny( "clus2GMass" );
  TH1I* clusPhi    = (TH1I*)gDirectory->FindObjectAny( "clusPhi" );
  TH1I* clusXYHigh = (TH1I*)gDirectory->FindObjectAny( "clusXYHigh" );
  TH1I* clusXYLow  = (TH1I*)gDirectory->FindObjectAny( "clusXYLow" );
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL Monitor", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide( 2, 2 );

  if( clus2GMass ){

    clus2GMass->SetStats( 0 );
    clus2GMass->SetFillColor( kBlue );
    c1->cd( 1 );
    clus2GMass->Draw();
  }
  
  if( clusPhi ){

    clusPhi->SetStats( 0 );
    clusPhi->SetFillColor( kBlue );
    c1->cd( 2 );
    clusPhi->Draw();
  }

  if( clusXYLow ){

    clusXYLow->SetStats( 0 );
    c1->cd( 3 );
    clusXYLow->Draw( "colz" );
  }

  if( clusXYHigh ){

    clusXYHigh->SetStats( 0 );
    c1->cd( 4 );
    clusXYHigh->Draw( "colz" );
  }

}
