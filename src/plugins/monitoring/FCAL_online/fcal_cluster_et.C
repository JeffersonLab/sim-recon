// hnamepath:  /fcal/clusN
// hnamepath:  /fcal/clusE
// hnamepath:  /fcal/clusT0
// hnamepath:  /fcal/clusTmT0

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("fcal");
  if(dir) dir->cd();

  TH1I* clusN    = (TH1I*)gDirectory->FindObjectAny( "clusN" );
  TH1I* clusE    = (TH1I*)gDirectory->FindObjectAny( "clusE" );
  TH1I* clusT0   = (TH1I*)gDirectory->FindObjectAny( "clusT0" );
  TH1I* clusTmT0 = (TH1I*)gDirectory->FindObjectAny( "clusTmT0" );
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL Monitor", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide( 2, 2 );

  if( clusN ){

    clusN->SetStats( 0 );
    clusN->SetFillColor( kBlue );
    c1->cd( 1 );
    clusN->Draw();
  }
  
  if( clusE ){

    clusE->SetStats( 0 );
    clusE->SetFillColor( kBlue );
    c1->cd( 2 );
    clusE->Draw();
  }

  if( clusT0 ){

    clusT0->SetStats( 0 );
    clusT0->SetFillColor( kBlue );
    c1->cd( 3 );
    clusT0->Draw();
  }

  if( clusTmT0 ){

    clusTmT0->SetStats( 0 );
    clusTmT0->SetFillColor( kBlue );
    c1->cd( 4 );
    clusTmT0->Draw();
  }

}
