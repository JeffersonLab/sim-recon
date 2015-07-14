
// hnamepath:  /fcal/digOcc2D
// hnamepath:  /fcal/digTmT02D
// hnamepath:  /fcal/digTmT0
// hnamepath:  /fcal/digT0
// hnamepath:  /fcal/digT

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("fcal");
  if(dir) dir->cd();

  TH2F *digOcc2D   = (TH2F*)gDirectory->FindObjectAny("digOcc2D");
  TH2F *digTmT02D  = (TH2F*)gDirectory->FindObjectAny("digTmT02D");
  TH1I *digTmT0    = (TH1I*)gDirectory->FindObjectAny("digTmT0");
  TH1I *digT0      = (TH1I*)gDirectory->FindObjectAny("digT0");
  TH1I *digT       = (TH1I*)gDirectory->FindObjectAny("digT");
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL Monitor", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide( 2, 2 );
  
  if( digT ){

    digT->SetStats( 0 );
    digT->SetFillColor( kRed );
    c1->cd( 1 );
    digT->Draw();
  }

  if( digT0 ){

    digT0->SetStats( 0 );
    digT0->SetFillColor( kRed );
    c1->cd( 2 );
    digT0->Draw();
  }

  if( digTmT0 ){

    digTmT0->SetStats( 0 );
    digTmT0->SetFillColor( kRed );
    c1->cd( 3 );
    digTmT0->Draw();
  }

  if( digTmT02D && digOcc2D ){

    TH2F* digTmT02DAvg = (TH2F*)digTmT02D->Clone( "digTmT02DAvg" );
    digTmT02DAvg->Divide( digOcc2D );

    for( int x = 1; x <= digTmT02DAvg->GetNbinsX(); ++x ){
      for( int y = 1; y <= digTmT02DAvg->GetNbinsY(); ++y ){

	if( digOcc2D->GetBinContent( x, y ) == 0 ){

	  // set this off scale low so unused blocks
	  // appear white in the plot
	  digTmT02DAvg->SetBinContent( x, y, -1E3 );
	}
      }
    }

    digTmT02DAvg->SetMinimum(-200 );
    digTmT02DAvg->SetMaximum( 200 );
    digTmT02DAvg->SetStats( 0 );
    c1->cd( 4 );
    digTmT02DAvg->Draw( "colz" );
  }
}
