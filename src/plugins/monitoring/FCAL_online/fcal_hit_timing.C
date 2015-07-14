
// hnamepath:  /fcal/hitOcc2D
// hnamepath:  /fcal/hitTmT02D
// hnamepath:  /fcal/hitTmT0Sq2D
// hnamepath:  /fcal/hitTmT0
// hnamepath:  /fcal/hitT0

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("fcal");
  if(dir) dir->cd();

  TH2F* hitOcc2D    = (TH2F*)gDirectory->FindObjectAny("hitOcc2D");
  TH2F* hitTmT02D   = (TH2F*)gDirectory->FindObjectAny("hitTmT02D");
  TH2F* hitTmT0Sq2D = (TH2F*)gDirectory->FindObjectAny("hitTmT0Sq2D");
  TH1I* hitTmT0     = (TH1I*)gDirectory->FindObjectAny("hitTmT0");
  TH1I* hitT0       = (TH1I*)gDirectory->FindObjectAny("hitT0");
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "FCAL Monitor", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide( 2, 2 );
  
  if( hitT0 ){

    hitT0->SetStats( 0 );
    hitT0->SetFillColor( kYellow );
    c1->cd( 1 );
    hitT0->Draw();
  }

  if( hitTmT0 ){

    hitTmT0->SetStats( 0 );
    hitTmT0->SetFillColor( kYellow );
    c1->cd( 2 );
    hitTmT0->Draw();
  }

  if( hitTmT02D && hitOcc2D ){

    TH2F* hitTmT02DAvg = (TH2F*)hitTmT02D->Clone( "hitTmT02DAvg" );
    hitTmT02DAvg->Divide( hitOcc2D );

    for( int x = 1; x <= hitTmT02DAvg->GetNbinsX(); ++x ){
      for( int y = 1; y <= hitTmT02DAvg->GetNbinsY(); ++y ){

	if( hitOcc2D->GetBinContent( x, y ) == 0 ){

	  // set this off scale low so unused blocks
	  // appear white in the plot
	  hitTmT02DAvg->SetBinContent( x, y, -1E3 );
	}
      }
    }

    hitTmT02DAvg->SetMinimum(-10 );
    hitTmT02DAvg->SetMaximum( 10 );
    hitTmT02DAvg->SetStats( 0 );
    c1->cd( 3 );
    hitTmT02DAvg->Draw( "colz" );
  }

  if( hitTmT02D && hitTmT0Sq2D && hitOcc2D ){

    TH2F* hitTmT02DAvg = (TH2F*)hitTmT02D->Clone( "hitTmT02DAvg" );
    TH2F* hitTmT02DRMS = (TH2F*)hitTmT02D->Clone( "hitTmT02DRMS" );
    TH2F* hitTmT0Sq2DAvg = (TH2F*)hitTmT0Sq2D->Clone( "hitTmT0Sq2DAvg" );
    hitTmT02DAvg->Divide( hitOcc2D );
    hitTmT0Sq2DAvg->Divide( hitOcc2D );

    hitTmT02DRMS->SetTitle( "FCAL Local Hit Time RMS [ns]" );

    for( int x = 1; x <= hitTmT02DAvg->GetNbinsX(); ++x ){
      for( int y = 1; y <= hitTmT02DAvg->GetNbinsY(); ++y ){

	if( hitOcc2D->GetBinContent( x, y ) != 0 ){

	  double var = hitTmT0Sq2DAvg->GetBinContent( x, y );
	  var -= ( hitTmT02DAvg->GetBinContent( x, y ) *
		   hitTmT02DAvg->GetBinContent( x, y ) );

	  hitTmT02DRMS->SetBinContent( x, y, TMath::Sqrt( var ) );
	}
	else{

	  // set this off scale low so unused blocks
	  // appear white in the plot
	  hitTmT02DRMS->SetBinContent( x, y, -1 );
	}
      }
    }

    hitTmT02DRMS->SetMinimum(  0 );
    hitTmT02DRMS->SetMaximum(  5 );
    hitTmT02DRMS->SetStats( 0 );
    c1->cd( 4 );
    hitTmT02DRMS->Draw( "colz" );
  }
}
