// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_300
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_500
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_700
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_900

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal_inv_mass");
  if(dir) dir->cd();

  TH1F* bcal_diphoton_mass_300 = (TH1F*)gDirectory->FindObjectAny("bcal_diphoton_mass_300");
  TH1F* bcal_diphoton_mass_500 = (TH1F*)gDirectory->FindObjectAny("bcal_diphoton_mass_500");
  TH1F* bcal_diphoton_mass_700 = (TH1F*)gDirectory->FindObjectAny("bcal_diphoton_mass_700");
  TH1F* bcal_diphoton_mass_900 = (TH1F*)gDirectory->FindObjectAny("bcal_diphoton_mass_900");
 
  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "BCAL_inv_mass_plot", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(2,2);

  if( bcal_diphoton_mass_300 ){

    bcal_diphoton_mass_300->SetStats(0);
    c1->cd(1);
    bcal_diphoton_mass_300->Draw();
  }
  if( bcal_diphoton_mass_500 ){

    bcal_diphoton_mass_500->SetStats(0);
    c1->cd(2);
    bcal_diphoton_mass_500->Draw();
  }
   if( bcal_diphoton_mass_700 ){

    bcal_diphoton_mass_700->SetStats(0);
    c1->cd(3);
    bcal_diphoton_mass_700->Draw();
  }
   if( bcal_diphoton_mass_900 ){

    bcal_diphoton_mass_900->SetStats(0);
    c1->cd(4);
    bcal_diphoton_mass_900->Draw();
  }


}
