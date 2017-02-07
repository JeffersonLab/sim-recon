// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_300
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_500
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_700
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_900
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_v_E
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_v_z_lowE
// hnamepath:  /bcal_inv_mass/bcal_diphoton_mass_v_z_highE

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal_inv_mass");
  if(dir) dir->cd();

  TH1I* bcal_diphoton_mass_300 = (TH1I*)gDirectory->FindObjectAny("bcal_diphoton_mass_300");
  TH1I* bcal_diphoton_mass_500 = (TH1I*)gDirectory->FindObjectAny("bcal_diphoton_mass_500");
  TH1I* bcal_diphoton_mass_700 = (TH1I*)gDirectory->FindObjectAny("bcal_diphoton_mass_700");
  TH1I* bcal_diphoton_mass_900 = (TH1I*)gDirectory->FindObjectAny("bcal_diphoton_mass_900");
  TH2I* bcal_diphoton_mass_v_E = (TH2I*)gDirectory->FindObjectAny("bcal_diphoton_mass_v_E");
  TH2I* bcal_diphoton_mass_v_z_lowE = (TH2I*)gDirectory->FindObjectAny("bcal_diphoton_mass_v_z_lowE");
  TH2I* bcal_diphoton_mass_v_z_highE = (TH2I*)gDirectory->FindObjectAny("bcal_diphoton_mass_v_z_highE");

  int polnumber = 3;
  float fit_low = 0.07;
  float fit_high = 0.20;
  float par_300[15]; 
  float par_500[15];
  float par_700[15];
  float par_900[15];

  TCanvas *c1 = NULL;
  if(gPad == NULL){
    c1 = new TCanvas( "c1", "BCAL_inv_mass_plot", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  } else {
      c1 = gPad->GetCanvas();
  }

  if( !gPad ) return;
  c1->Divide(2,2);

   if( bcal_diphoton_mass_v_E ){
    c2->cd(1)->SetLogz(1);
    bcal_diphoton_mass_v_E->Draw("colz");
   }
	
   if( bcal_diphoton_mass_v_z_lowE ){
    c2->cd(2)->SetLogz(1);
    bcal_diphoton_mass_v_z_lowE->Draw("colz");
   }

   if( bcal_diphoton_mass_v_z_highE ){
    c2->cd(3)->SetLogz(1);
    bcal_diphoton_mass_v_z_highE->Draw("colz");
   }

}
