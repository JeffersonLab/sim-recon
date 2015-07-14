// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /FDC/Package_1/fdc_pack2_chamber1_wire_occ
// hnamepath: /FDC/Package_1/fdc_pack2_chamber2_wire_occ
// hnamepath: /FDC/Package_1/fdc_pack2_chamber3_wire_occ
// hnamepath: /FDC/Package_1/fdc_pack2_chamber4_wire_occ
// hnamepath: /FDC/Package_1/fdc_pack2_chamber5_wire_occ
// hnamepath: /FDC/Package_1/fdc_pack2_chamber6_wire_occ
{
  gDirectory->cd();
  
  if(gPad == NULL){
      
   TCanvas *c1 = new TCanvas( "cw2", "FDC Monitor", 800, 800 );
   c1->cd(0);
   c1->Draw();
   c1->Update();
  }
  
  if( !gPad ) return;  

  TCanvas *c1=gPad->GetCanvas();
  if (c1){
    c1->cd();
    c1->Divide(3,2);
  
    char myhistoname[80];
    for (unsigned int chamber=1;chamber<=6;chamber++){
      sprintf(myhistoname,"fdc_pack2_chamber%d_wire_occ",chamber);
      TH1I *myhisto=(TH1I *)gDirectory->FindObjectAny(myhistoname);
      c1->cd(chamber);
      if (myhisto) myhisto->Draw();
    }
  }
}
