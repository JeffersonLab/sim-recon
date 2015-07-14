// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /TPOL/Hit/Hit_NHits

{  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("TPOL");
  if(dir) dir->cd();

  TH1I* tpol_num_events = (TH1I*)gDirectory->FindObjectAny("tpol_num_events");
  TH1I* tpol_hitMultiplicity = (TH1I*)gDirectory->FindObjectAny("tpol_hitMultiplicity");

  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","TPOL Hit Monitor",150,10,990,660);
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if(!gPad) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(1,2);

  double tsize = 0.0475;  
  gStyle->SetOptStat("emr");
  if(tpol_num_events){
    tpol_num_events->SetFillColor(kBlue);
    c1->cd(1);
    tpol_num_events->SetTitleSize(tsize,"xy");
    tpol_num_events->Draw();
  }

  if(tpol_hitMultiplicity){
    tpol_hitMultiplicity->SetFillColor(kBlue);
    c1->cd(2);
    tpol_hitMultiplicity->SetTitleSize(tsize,"xy");
    tpol_hitMultiplicity->GetXaxis()->SetRange(1,tpol_hitMultiplicity->FindLastBinAbove(1.0));
    tpol_hitMultiplicity->Draw();
  }

}
