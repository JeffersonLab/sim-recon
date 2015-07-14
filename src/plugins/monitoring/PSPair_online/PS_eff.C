// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSPair/PSC/PSC_PSCIDLeftVsIDRight
// hnamepath: /PSPair/PSC_PS/PS_PSCIDLeftVsIDRight

{
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSPair");
  if(dir) dir->cd();
  
  TH2F *id_LR_D =(TH2F*)gDirectory->Get("PSC/PSC_PSCIDLeftVsIDRight");
  TH2F *id_LR_N =(TH2F*)gDirectory->Get("PSC_PS/PS_PSCIDLeftVsIDRight");
    
  TH2F *eff = (TH2F*)id_LR_N->Clone(); 
  eff->Sumw2();
  eff->Divide(id_LR_N,id_LR_D);

  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","Fine PS Efficiency",150,10,990,660);
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }
    
  if(!gPad) return;
  TCanvas* c1 = gPad->GetCanvas();

  double tsize = 0.045;
  gStyle->SetOptStat("");
  eff->SetTitle("Fine PS Efficiency: N(PSC,PS) / N(PSC)");
  eff->SetTitleSize(tsize);
  eff->SetTitleSize(tsize,"xy");
  eff->Draw("colz");
    
}
