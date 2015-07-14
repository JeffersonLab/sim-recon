// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSC/Hit/LeftArm/Hit_HasTDCvsHasADC_LeftArm
// hnamepath: /PSC/Hit/RightArm/Hit_HasTDCvsHasADC_RightArm
// hnamepath: /PSC/Hit/LeftArm/Hit_tdcadcTimeDiffVsModule_LeftArm
// hnamepath: /PSC/Hit/RightArm/Hit_tdcadcTimeDiffVsModule_RightArm

{  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSC");
  if(dir) dir->cd();

  TH2F* hStat_l = (TH2F*)gDirectory->Get("Hit/LeftArm/Hit_HasTDCvsHasADC_LeftArm");
  TH2F* hStat_r = (TH2F*)gDirectory->Get("Hit/RightArm/Hit_HasTDCvsHasADC_RightArm");
  TH2I* hTl = (TH2I*)gDirectory->Get("Hit/LeftArm/Hit_tdcadcTimeDiffVsModule_LeftArm");
  TH2I* hTr = (TH2I*)gDirectory->Get("Hit/RightArm/Hit_tdcadcTimeDiffVsModule_RightArm");

  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","Coarse PS Hit Monitor III",150,10,990,660);
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if(!gPad) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(2,2);

  double tsize = 0.0475;  
  gStyle->SetOptStat("emr");
  if(hStat_l){
    c1->cd(1);
    hStat_l->SetTitleSize(tsize,"xy");
    hStat_l->GetYaxis()->SetRange(hStat_l->FindFirstBinAbove(10.0,2),hStat_l->FindLastBinAbove(10.0,2));
    hStat_l->SetMarkerColor(kRed);
    hStat_l->SetMarkerSize(2.0);
    hStat_l->DrawNormalized("text");
  }
 
  if(hStat_r){
    c1->cd(2);
    hStat_r->SetTitleSize(tsize,"xy");
    hStat_r->GetYaxis()->SetRange(hStat_r->FindFirstBinAbove(10.0,2),hStat_r->FindLastBinAbove(10.0,2));
    hStat_r->SetMarkerColor(kRed);
    hStat_r->SetMarkerSize(2.0);
    hStat_r->DrawNormalized("text");
  }
 
  if(hTl){
    c1->cd(3);
    hTl->SetTitleSize(tsize,"xy");
    hTl->GetYaxis()->SetRange(hTl->FindFirstBinAbove(10.0,2),hTl->FindLastBinAbove(10.0,2));
    hTl->Draw("colz");
  }

  if(hTr){
    c1->cd(4);
    hTr->SetTitleSize(tsize,"xy");
    hTr->GetYaxis()->SetRange(hTr->FindFirstBinAbove(10.0,2),hTr->FindLastBinAbove(10.0,2));
    hTr->Draw("colz");
  }

}
