// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSC/Hit/LeftArm/Hit_Occupancy_LeftArm
// hnamepath: /PSC/Hit/RightArm/Hit_Occupancy_RightArm
// hnamepath: /PSC/Hit/LeftArm/Hit_TimeVsModule_LeftArm
// hnamepath: /PSC/Hit/RightArm/Hit_TimeVsModule_RightArm

{  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSC");
  if(dir) dir->cd();

  TH1I* hOcc_l = (TH1I*)gDirectory->Get("Hit/LeftArm/Hit_Occupancy_LeftArm");
  TH1I* hOcc_r = (TH1I*)gDirectory->Get("Hit/RightArm/Hit_Occupancy_RightArm");
  TH2I* hTl = (TH2I*)gDirectory->Get("Hit/LeftArm/Hit_TimeVsModule_LeftArm");
  TH2I* hTr = (TH2I*)gDirectory->Get("Hit/RightArm/Hit_TimeVsModule_RightArm");

  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","Coarse PS Hit Monitor I",150,10,990,660);
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if(!gPad) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(2,2);

  double tsize = 0.0475;  
  gStyle->SetOptStat("emr");
  if(hOcc_l){
    hOcc_l->SetFillColor(kBlue);
    c1->cd(1);
    hOcc_l->SetTitleSize(tsize,"xy");
    hOcc_l->Draw();
  }
 
  if(hOcc_r){
    hOcc_r->SetFillColor(kBlue);
    c1->cd(2);
    hOcc_r->SetTitleSize(tsize,"xy");
    hOcc_r->Draw();
  }
 
  if(hTl){
    hTl->SetFillColor(kBlue);
    c1->cd(3);
    hTl->SetTitleSize(tsize,"xy");
    hTl->GetYaxis()->SetRange(hTl->FindFirstBinAbove(10.0,2),hTl->FindLastBinAbove(10.0,2));
    hTl->Draw("colz");
  }

  if(hTr){
    hTr->SetFillColor(kBlue);
    c1->cd(4);
    hTr->SetTitleSize(tsize,"xy");
    hTr->GetYaxis()->SetRange(hTr->FindFirstBinAbove(10.0,2),hTr->FindLastBinAbove(10.0,2));
    hTr->Draw("colz");
  }

}
