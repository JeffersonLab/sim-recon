// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSC/Hit/LeftArm/Hit_IntegralVsModule_LeftArm
// hnamepath: /PSC/Hit/RightArm/Hit_IntegralVsModule_RightArm
// hnamepath: /PSC/Hit/LeftArm/Hit_fadcTimeVsModule_LeftArm
// hnamepath: /PSC/Hit/RightArm/Hit_fadcTimeVsModule_RightArm

{  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSC");
  if(dir) dir->cd();

  TH2I* hPI_l = (TH2I*)gDirectory->Get("Hit/LeftArm/Hit_IntegralVsModule_LeftArm");
  TH2I* hPI_r = (TH2I*)gDirectory->Get("Hit/RightArm/Hit_IntegralVsModule_RightArm");
  TH2I* hTl = (TH2I*)gDirectory->Get("Hit/LeftArm/Hit_fadcTimeVsModule_LeftArm");
  TH2I* hTr = (TH2I*)gDirectory->Get("Hit/RightArm/Hit_fadcTimeVsModule_RightArm");

  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","Coarse PS Hit Monitor II",150,10,990,660);
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if(!gPad) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(2,2);

  double tsize = 0.0475;  
  gStyle->SetOptStat("emr");
  if(hPI_l){
    c1->cd(1);
    hPI_l->SetTitleSize(tsize,"xy");
    hPI_l->GetYaxis()->SetRange(hPI_l->FindFirstBinAbove(10.0,2),hPI_l->FindLastBinAbove(10.0,2));
    hPI_l->Draw("colz");
  }
 
  if(hPI_r){
    c1->cd(2);
    hPI_r->SetTitleSize(tsize,"xy");
    hPI_r->GetYaxis()->SetRange(hPI_r->FindFirstBinAbove(10.0,2),hPI_r->FindLastBinAbove(10.0,2));
    hPI_r->Draw("colz");
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
