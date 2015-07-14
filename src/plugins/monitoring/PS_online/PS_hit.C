// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PS/Hit/LeftArm/Hit_Occupancy_LeftArm
// hnamepath: /PS/Hit/RightArm/Hit_Occupancy_RightArm
// hnamepath: /PS/Hit/LeftArm/Hit_Energy_LeftArm
// hnamepath: /PS/Hit/RightArm/Hit_Energy_RightArm

{  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PS");
  if(dir) dir->cd();

  TH1I* hOcc_l = (TH1I*)gDirectory->Get("Hit/LeftArm/Hit_Occupancy_LeftArm");
  TH1I* hOcc_r = (TH1I*)gDirectory->Get("Hit/RightArm/Hit_Occupancy_RightArm");
  TH1I* hEl = (TH1I*)gDirectory->Get("Hit/LeftArm/Hit_Energy_LeftArm");
  TH1I* hEr = (TH1I*)gDirectory->Get("Hit/RightArm/Hit_Energy_RightArm");

  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","Fine PS Hit Monitor I",150,10,990,660);
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
    c1->cd(3);
    hOcc_l->SetTitleSize(tsize,"xy");
    hOcc_l->Draw();
  }
 
  if(hOcc_r){
    hOcc_r->SetFillColor(kBlue);
    c1->cd(4);
    hOcc_r->SetTitleSize(tsize,"xy");
    hOcc_r->Draw();
  }
 
  if(hEl){
    hEl->SetFillColor(kBlue);
    c1->cd(1);
    hEl->SetTitleSize(tsize,"xy");
    hEl->Draw();
  }

  if(hEr){
    hEr->SetFillColor(kBlue);
    c1->cd(2);
    hEr->SetTitleSize(tsize,"xy");
    hEr->Draw();
  }

}
