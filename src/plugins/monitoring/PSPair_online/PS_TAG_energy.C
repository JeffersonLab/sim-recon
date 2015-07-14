// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//

// hnamepath: /PSPair/PSC_PS_TAGH/PSTAGH_EdiffVsEtagh
// hnamepath: /PSPair/PSC_PS_TAGH/PSTAGH_EVsEtagh
// hnamepath: /PSPair/PSC_PS_TAGM/PSTAGM_EdiffVsEtagm
// hnamepath: /PSPair/PSC_PS_TAGM/PSTAGM_EVsEtagm

{  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("PSPair");
  if(dir) dir->cd();

  TH2I* hEdiff_tagh = (TH2I*)gDirectory->Get("PSC_PS_TAGH/PSTAGH_EdiffVsEtagh");
  TH2I* hE_tagh = (TH2I*)gDirectory->Get("PSC_PS_TAGH/PSTAGH_EVsEtagh");
  TH2I* hEdiff_tagm = (TH2I*)gDirectory->Get("PSC_PS_TAGM/PSTAGM_EdiffVsEtagm");
  TH2I* hE_tagm = (TH2I*)gDirectory->Get("PSC_PS_TAGM/PSTAGM_EVsEtagm");

  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","PS/Tagger Energy Correlation Monitor",150,10,990,660);
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if(!gPad) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(2,2);

  double tsize = 0.0475;  
  gStyle->SetOptStat("emr");
  if(hE_tagh){
    c1->cd(2);
    hE_tagh->SetTitleSize(tsize,"xy");
    hE_tagh->GetXaxis()->SetRange(hE_tagh->FindFirstBinAbove(10.0),hE_tagh->FindLastBinAbove(10.0));
    hE_tagh->GetYaxis()->SetRange(hE_tagh->FindFirstBinAbove(10.0,2),hE_tagh->FindLastBinAbove(10.0,2));
    hE_tagh->Draw("colz");
  }
 
  if(hE_tagm){
    c1->cd(1);
    hE_tagm->SetTitleSize(tsize,"xy");
    hE_tagm->GetXaxis()->SetRange(hE_tagm->FindFirstBinAbove(10.0),hE_tagm->FindLastBinAbove(10.0));
    hE_tagm->GetYaxis()->SetRange(hE_tagm->FindFirstBinAbove(10.0,2),hE_tagm->FindLastBinAbove(10.0,2));
    hE_tagm->Draw("colz");
  }
 
  if(hEdiff_tagh){
    c1->cd(4);
    hEdiff_tagh->SetTitleSize(tsize,"xy");
    hEdiff_tagh->GetXaxis()->SetRange(hEdiff_tagh->FindFirstBinAbove(10.0),hEdiff_tagh->FindLastBinAbove(10.0));
    hEdiff_tagh->GetYaxis()->SetRange(hEdiff_tagh->FindFirstBinAbove(10.0,2),hEdiff_tagh->FindLastBinAbove(10.0,2));
    hEdiff_tagh->Draw("colz");
  }

  if(hEdiff_tagm){
    c1->cd(3);
    hEdiff_tagm->SetTitleSize(tsize,"xy");
    hEdiff_tagm->GetXaxis()->SetRange(hEdiff_tagm->FindFirstBinAbove(10.0),hEdiff_tagm->FindLastBinAbove(10.0));
    hEdiff_tagm->GetYaxis()->SetRange(hEdiff_tagm->FindFirstBinAbove(10.0,2),hEdiff_tagm->FindLastBinAbove(10.0,2));
    hEdiff_tagm->Draw("colz");
  }

}
