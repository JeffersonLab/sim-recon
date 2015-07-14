
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /CDC_roc_hits/cdc_amp_roc25
// hnamepath: /CDC_roc_hits/cdc_amp_roc26
// hnamepath: /CDC_roc_hits/cdc_amp_roc27
// hnamepath: /CDC_roc_hits/cdc_amp_roc28


{
  // Get number of events
  double Nevents = 1.0;
  TDirectory *CDCdir = (TDirectory*)gDirectory->FindObjectAny("CDC_roc_hits");

  if(!CDCdir) return;

  CDCdir->cd();

  TH1I *cdc_nevents = (TH1I*)CDCdir->Get("cdc_nevents");
  if(cdc_nevents) Nevents = (double)cdc_nevents->GetBinContent(1);

  // Just for testing
  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1");
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }
	
  if(!gPad) return;

  TCanvas *c1 = gPad->GetCanvas();
  c1->Divide(1,4);

  for(unsigned int iroc=25; iroc<=28; iroc++){
    c1->cd(iroc-24);
    char hname[256];
    sprintf(hname, "cdc_amp_roc%d", iroc);
    TH2 *h = (TH2*)(CDCdir->Get(hname));

    if(h){
      sprintf(hname, "cdc_amp_roc%d_norm", iroc);
      TH2 *hh = (TH2*)h->Clone(hname);

      hh->Scale(1.0/Nevents);
      hh->SetStats(0);
      hh->Draw("colz");  // draw remaining histos without overwriting color palette
      c1->Update();  

      hh->GetXaxis()->SetRangeUser(290,1780);
      hh->GetYaxis()->SetRangeUser(100,500);
      hh->GetYaxis()->SetNdivisions(210,kFALSE);
      hh->GetYaxis()->SetTickLength(0.001);

      hh->GetYaxis()->SetLabelSize(0.05);
      hh->GetYaxis()->SetTitleSize(0.05);
      hh->GetXaxis()->SetLabelSize(0.05);
      hh->GetXaxis()->SetTitleSize(0.05);

      gPad->Modified();


    }
  }
}
