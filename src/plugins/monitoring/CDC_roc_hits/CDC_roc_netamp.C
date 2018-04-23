
// The following are special comments used by RootSpy to know
// which histograms to fetch for the macro.
//
// hnamepath: /CDC_roc_hits/cdc_amp_roc25
// hnamepath: /CDC_roc_hits/cdc_amp_roc26
// hnamepath: /CDC_roc_hits/cdc_amp_roc27
// hnamepath: /CDC_roc_hits/cdc_amp_roc28


{

  gStyle->SetOptStat(kFALSE);
  gStyle->SetFillStyle(0);
  

  // Get number of events
  double Nevents = 1.0;
  TDirectory *CDCdir = (TDirectory*)gDirectory->FindObjectAny("CDC_roc_hits");

  if(!CDCdir) return;

  CDCdir->cd();

  TH1I *cdc_nevents = (TH1I*)CDCdir->Get("cdc_nevents");
  if(cdc_nevents) Nevents = (double)cdc_nevents->GetBinContent(1);

  // Just for testing
  if(gPad == NULL){
    TCanvas *c1 = new TCanvas("c1","c1 CDC_roc_hits",1600,800);
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }
	
  if(!gPad) return;
  TCanvas *c1 = gPad->GetCanvas();
  c1->Divide(1,4);

  for(unsigned int iroc=25; iroc<=28; iroc++){
    c1->cd(iroc-24);

    gPad->SetLeftMargin(0.001);
    gPad->SetRightMargin(0.001);
    gPad->SetBottomMargin(0.02);
    gPad->SetTopMargin(0.02);

    char hname[256];
    sprintf(hname, "cdc_netamp_roc%d", iroc);
    TH2 *h = (TH2*)(CDCdir->Get(hname));

    if(h){
      //      sprintf(hname, "cdc_netamp_roc%d_norm", iroc);
      // TH2 *hh = (TH2*)h->Clone(hname);
      TH2 *hh = (TH2*)h->Clone();

      //hh->Scale(1.0/Nevents);
      //      hh->SetStats(0);

      hh->SetTitle("");   
      hh->GetXaxis()->SetTitle("");
      hh->GetYaxis()->SetTitle("");

      hh->GetXaxis()->SetRangeUser(290,1780);
      hh->GetYaxis()->SetRangeUser(0,4096);

      //hh->GetYaxis()->SetNdivisions(210,kFALSE);
      hh->GetYaxis()->SetNdivisions(110);
      hh->GetXaxis()->SetNdivisions(0);
      hh->GetYaxis()->SetTickLength(0.001);

      hh->GetYaxis()->SetLabelSize(0);  //0.05
      hh->GetYaxis()->SetTitleSize(0);
      hh->GetXaxis()->SetLabelSize(0);  //0.05
      hh->GetXaxis()->SetTitleSize(0);

      hh->Draw("col");  // draw remaining histos without overwriting color palette
      c1->Update();  

      int lastslot = 17;
      if (iroc==26 || iroc==27) lastslot = 16;     

      TPaveText *pt = new TPaveText(0.5,0.15,0.7,0.85,"NDC");
      pt->AddText(Form("ROC %i",iroc));
      pt->AddText("");
      pt->AddText("Peak height - pedestal");
      pt->AddText("vs channel");
      pt->AddText("");
      pt->AddText(Form("slots 3 to %i",lastslot));
      pt->Draw();
      pt->SetBorderSize(0);

      gPad->Modified();


    }
  }


}
