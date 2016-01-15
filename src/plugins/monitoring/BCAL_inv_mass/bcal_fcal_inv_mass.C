// hnamepath:  /bcal_inv_mass/bcal_fcal_diphoton_mass_300
// hnamepath:  /bcal_inv_mass/bcal_fcal_diphoton_mass_500
// hnamepath:  /bcal_inv_mass/bcal_fcal_diphoton_mass_700
// hnamepath:  /bcal_inv_mass/bcal_fcal_diphoton_mass_900

{
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("bcal_inv_mass");
  if(dir) dir->cd();

  TH1F* bcal_fcal_diphoton_mass_300 = (TH1F*)gDirectory->FindObjectAny("bcal_fcal_diphoton_mass_300");
  TH1F* bcal_fcal_diphoton_mass_500 = (TH1F*)gDirectory->FindObjectAny("bcal_fcal_diphoton_mass_500");
  TH1F* bcal_fcal_diphoton_mass_700 = (TH1F*)gDirectory->FindObjectAny("bcal_fcal_diphoton_mass_700");
  TH1F* bcal_fcal_diphoton_mass_900 = (TH1F*)gDirectory->FindObjectAny("bcal_fcal_diphoton_mass_900");

  int polnumber = 3;
  float fit_low = 0.04;
  float fit_high = 0.20;
  float par_300[15];
  float par_500[15];
  float par_700[15];
  float par_900[15];

  if(gPad == NULL){

    TCanvas *c1 = new TCanvas( "c1", "BCAL_FCAL_inv_mass_plot", 800, 800 );
    c1->cd(0);
    c1->Draw();
    c1->Update();
  }

  if( !gPad ) return;
  TCanvas* c1 = gPad->GetCanvas();
  c1->Divide(2,2);

  if( bcal_fcal_diphoton_mass_300 ){

    bcal_fcal_diphoton_mass_300->SetStats(0);
    double max_300 = bcal_diphoton_mass_300->GetMaximum();
    c1->cd(1);
    bcal_fcal_diphoton_mass_300->Draw();
    bcal_fcal_diphoton_mass_300->SetLineWidth(2);
    TF1 *fitfunc_300 = new TF1("fitfunc_300",Form("gaus(0)+pol%i(3)",polnumber),fit_low,fit_high);
    fitfunc_300->SetParameters(max_300/2, 0.1, 0.009);
    for (int i=0; i<=polnumber; i++) {
         fitfunc_300->SetParameter(3+i,par_300[i]);
    }
    fitfunc_300->SetParNames("height","mean","sigma");
    fitfunc_300->SetParLimits(0,max_300/5.,max_300+max_300/10.);
    fitfunc_300->SetParLimits(1,0.06,0.18);
    fitfunc_300->SetParLimits(2,0.006,0.03);
    fitfunc_300->SetLineWidth(2);
    bcal_fcal_diphoton_mass_300->Fit(fitfunc_300,"RQ");
    for (int i=0; i<9; i++) {
         par_300[i] = fitfunc_300->GetParameter(i);
    }
    fitfunc_300->Draw("same");

    TPaveText *pt_300 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_300->SetFillColor(0);
    pt_300->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_300[1]*1000));
    pt_300->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_300[2]*1000));
    pt_300->AddText(Form("#sigma/M = %.3f %%",(par_300[2]/par_300[1])*100));
    pt_300->Draw();
  }
  if( bcal_fcal_diphoton_mass_500 ){

    bcal_fcal_diphoton_mass_500->SetStats(0);
    double max_500 = bcal_fcal_diphoton_mass_500->GetMaximum();
    c1->cd(2);
    bcal_fcal_diphoton_mass_500->Draw();
    bcal_fcal_diphoton_mass_500->SetLineWidth(2);
    TF1 *fitfunc_500 = new TF1("fitfunc_500",Form("gaus(0)+pol%i(3)",polnumber),fit_low,fit_high);
    fitfunc_500->SetParameters(max_500/5, 0.135, 0.01);
    for (int i=0; i<=polnumber; i++) {
         fitfunc_500->SetParameter(3+i,par_500[i]);
    }
    fitfunc_500->SetParNames("height","mean","sigma");
    fitfunc_500->SetParLimits(0,0.,max_500+max_500/10.);
    fitfunc_500->SetParLimits(1,0.02,0.30);
    fitfunc_500->SetParLimits(2,0.001,0.050);
    fitfunc_500->SetLineWidth(2);
    bcal_fcal_diphoton_mass_500->Fit(fitfunc_500,"RQ");
    for (int i=0; i<9; i++) {
         par_500[i] = fitfunc_500->GetParameter(i);
    }
    fitfunc_500->Draw("same");

    TPaveText *pt_500 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_500->SetFillColor(0);
    pt_500->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_500[1]*1000));
    pt_500->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_500[2]*1000));
    pt_500->AddText(Form("#sigma/M = %.3f %%",(par_500[2]/par_500[1])*100));
    pt_500->Draw();    

  }
   if( bcal_fcal_diphoton_mass_700 ){
   
    bcal_fcal_diphoton_mass_700->SetStats(0);
    double max_700 = bcal_fcal_diphoton_mass_700->GetMaximum();
    c1->cd(3);
    bcal_fcal_diphoton_mass_700->Draw();
    bcal_fcal_diphoton_mass_700->SetLineWidth(2);
    TF1 *fitfunc_700 = new TF1("fitfunc_700",Form("gaus(0)+pol%i(3)",polnumber),fit_low,fit_high);
    fitfunc_700->SetParameters(max_700/5, 0.135, 0.01);
    for (int i=0; i<=polnumber; i++) {
         fitfunc_700->SetParameter(3+i,par_700[i]);
    }
    fitfunc_700->SetParNames("height","mean","sigma");
    fitfunc_700->SetParLimits(0,0.,max_700+max_700/10.);
    fitfunc_700->SetParLimits(1,0.02,0.30);
    fitfunc_700->SetParLimits(2,0.001,0.050);
    fitfunc_700->SetLineWidth(2);
    bcal_fcal_diphoton_mass_700->Fit(fitfunc_700,"RQ");
    for (int i=0; i<9; i++) {
         par_700[i] = fitfunc_700->GetParameter(i);
    }
    fitfunc_700->Draw("same");

    TPaveText *pt_700 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_700->SetFillColor(0);
    pt_700->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_700[1]*1000));
    pt_700->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_700[2]*1000));
    pt_700->AddText(Form("#sigma/M = %.3f %%",(par_700[2]/par_700[1])*100));
    pt_700->Draw();
  }
   if( bcal_fcal_diphoton_mass_900 ){

    bcal_fcal_diphoton_mass_900->SetStats(0);
    double max_900 = bcal_fcal_diphoton_mass_900->GetMaximum();
    c1->cd(4);
    bcal_fcal_diphoton_mass_900->Draw();
    bcal_fcal_diphoton_mass_900->SetLineWidth(2);
    TF1 *fitfunc_900 = new TF1("fitfunc_900",Form("gaus(0)+pol%i(3)",polnumber),fit_low,fit_high);
    fitfunc_900->SetParameters(max_900/5, 0.135, 0.01);
    for (int i=0; i<=polnumber; i++) {
         fitfunc_900->SetParameter(3+i,par_900[i]);
    }
    fitfunc_900->SetParNames("height","mean","sigma");
    fitfunc_900->SetParLimits(0,0.,max_900+max_900/10.);
    fitfunc_900->SetParLimits(1,0.02,0.30);
    fitfunc_900->SetParLimits(2,0.001,0.050);
    fitfunc_900->SetLineWidth(2);
    bcal_fcal_diphoton_mass_900->Fit(fitfunc_900,"RQ");
    for (int i=0; i<9; i++) {
         par_900[i] = fitfunc_900->GetParameter(i);
    }
    fitfunc_900->Draw("same");

    TPaveText *pt_900 = new TPaveText(0.6, 0.65, 0.99, 0.89, "NDC");
    pt_900->SetFillColor(0);
    pt_900->AddText(Form("M_{#pi^{0}} = %.3f MeV",par_900[1]*1000));
    pt_900->AddText(Form("#sigma_{#pi^{0}} = %.3f MeV",par_900[2]*1000));
    pt_900->AddText(Form("#sigma/M = %.3f %%",(par_900[2]/par_900[1])*100));
    pt_900->Draw();
  }


}
