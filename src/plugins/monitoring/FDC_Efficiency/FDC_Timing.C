
void FDC_Timing(bool save = 0){
  
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("FDC_Efficiency");
  if(!dir) return;
  dir->cd();

  gDirectory->cd("Residuals");

  TCanvas *cWireTiming = new TCanvas("cWireTiming", "WireTiming", 1000, 800);
  cWireTiming->Divide(6,4);
  
  double cell[24];
  double cell_err[24];
  double wire[24];
  double wire_err[24];

  double cath[24];
  double cath_err[24];
    
  double delta[24];
  double delta_err[24];
    
  double pull[24];
  double pull_err[24];

  for(unsigned int icell=1; icell<=24; icell++){
    cWireTiming->cd(icell);

    cell[icell-1] = icell;
    cell_err[icell-1] = 0;
    
    char hname[256];
    sprintf(hname, "hWireTime_cell[%d]", icell);
    TH1 *hWire = (TH1*)(gDirectory->Get(hname));
      
    hWire->GetXaxis()->SetTitle("Wire Time (ns)");
    //hWire->Draw();

    int tzero_bin = hWire->GetMaximumBin();
    double tzero = hWire->GetXaxis()->GetBinCenter(tzero_bin);

    TF1 *fwire = new TF1("fwire", "gaus(0)", tzero - 10, tzero + 10);
    fwire->SetLineColor(2);
    fwire->SetNpx(600);
    fwire->SetParameter(1,tzero);
    hWire->Fit("fwire","q0r");
    wire[icell-1] = fwire->GetParameter(1);
    wire_err[icell-1] = fwire->GetParError(1);

    sprintf(hname, "hCathodeTime_cell[%d]", icell);
    TH1 *hCathode = (TH1*)(gDirectory->Get(hname));
      
    hCathode->GetXaxis()->SetTitle("Cathode Time (ns)");
    //hCathode->Draw();

    tzero_bin = hCathode->GetMaximumBin();
    tzero = hCathode->GetXaxis()->GetBinCenter(tzero_bin);
    hCathode->GetXaxis()->SetRangeUser(tzero-40,tzero+40);

    TF1 *fcath = new TF1("fcath", "gaus(0)", tzero - 20, tzero + 20);
    fcath->SetLineColor(2);
    fcath->SetNpx(600);
    fcath->SetParameter(1,tzero);
    hCathode->Fit("fcath","qr");
    cath[icell-1] = fcath->GetParameter(1);
    cath_err[icell-1] = fcath->GetParError(1);

    sprintf(hname, "hDeltaTime_cell[%d]", icell);
    TH1 *hDelta = (TH1*)(gDirectory->Get(hname));
      
    hDelta->GetXaxis()->SetTitle("Wire Time - Cathode Time (ns)");
    hDelta->Draw();

    tzero_bin = hDelta->GetMaximumBin();
    tzero = hDelta->GetXaxis()->GetBinCenter(tzero_bin);

    TF1 *fdelta = new TF1("fdelta", "gaus(0)", tzero - 5, tzero + 5);
    fdelta->SetLineColor(2);
    fdelta->SetNpx(600);
    fdelta->SetParameter(1,tzero);
    hDelta->Fit("fdelta","qr");
    delta[icell-1] = fdelta->GetParameter(1);
    delta_err[icell-1] = fdelta->GetParError(1);

    sprintf(hname, "hPullTime_cell[%d]", icell);
    TH1 *hPull = (TH1*)(gDirectory->Get(hname));

    hPull->GetXaxis()->SetTitle("Pseudo Time (from pull) (ns)");
    //hPull->Draw();

    tzero_bin = hPull->GetMaximumBin();
    tzero = hPull->GetXaxis()->GetBinCenter(tzero_bin);
    hPull->GetXaxis()->SetRangeUser(tzero-10,tzero+10);

    TF1 *fpull = new TF1("fpull", "gaus(0)", tzero - 5, tzero + 3);
    fpull->SetLineColor(2);
    fpull->SetNpx(600);
    fpull->SetParameter(1,tzero);
    hPull->Fit("fpull","q0r");
    pull[icell-1] = fpull->GetParameter(1);
    pull_err[icell-1] = fpull->GetParError(1);

  }
  
  TCanvas *cTiming = new TCanvas("cTiming", "Timing", 1400, 1000);
  cTiming->Divide(2,2);
  cTiming->cd(1);
  
  TGraphErrors *gWireTiming = new TGraphErrors(24, cell, wire, cell_err, wire_err);
  gWireTiming->SetTitle("; Cell # ; Wire Timing (ns)");
  gWireTiming->SetMarkerColor(1);
  gWireTiming->SetMarkerStyle(8);
  gWireTiming->Draw("AP");

  cTiming->cd(2);
  TGraphErrors *gCathTiming = new TGraphErrors(24, cell, cath, cell_err, cath_err);
  gCathTiming->SetTitle("; Cell # ; Cathode Timing (ns)");
  gCathTiming->SetMarkerColor(2);
  gCathTiming->SetMarkerStyle(8);
  gCathTiming->Draw("AP");

  cTiming->cd(3);
  TGraphErrors *gDelta = new TGraphErrors(24, cell, delta, cell_err, delta_err);
  gDelta->SetTitle("; Cell # ; Wire - Cathode Timing (ns)");
  gDelta->SetMarkerColor(4);
  gDelta->SetMarkerStyle(8);
  gDelta->Draw("AP");

  cTiming->cd(4);
  TGraphErrors *gPull = new TGraphErrors(24, cell, pull, cell_err, pull_err);
  gPull->SetTitle("; Cell # ; Pseudo Time from Pull (ns)");
  gPull->SetMarkerColor(6);
  gPull->SetMarkerStyle(8);
  gPull->Draw("AP");
}
