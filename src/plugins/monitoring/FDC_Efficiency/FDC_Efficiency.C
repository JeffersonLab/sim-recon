
void FDC_Efficiency(bool save = 0){
  
  gStyle->SetPalette(1);

  //gDirectory->cd();
  //TDirectory *gDirectory = TFile::Open("hd_root.root");
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("FDC_Efficiency");
  if(!dir) return;
  dir->cd();

  gDirectory->cd("Track_Quality");
  TCanvas *cTrackQuality = new TCanvas("cTrackQuality", "Track Quality Histograms", 900, 800);
  cTrackQuality->Divide(3,2);

  cTrackQuality->cd(1);
  TH1 *tmom = (TH1*) gDirectory->Get("hMom");
  TH1 *tmom_acc = (TH1*) gDirectory->Get("hMom_accepted");
  tmom->Draw();
  tmom_acc->SetLineColor(2);
  tmom_acc->Draw("same");

  cTrackQuality->cd(2);
  TH1 *ttheta = (TH1*) gDirectory->Get("hTheta");
  TH1 *ttheta_acc = (TH1*) gDirectory->Get("hTheta_accepted");
  ttheta->Draw();
  ttheta_acc->SetLineColor(2);
  ttheta_acc->Draw("same");

  cTrackQuality->cd(3);
  TH1 *tphi = (TH1*) gDirectory->Get("hPhi");
  TH1 *tphi_acc = (TH1*) gDirectory->Get("hPhi_accepted");
  tphi->SetMinimum(0.0);
  tphi->Draw();
  tphi_acc->SetLineColor(2);
  tphi_acc->Draw("same");

  cTrackQuality->cd(4);
  TH1 *tchi = (TH1*) gDirectory->Get("hChi2OverNDF");
  TH1 *tchi_acc = (TH1*) gDirectory->Get("hChi2OverNDF_accepted");
  tchi->Draw();
  tchi_acc->SetLineColor(2);
  tchi_acc->Draw("same");

  cTrackQuality->cd(5);
  TH1 *tcells = (TH1*) gDirectory->Get("hCellsHit");
  TH1 *tcells_acc = (TH1*) gDirectory->Get("hCellsHit_accepted");
  tcells->Draw();
  tcells_acc->SetLineColor(2);
  tcells_acc->Draw("same");

  cTrackQuality->cd(6);
  TH1 *trings = (TH1*) gDirectory->Get("hRingsHit");
  TH1 *trings_acc = (TH1*) gDirectory->Get("hRingsHit_accepted");
  trings->Draw();
  trings_acc->SetLineColor(2);
  trings_acc->Draw("same");

  Float_t minScale = 0.8;
  Float_t maxScale = 1.0;    

  gDirectory->cd("../FDC_View");
  TCanvas *cEfficiency_Wire = new TCanvas("cEfficiency_Wire", "Wire Efficiency", 1200, 800);
  cEfficiency_Wire->Divide(6,4);

  TGraphAsymmErrors *EffWire[24];

  for(unsigned int icell=1; icell<=24; icell++){
    cEfficiency_Wire->cd(icell);
    char hname1[256];
    sprintf(hname1, "fdc_wire_measured_cell[%d]", icell);
    TH1 *h1 = (TH1*)(gDirectory->Get(hname1));
    char hname2[256];
    sprintf(hname2, "fdc_wire_expected_cell[%d]", icell);
    TH1 *h2 = (TH1*)(gDirectory->Get(hname2));
      
    if(h1 && h2){
      EffWire[icell-1] = new  TGraphAsymmErrors(h1, h2, "ac");
      EffWire[icell-1]->Draw("ap");
      EffWire[icell-1]->SetMinimum(minScale);
      EffWire[icell-1]->SetMaximum(maxScale);
      EffWire[icell-1]->SetTitle("");
      EffWire[icell-1]->GetXaxis()->SetTitle("Wire Number");
      EffWire[icell-1]->GetYaxis()->SetTitle("Efficiency");
      EffWire[icell-1]->GetYaxis()->SetTitleOffset(1.15);
      EffWire[icell-1]->SetMinimum(minScale);
      EffWire[icell-1]->SetMaximum(maxScale);
    }
  }
  if (save) cEfficiency_Wire->SaveAs("cEfficiencyWire.pdf");

  TCanvas *cEfficiency_Pseudo = new TCanvas("cEfficiency_Pseudo", "Pseudo Hit Efficiency", 1200, 800);
  cEfficiency_Pseudo->Divide(6,4);
  double eff[24];
  double tot = 0;
    
  for(unsigned int icell=1; icell<=24; icell++){
    cEfficiency_Pseudo->cd(icell);
    char hname3[256];
    sprintf(hname3, "fdc_pseudo_measured_cell[%d]", icell);
    TH2 *h3 = (TH2*)(gDirectory->Get(hname3));
    char hname4[256];
    sprintf(hname4, "fdc_pseudo_expected_cell[%d]", icell);
    TH2 *h4 = (TH2*)(gDirectory->Get(hname4));
      
    if(h3 && h4){
      eff[icell-1] = (double) h3->GetEntries();
      eff[icell-1] /= (double) h4->GetEntries();
      
      h3->Divide(h4);
      h3->SetMinimum(minScale);
      h3->SetMaximum(maxScale);
      h3->GetXaxis()->SetTitle("X Position (cm)");
      h3->GetYaxis()->SetTitle("Y Position (cm)");
      h3->SetStats(0);
      h3->Draw("colz");

      cout << "Pseudo-Efficiency of cell " << icell << " : " << eff[icell-1] << endl;
      tot += eff[icell-1];
    }
  }

  cout << "Total : " << tot/24. << endl;

  dir->cd();

  TCanvas *cEfficiency_vs = new TCanvas("cEfficiency_vs", "Efficiency vs Things", 900, 600);
  cEfficiency_vs->Divide(3,2);

  cEfficiency_vs->cd(1);
  TH1I *MeasDOCA = (TH1I*)(gDirectory->Get("Offline/Measured Hits Vs DOCA"));
  TH1I *ExpDOCA = (TH1I*)(gDirectory->Get("Offline/Expected Hits Vs DOCA"));
  if(MeasDOCA && ExpDOCA){
    TGraphAsymmErrors *EffDOCA = new  TGraphAsymmErrors(MeasDOCA, ExpDOCA, "ac");
    EffDOCA->Draw("ap");
    EffDOCA->SetMinimum(minScale);
    EffDOCA->SetMaximum(maxScale);
    EffDOCA->SetTitle("FDC Per Wire Efficiency Vs. DOCA");
    EffDOCA->GetXaxis()->SetTitle("Closest distance between track and wire [cm]");
    EffDOCA->GetYaxis()->SetTitle("Efficiency");
    TLine *lcut = new TLine(0.5,minScale,0.5,maxScale);
    lcut->SetLineColor(2);
    lcut->SetLineStyle(2);
    lcut->SetLineWidth(2);
    lcut->Draw();
  }

  cEfficiency_vs->cd(2);
  TH1I *MeasTrackingFOM = (TH1I*)(gDirectory->Get("Offline/Measured Hits Vs Tracking FOM"));
  TH1I *ExpTrackingFOM = (TH1I*)(gDirectory->Get("Offline/Expected Hits Vs Tracking FOM"));
  if(MeasTrackingFOM && ExpTrackingFOM){
    TGraphAsymmErrors *EffTrackingFOM = new  TGraphAsymmErrors(MeasTrackingFOM, ExpTrackingFOM, "ac");
    EffTrackingFOM->Draw("ap");
    EffTrackingFOM->SetMinimum(minScale);
    EffTrackingFOM->SetMaximum(maxScale);
    EffTrackingFOM->SetTitle("FDC Per Wire Efficiency Vs. Tracking FOM");
    EffTrackingFOM->GetXaxis()->SetTitle("Tracking FOM");
    EffTrackingFOM->GetYaxis()->SetTitle("Efficiency");
  }

  cEfficiency_vs->cd(3);
  TH1I *Meastheta = (TH1I*)(gDirectory->Get("Offline/Measured Hits Vs theta"));
  TH1I *Exptheta = (TH1I*)(gDirectory->Get("Offline/Expected Hits Vs theta"));
  if(Meastheta && Exptheta){
    TGraphAsymmErrors *Efftheta = new  TGraphAsymmErrors(Meastheta, Exptheta, "ac");
    Efftheta->Draw("ap");
    Efftheta->SetMinimum(minScale);
    Efftheta->SetMaximum(maxScale);
    Efftheta->SetTitle("FDC Per Wire Efficiency Vs. #theta");
    Efftheta->GetXaxis()->SetTitle("Track #theta [deg.]");
    Efftheta->GetYaxis()->SetTitle("Efficiency");
  }

  cEfficiency_vs->cd(4);
  TH1I *Measphi = (TH1I*)(gDirectory->Get("Offline/Measured Hits Vs phi"));
  TH1I *Expphi = (TH1I*)(gDirectory->Get("Offline/Expected Hits Vs phi"));
  if(Measphi && Expphi){
    TGraphAsymmErrors *Effphi = new  TGraphAsymmErrors(Measphi, Expphi, "ac");
    Effphi->Draw("ap");
    Effphi->SetMinimum(minScale);
    Effphi->SetMaximum(maxScale);
    Effphi->SetTitle("FDC Per Wire Efficiency Vs. #phi");
    Effphi->GetXaxis()->SetTitle("Track #phi [deg.]");
    Effphi->GetYaxis()->SetTitle("Efficiency");
  }

  cEfficiency_vs->cd(5);
  TH1I *Measp = (TH1I*)(gDirectory->Get("Offline/Measured Hits Vs p"));
  TH1I *Expp = (TH1I*)(gDirectory->Get("Offline/Expected Hits Vs p"));
  if(Measp && Expp){
    TGraphAsymmErrors *Effp = new  TGraphAsymmErrors(Measp, Expp, "ac");
    Effp->Draw("ap");
    Effp->SetMinimum(minScale);
    Effp->SetMaximum(maxScale);
    Effp->SetTitle("FDC Per Wire Efficiency Vs. p");
    Effp->GetXaxis()->SetTitle("Track Momentum [GeV/c]");
    Effp->GetYaxis()->SetTitle("Efficiency");
  }

  cEfficiency_vs->cd(6);
  TH1I *Meascells = (TH1I*)(gDirectory->Get("Offline/Measured Hits Vs Hit Cells"));
  TH1I *Expcells = (TH1I*)(gDirectory->Get("Offline/Expected Hits Vs Hit Cells"));
  if(Meascells && Expcells){
    TGraphAsymmErrors *Effcells = new  TGraphAsymmErrors(Meascells, Expcells, "ac");
    Effcells->Draw("ap");
    Effcells->SetMinimum(minScale);
    Effcells->SetMaximum(maxScale);
    Effcells->SetTitle("FDC Per Wire Efficiency Vs. Contributing Cells");
    Effcells->GetXaxis()->SetTitle("FDC Cells");
    Effcells->GetYaxis()->SetTitle("Efficiency");
  }

  if (save) cEfficiency_vs->SaveAs("cEfficiency_vs.pdf");

  dir->cd();
  gDirectory->cd("Residuals");


  TCanvas *cResidual_Pseudo = new TCanvas("cResidual_Pseudo", "Pseudo Hit Resolution", 1000, 800);
  cResidual_Pseudo->Divide(6,4);

  double meanV[24];
  double meanV_err[24];
  double meanU[24];
  double meanU_err[24];
    
  for(unsigned int icell=1; icell<=24; icell++){
    cResidual_Pseudo->cd(icell);
    char hname5[256];
    sprintf(hname5, "hPseudoResV_cell[%d]", icell);
    TH1 *h5 = (TH1*)(gDirectory->Get(hname5));
    char hname6[256];
    sprintf(hname6, "hPseudoResU_cell[%d]", icell);
    TH1 *h6 = (TH1*)(gDirectory->Get(hname6));
      
    h5->GetXaxis()->SetTitle("Position Resolution (cm)");
    h5->GetXaxis()->SetRangeUser(-1,1);
    h5->Draw();
    h5->Fit("gaus","q0");
    TF1 *fgaus5 = h5->GetFunction("gaus");
    meanV[icell-1] = fgaus5->GetParameter(1);
    //mean_err[icell-1] = fgaus5->GetParError(1) * fgaus->GetChisquare()/(fgaus->GetNDF()-1);
    meanV_err[icell-1] = fgaus5->GetParameter(2);

    h6->SetLineColor(2);
    h6->Draw("same");
    h6->Fit("gaus","q0");
    TF1 *fgaus6 = h6->GetFunction("gaus");
    meanU[icell-1] = fgaus6->GetParameter(1);
    //mean_err[icell-1] = fgaus6->GetParError(1) * fgaus->GetChisquare()/(fgaus->GetNDF()-1);
    meanU_err[icell-1] = fgaus6->GetParameter(2);
    
  }

  const unsigned int rad = 1;

  TCanvas *cResidual_Pseudo2 = new TCanvas("cResidual_Pseudo2", "Pseudo Hit Resolution 2D", 1000, 800);
  cResidual_Pseudo2->Divide(6,4);

  for(unsigned int icell=1; icell<=24; icell++){
    cResidual_Pseudo2->cd(icell);
    for (unsigned int r=0; r<rad; r++){
      char hname7[256];
      sprintf(hname7, "hPseudoResUvsV_cell[%d]_radius[%d]", icell, (r+1)*(45/rad));
      TH2 *h7 = (TH2*)(gDirectory->Get(hname7));
      
      if (h7){
	h7->GetYaxis()->SetTitle("Position reconstructed along Wire (cm)");
	h7->GetYaxis()->SetRangeUser(-0.5,0.5);
	h7->GetXaxis()->SetTitle("Position perp. to Wire (cm)");
	h7->GetXaxis()->SetRangeUser(-0.5,0.5);
	h7->Draw("col");
	h7->ProfileY();
      }
    }
  }

  double slope[rad][24];
  double slope_err[rad][24];
  double cell[rad][24];
  double cell_err[rad][24];
  TCanvas *cResidual_Profile = new TCanvas("cResidual_Profile", "Pseudo Resolution Profile", 1000, 800);
  cResidual_Profile->Divide(6,4);
  
  for(unsigned int icell=1; icell<=24; icell++){
    cResidual_Profile->cd(icell);
    for (unsigned int r=0; r<rad; r++){
      cell[r][icell-1] = icell + 0.1*r;
      cell_err[r][icell-1] = 0;
      char hname8[256];
      sprintf(hname8, "hPseudoResUvsV_cell[%d]_radius[%d]_pfy", icell, (r+1)*(45/rad));
      TH1 *h8 = (TH1*)(gDirectory->Get(hname8));
      
      if (h8){
	h8->GetXaxis()->SetTitle("Position reconstructed along Wire (cm)");
	h8->GetXaxis()->SetRangeUser(-0.15,0.15);
	h8->GetYaxis()->SetTitle("Position perp. to Wire (cm)");
	h8->GetYaxis()->SetRangeUser(-0.5,0.5);
	h8->GetYaxis()->SetTitleFont(42);
	h8->GetYaxis()->SetTitleSize(0.035);
	h8->GetYaxis()->SetLabelFont(42);
	h8->GetYaxis()->SetLabelSize(0.035);
	h8->SetLineColor(r+1);
	if (r == 0)
	  h8->Draw();
	else 
	  h8->Draw("same");
	
	TF1 *fp1 = new TF1("fp1","[0]+x*[1]",-0.08,0.08);
	fp1->SetLineColor(r+1);
	h8->Fit("fp1","qr");
	slope[r][icell-1] = fp1->GetParameter(1);
	slope_err[r][icell-1] = fp1->GetParError(1);
      }
      else {
	// histograms not present
	slope[r][icell-1] = 0;
	slope_err[r][icell-1] = 0;
      }
    }
  }

  TCanvas *cAlignment = new TCanvas("cAlignment", "Alignment", 1200, 1000);
  cAlignment->Divide(2,2);
  
  cAlignment->cd(1);
  TGraphErrors *galignV = new TGraphErrors(24, cell[0], meanV, cell_err[0], meanV_err);
  galignV->SetTitle("Alignment along Wire (Error Bars: Sigma); Cell # ; V_{Hit} - V_{Track} (cm)");
  galignV->SetMarkerColor(1);
  galignV->SetMarkerStyle(2);
  galignV->Draw("AP");

  TLine *lpack12 = new TLine(6.5, -0.1, 6.5, 0.1);
  lpack12->SetLineColor(2);
  lpack12->SetLineStyle(2);
  lpack12->SetLineWidth(2);
  lpack12->Draw();
  TLine *lpack23 = new TLine(12.5, -0.1, 12.5, 0.1);
  lpack23->SetLineColor(2);
  lpack23->SetLineStyle(2);
  lpack23->SetLineWidth(2);
  lpack23->Draw();
  TLine *lpack34 = new TLine(18.5, -0.1, 18.5, 0.1);
  lpack34->SetLineColor(2);
  lpack34->SetLineStyle(2);
  lpack34->SetLineWidth(2);
  lpack34->Draw();

  cAlignment->cd(2);
  TGraphErrors *galignU = new TGraphErrors(24, cell[0], meanU, cell_err[0], meanU_err);
  galignU->SetTitle("Alignment perp. to Wire (Error Bars: Sigma); Cell # ; U_{Hit} - U_{Track} (cm)");
  galignU->SetMarkerColor(1);
  galignU->SetMarkerStyle(2);
  galignU->Draw("AP");

  lpack12->Draw();
  lpack23->Draw();
  lpack34->Draw();

  cAlignment->cd(3);
  TGraphErrors *gmagnet[rad];
  for (unsigned int r=0; r<rad; r++){
    gmagnet[r] = new TGraphErrors(24, cell[r], slope[r], cell_err[r], slope_err[r]);
    gmagnet[r]->SetTitle("Value for slope of magnetic deflection; Cell # ; Slope (cm^{-1})");
    gmagnet[r]->SetMarkerColor(r+1);
    gmagnet[r]->SetMarkerStyle(2);
    if (r==0)
      gmagnet[r]->Draw("AP");
    else 
      gmagnet[r]->Draw("Psame");
  }
}
