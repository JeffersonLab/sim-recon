
void FDC_Residuals(bool save = 0){
  
  //gDirectory->cd();
  //TDirectory *gDirectory = TFile::Open("hd_root.root");
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("FDC_Residuals");
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
	h7->ProfileX();
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
      sprintf(hname8, "hPseudoResUvsV_cell[%d]_radius[%d]_pfx", icell, (r+1)*(45/rad));
      TH1 *h8 = (TH1*)(gDirectory->Get(hname8));
      
      if (h8){
	h8->GetXaxis()->SetTitle("Position reconstructed along Wire (cm)");
	h8->GetXaxis()->SetRangeUser(-0.6,0.6);
	h8->GetYaxis()->SetTitle("Position perp. to Wire (cm)");
	h8->GetYaxis()->SetRangeUser(-0.2,0.2);
	h8->GetYaxis()->SetTitleFont(42);
	h8->GetYaxis()->SetTitleSize(0.035);
	h8->GetYaxis()->SetLabelFont(42);
	h8->GetYaxis()->SetLabelSize(0.035);
	h8->SetLineColor(r+1);
	if (r == 0)
	  h8->Draw();
	else 
	  h8->Draw("same");
	
	TF1 *fp1 = new TF1("fp1","[0]+x*[1]",-0.45,0.45);
	fp1->SetLineColor(r+1);
	h8->Fit("fp1","qr");
	slope[r][icell-1] = fp1->GetParameter(1);
	slope_err[r][icell-1] = fp1->GetParError(1);

	//cout << slope[r][icell-1] << endl;
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

  TCanvas *cUncertain = new TCanvas("cUncertain", "Uncertainties", 1200, 600);
  cUncertain->Divide(2);
  cUncertain->cd(1);
  
  TH2 *hPseudoResSvsQ = (TH2*) (gDirectory->Get("hPseudoResSvsQ"));
  hPseudoResSvsQ->Draw("colz");

  cUncertain->cd(2);

  nslices = 150;
  TH1D *ResS[nslices];
  double q[nslices];
  double q_err[nslices];
  double sigma[nslices];
  double sigma_err[nslices];
  char hname[256];
  for (int i=0; i < nslices; i++) {
    q[i] = 2*i/10. + 1./10.  ;
    q_err[i] = 0;

    sprintf(hname, "ResS_[%d]", i);
    
    ResS[i] = hPseudoResSvsQ->ProjectionY(hname,  2*i, 2*(i+1) );
    ResS[i]->GetYaxis()->SetLabelFont(132);
    ResS[i]->GetYaxis()->SetTitleFont(132);
    ResS[i]->GetXaxis()->SetLabelFont(132);
    ResS[i]->GetXaxis()->SetTitleFont(132);
    ResS[i]->GetXaxis()->SetLabelSize(0.04);
    ResS[i]->GetXaxis()->SetTitleSize(0.04);
    
    //ResS[i]->Draw();
    
    if (ResS[i]->GetEntries() > 1000){
      TF1 *fgaus = new TF1("fgaus", "gaus", -0.05, 0.05);
      fgaus->SetLineColor(2);
      fgaus->SetNpx(360);
      ResS[i]->Fit("fgaus", "QR0");


      sigma[i] = fgaus->GetParameter(2);
      sigma_err[i] = fgaus->GetParError(2);

    }
    else {
      sigma[i] = 0;
      sigma_err[i] = 0;
    }
  }

  TH1D *hsigma = new TH1D("hsigma", "", nslices, 0, 30);
  hsigma->GetYaxis()->SetRangeUser(0,0.2);
  hsigma->GetYaxis()->SetTitleFont(132);
  hsigma->GetXaxis()->SetTitleFont(132);
  hsigma->GetYaxis()->SetLabelFont(132);
  hsigma->GetXaxis()->SetLabelFont(132);
  // hsigma_para->GetXaxis()->SetTitle("Beam Energy (GeV)");
  // hsigma_para->GetYaxis()->SetTitle("P#Sigma");
  hsigma->SetStats(0);
  hsigma->Draw();

  TGraph *gsigma = new TGraphErrors(nslices, q, sigma, q_err, sigma_err);
  gsigma->SetMarkerColor(1);
  gsigma->SetMarkerStyle(8);
  gsigma->Draw("Psame");

  TF1 *funct=new TF1("funct", "0.01/x + 0.01", 0 ,50);
  funct->SetLineColor(2);
  funct->Draw("same");

  TF1 *fit=new TF1("fit", "[0]/x + [1]", 0.4 , 3.5);
  fit->SetLineColor(4);
  gsigma->Fit("fit","R");
 
  TF1 *fit_full=new TF1("fit_full", "[0]/x + [1]", 0 , 50);
  fit_full->SetLineColor(4);
  fit_full->SetLineStyle(2);
  fit_full->SetParameter(0, fit->GetParameter(0));
  fit_full->SetParameter(1, fit->GetParameter(1));
  fit_full->Draw("same");


    
}
