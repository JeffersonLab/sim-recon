{
  TFile *file = TFile::Open("fitter.root");

  nFitter->Draw("1.0/ptinv >> hPt(200, 0.0, 4.0)");
  TH1F *hPt = (TH1F*)gDirectory->Get("hPt");
  c1->SetLeftMargin(0.15);
  hPt->SetTitle("transverse momentum (p_{t})");
  hPt->GetXaxis()->SetTitle("p_{t} (GeV/c)");
  hPt->GetYaxis()->SetTitle("Events / 2 MeV/c");
  hPt->GetYaxis()->SetTitleOffset(1.5);
  hPt->Draw();
  c1->Print("fitter_p.png", "png");

  nFitter->Draw("theta >> hTheta");
  TH1F *hTheta = (TH1F*)gDirectory->Get("hTheta");
  hTheta->SetTitle("polar angle (#theta)");
  hTheta->GetXaxis()->SetTitle("#theta (radians)");
  hTheta->Draw();
  c1->Print("fitter_theta.png", "png");

  nFitter->Draw("phi >> hPhi");
  TH1F *hPhi = (TH1F*)gDirectory->Get("hPhi");
  hPhi->SetTitle("azimuthal angle (#phi)");
  hPhi->GetXaxis()->SetTitle("#phi (radians)");
  hPhi->Draw();
  c1->Print("fitter_phi.png", "png");

  nFitter->Draw("xp0 >> hXp0(100, -0.5, 0.5)");
  TH1F *hXp0 = (TH1F*)gDirectory->Get("hXp0");
  hXp0->SetTitle("x'_{0}");
  hXp0->GetXaxis()->SetTitle("x'_{0} (cm)");
  hXp0->GetYaxis()->SetTitle("Events / 100 #mum");
  hXp0->GetYaxis()->SetTitleOffset(1.5);
  hXp0->Draw();
  c1->Print("fitter_xp0.png", "png");

  nFitter->Draw("z0 >> hZ0(100, 62.5, 67.5)");
  hZ0->SetTitle("z_{0}");
  hZ0->GetXaxis()->SetTitle("z_{0} (cm)");
  hZ0->GetYaxis()->SetTitle("Events / 500 #mum");
  hZ0->GetYaxis()->SetTitleOffset(1.5);
  hZ0->Draw();
  c1->Print("fitter_z0.png", "png");

  nFitter->Draw("nCDC:nFDC >> h2NHits(30, -0.5, 29.5, 30, -0.5, 29.5)");
  TH1F *h2NHits = (TH1F*)gDirectory->Get("h2NHits");
  h2NHits->SetTitle("n_{CDC} vs. n_{FDC}");
  h2NHits->GetXaxis()->SetTitle("n_{FDC}");
  h2NHits->GetYaxis()->SetTitle("n_{CDC}");
  h2NHits->Draw("BOX");
  c1->Print("fitter_f_vs_c.png", "png");

  nFitter->Draw("nCDC:nFDC >> h2NHits2(29, 0.5, 29.5, 29, 0.5, 29.5)");
  TH1F *h2NHits2 = (TH1F*)gDirectory->Get("h2NHits2");
  h2NHits2->SetTitle("n_{CDC} vs. n_{FDC}, mixed events");
  h2NHits2->GetXaxis()->SetTitle("n_{FDC}");
  h2NHits2->GetYaxis()->SetTitle("n_{CDC}");
  h2NHits2->Draw("BOX");
  c1->Print("fitter_f_vs_c2.png", "png");

  nFitter->Draw("1.0/ptinv:theta >> h2PvsTheta(60, 0.0, 3.0, 50, 1.5, 2.5)");
  gStyle->SetPalette(1);
  h2PvsTheta->SetTitle("transverse momentum vs. polar angle");
  h2PvsTheta->GetXaxis()->SetTitle("#theta (radians)");
  h2PvsTheta->GetXaxis()->SetTitleOffset(1.5);
  h2PvsTheta->GetYaxis()->SetTitle("p_{t} (GeV/c)");
  h2PvsTheta->GetYaxis()->SetTitleOffset(2.0);
  h2PvsTheta->Draw("LEGO1");
  c1->Print("fitter_p_vs_theta.png", "png");

}
