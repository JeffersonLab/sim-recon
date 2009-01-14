{
  gStyle->SetOptFit();

  TFile *file = TFile::Open("fitter.root");

  cPtot = new TCanvas("cPtot","momentum residual",900,700);
  nFitter->Draw("1.0/ptinv - 1.0/ptinvTrue >> hPtot(200, -0.1, 0.1)", "nFDC == 0");
  hPtot->Fit("gaus");
  Double_t ptotSigma = 2.74990e-02;

  cTheta = new TCanvas("cTheta","theta residual",900,700);
  nFitter->Draw("theta - thetaTrue >> hTheta(200, -0.02, 0.02)", "nFDC == 0");
  hTheta->Fit("gaus");
  Double_t thetaSigma = 2.10375e-03;

  cPhi = new TCanvas("cPhi","phi residual",900,700);
  nFitter->Draw("phi - phiTrue >> hPhi(200, -0.02, 0.02)", "nFDC == 0");
  hPhi->Fit("gaus");
  Double_t phiSigma = 1.38125e-03;

  cProb = new TCanvas("prob","probability",900,700);
  nFitter->Draw("TMath::Prob(((1.0/ptinv - 1.0/ptinvTrue)/0.0126)^2 + ((theta - thetaTrue)/0.003273)^2 + ((phi - phiTrue)/0.005703)^2, 3) >> hProb(100, 0.0, 1.0)", "nFDC == 0");

  cPtot->Print("fitter_pt.png", "png");
  cTheta->Print("fitter_theta.png", "png");
  cPhi->Print("fitter_phi.png", "png");
  cProb->Print("fitter_prob.png", "png");

}
