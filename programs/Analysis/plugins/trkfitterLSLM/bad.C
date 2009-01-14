{
  TFile *file = TFile::Open("fitter.root");
  nFitter->Draw("phi - phiTrue");
  c1->Print("badphi.png", "png");
  nFitter->Draw("theta");
  c1->Print("badtheta.png", "png");
  nFitter->Draw("1.0/ptinv/TMath::Sin(theta)");
  c1->Print("badptot.png", "png");
}
