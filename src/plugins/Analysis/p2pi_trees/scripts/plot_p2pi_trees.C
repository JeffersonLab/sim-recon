void plot_p2pi_trees(void)
{
// File: p2pi_trees.C
//

  gStyle->SetPalette(1,0);
  gStyle->SetOptStat(111111);
  gStyle->SetOptFit(111111);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadBottomMargin(0.15);
    
    TString filename = "DSelector_p2pi_trees_geant3";
    
    TFile *f = new TFile(filename+".root","read");
    
    cout << "Opening root file: " << (filename+".root").Data() << endl;
    
    TH1F *NumEventsSurvivedAction = (TH1F*)f->Get("NumEventsSurvivedAction");
    TH2F *NumCombosSurvivedAction = (TH2F*)f->Get("NumCombosSurvivedAction");
    TH1F *BeamEnergy = (TH1F*)f->Get("BeamEnergy");
    TH1F *pMomentumMeasured = (TH1F*)f->Get("pMomentumMeasured");
    TH1F *piplusMomentumMeasured = (TH1F*)f->Get("piplusMomentumMeasured");
    TH1F *piminusMomentumMeasured = (TH1F*)f->Get("piminusMomentumMeasured");
    
    TH1F *KinFitChiSq = (TH1F*)f->Get("KinFitChiSq");
    TH1F *KinFitCL = (TH1F*)f->Get("KinFitCL");
    TH1F *MissingMassSquared = (TH1F*)f->Get("MissingMassSquared");
    TH1F *M2pi = (TH1F*)f->Get("M2pi");
    TH1F *t = (TH1F*)f->Get("t");
    TH2F *CosTheta_Psi = (TH2F*)f->Get("CosTheta_Psi");
    
    TH1F *pDeltap_Measured = (TH1F*)f->Get("pDeltap_Measured");
    TH1F *pipDeltap_Measured = (TH1F*)f->Get("pipDeltap_Measured");
    TH1F *pimDeltap_Measured = (TH1F*)f->Get("pimDeltap_Measured");
    TH1F *pDeltap = (TH1F*)f->Get("pDeltap");
    TH1F *pipDeltap = (TH1F*)f->Get("pipDeltap");
    TH1F *pimDeltap = (TH1F*)f->Get("pimDeltap");
    
   TCanvas *c0 = new TCanvas("c0", "c0",200,10,1000,700);

   c0->Divide(3,2);
    c0->cd(1);
    // gPad->SetLogy();
    Double_t xmin = 0;
    Double_t xmax = 2;
    Double_t ymin = 100;
    Double_t ymax = 10000;
        
    NumEventsSurvivedAction->SetTitle(filename);
    // NumEventsSurvivedAction->GetXaxis()->SetRangeUser(xmin,xmax);
    // NumEventsSurvivedAction->GetYaxis()->SetRangeUser(ymin,ymax);
    NumEventsSurvivedAction->GetXaxis()->SetTitleSize(0.05);
    NumEventsSurvivedAction->GetYaxis()->SetTitleSize(0.05);
    // NumEventsSurvivedAction->GetXaxis()->SetTitle("Events");
    NumEventsSurvivedAction->SetMarkerColor(4);
    NumEventsSurvivedAction->Draw();
    
    c0->cd(2);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 5000;
    
    NumCombosSurvivedAction->SetTitle(filename);
    // NumCombosSurvivedAction->GetXaxis()->SetRangeUser(xmin,xmax);
    // NumCombosSurvivedAction->GetYaxis()->SetRangeUser(ymin,ymax);
    NumCombosSurvivedAction->GetXaxis()->SetTitleSize(0.05);
    NumCombosSurvivedAction->GetYaxis()->SetTitleSize(0.05);
    // NumCombosSurvivedAction->GetXaxis()->SetTitle("Events");
    // NumCombosSurvivedAction->SetMarkerColor(4);
    NumCombosSurvivedAction->Draw("colz");
    
    c0->cd(3);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    BeamEnergy->SetTitle(filename);
    // BeamEnergy->GetXaxis()->SetRangeUser(xmin,xmax);
    // BeamEnergy->GetYaxis()->SetRangeUser(ymin,ymax);
    BeamEnergy->GetXaxis()->SetTitleSize(0.05);
    BeamEnergy->GetYaxis()->SetTitleSize(0.05);
    BeamEnergy->GetXaxis()->SetTitle("Energy (GeV)");
    BeamEnergy->SetMarkerColor(4);
    BeamEnergy->Draw();
    
    c0->cd(4);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    pMomentumMeasured->SetTitle(filename);
    // pMomentumMeasured->GetXaxis()->SetRangeUser(xmin,xmax);
    // pMomentumMeasured->GetYaxis()->SetRangeUser(ymin,ymax);
    pMomentumMeasured->GetXaxis()->SetTitleSize(0.05);
    pMomentumMeasured->GetYaxis()->SetTitleSize(0.05);
    pMomentumMeasured->GetXaxis()->SetTitle("Proton Momentum (GeV)");
    pMomentumMeasured->SetMarkerColor(4);
    pMomentumMeasured->Draw();
    
    c0->cd(5);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    piplusMomentumMeasured->SetTitle(filename);
    // piplusMomentumMeasured->GetXaxis()->SetRangeUser(xmin,xmax);
    // piplusMomentumMeasured->GetYaxis()->SetRangeUser(ymin,ymax);
    piplusMomentumMeasured->GetXaxis()->SetTitleSize(0.05);
    piplusMomentumMeasured->GetYaxis()->SetTitleSize(0.05);
    piplusMomentumMeasured->GetXaxis()->SetTitle("PiPlus Momentum (GeV)");
    piplusMomentumMeasured->SetMarkerColor(4);
    piplusMomentumMeasured->Draw();
    
    c0->cd(6);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    piminusMomentumMeasured->SetTitle(filename);
    // piminusMomentumMeasured->GetXaxis()->SetRangeUser(xmin,xmax);
    // piminusMomentumMeasured->GetYaxis()->SetRangeUser(ymin,ymax);
    piminusMomentumMeasured->GetXaxis()->SetTitleSize(0.05);
    piminusMomentumMeasured->GetYaxis()->SetTitleSize(0.05);
    piminusMomentumMeasured->GetXaxis()->SetTitle("PiMinus Momentum (GeV)");
    piminusMomentumMeasured->SetMarkerColor(4);
    piminusMomentumMeasured->Draw();
    
    
    TCanvas *c1 = new TCanvas("c1", "c1",200,10,1000,700);
    
    c1->Divide(3,2);
    
    
    c1->cd(1);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    KinFitChiSq->SetTitle(filename);
    // KinFitChiSq->GetXaxis()->SetRangeUser(xmin,xmax);
    // KinFitChiSq->GetYaxis()->SetRangeUser(ymin,ymax);
    KinFitChiSq->GetXaxis()->SetTitleSize(0.05);
    KinFitChiSq->GetYaxis()->SetTitleSize(0.05);
    KinFitChiSq->GetXaxis()->SetTitle("#chi^{2}");
    KinFitChiSq->SetMarkerColor(4);
    KinFitChiSq->Draw();
    
    c1->cd(2);
    gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    KinFitCL->SetTitle(filename);
    // KinFitCL->GetXaxis()->SetRangeUser(xmin,xmax);
    // KinFitCL->GetYaxis()->SetRangeUser(ymin,ymax);
    KinFitCL->GetXaxis()->SetTitleSize(0.05);
    KinFitCL->GetYaxis()->SetTitleSize(0.05);
    KinFitCL->GetXaxis()->SetTitle("Probability");
    KinFitCL->SetMarkerColor(4);
    KinFitCL->Draw();
    
    c1->cd(3);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    MissingMassSquared->SetTitle(filename);
    // MissingMassSquared->GetXaxis()->SetRangeUser(xmin,xmax);
    // MissingMassSquared->GetYaxis()->SetRangeUser(ymin,ymax);
    MissingMassSquared->GetXaxis()->SetTitleSize(0.05);
    MissingMassSquared->GetYaxis()->SetTitleSize(0.05);
    MissingMassSquared->GetXaxis()->SetTitle("Missing Mass (GeV2)");
    MissingMassSquared->SetMarkerColor(4);
    MissingMassSquared->Draw();
    
    c1->cd(4);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    M2pi->SetTitle(filename);
    // M2pi->GetXaxis()->SetRangeUser(xmin,xmax);
    // M2pi->GetYaxis()->SetRangeUser(ymin,ymax);
    M2pi->GetXaxis()->SetTitleSize(0.05);
    M2pi->GetYaxis()->SetTitleSize(0.05);
    M2pi->GetXaxis()->SetTitle("Mass 2pi (GeV)");
    M2pi->SetMarkerColor(4);
    M2pi->Draw();
    
    c1->cd(5);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    t->SetTitle(filename);
    // t->GetXaxis()->SetRangeUser(xmin,xmax);
    // t->GetYaxis()->SetRangeUser(ymin,ymax);
    t->GetXaxis()->SetTitleSize(0.05);
    t->GetYaxis()->SetTitleSize(0.05);
    t->GetXaxis()->SetTitle("-t (GeV2)");
    t->SetMarkerColor(4);
    t->Draw();
    
    c1->cd(6);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    CosTheta_Psi->SetTitle(filename);
    // CosTheta_Psi->GetXaxis()->SetRangeUser(xmin,xmax);
    // CosTheta_Psi->GetYaxis()->SetRangeUser(ymin,ymax);
    CosTheta_Psi->GetXaxis()->SetTitleSize(0.05);
    CosTheta_Psi->GetYaxis()->SetTitleSize(0.05);
    CosTheta_Psi->GetYaxis()->SetTitle("CosTheta");
    CosTheta_Psi->GetXaxis()->SetTitle("#Psi (Degrees)");
    CosTheta_Psi->SetMarkerColor(4);
    CosTheta_Psi->Draw("colz");
    
    
    TCanvas *c2 = new TCanvas("c2", "c2",200,10,1000,700);
    
    c2->Divide(3,2);
    
    c2->cd(1);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    pDeltap_Measured->SetTitle(filename);
    // pDeltap_Measured->GetXaxis()->SetRangeUser(xmin,xmax);
    // pDeltap_Measured->GetYaxis()->SetRangeUser(ymin,ymax);
    pDeltap_Measured->GetXaxis()->SetTitleSize(0.05);
    pDeltap_Measured->GetYaxis()->SetTitleSize(0.05);
    // pDeltap_Measured->GetXaxis()->SetTitle("");
    pDeltap_Measured->SetMarkerColor(4);
    pDeltap_Measured->Draw();
    
    c2->cd(2);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    pipDeltap_Measured->SetTitle(filename);
    // pipDeltap_Measured->GetXaxis()->SetRangeUser(xmin,xmax);
    // pipDeltap_Measured->GetYaxis()->SetRangeUser(ymin,ymax);
    pipDeltap_Measured->GetXaxis()->SetTitleSize(0.05);
    pipDeltap_Measured->GetYaxis()->SetTitleSize(0.05);
    // pipDeltap_Measured->GetXaxis()->SetTitle("#chi^{2}");
    pipDeltap_Measured->SetMarkerColor(4);
    pipDeltap_Measured->Draw();
    
    c2->cd(3);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    pimDeltap_Measured->SetTitle(filename);
    // pimDeltap_Measured->GetXaxis()->SetRangeUser(xmin,xmax);
    // pimDeltap_Measured->GetYaxis()->SetRangeUser(ymin,ymax);
    pimDeltap_Measured->GetXaxis()->SetTitleSize(0.05);
    pimDeltap_Measured->GetYaxis()->SetTitleSize(0.05);
    // pimDeltap_Measured->GetXaxis()->SetTitle("#chi^{2}");
    pimDeltap_Measured->SetMarkerColor(4);
    pimDeltap_Measured->Draw();
    
    c2->cd(4);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    pDeltap->SetTitle(filename);
    // pDeltap->GetXaxis()->SetRangeUser(xmin,xmax);
    // pDeltap->GetYaxis()->SetRangeUser(ymin,ymax);
    pDeltap->GetXaxis()->SetTitleSize(0.05);
    pDeltap->GetYaxis()->SetTitleSize(0.05);
    // pDeltap->GetXaxis()->SetTitle("");
    pDeltap->SetMarkerColor(4);
    pDeltap->Draw();
    
    c2->cd(5);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    pipDeltap->SetTitle(filename);
    // pipDeltap->GetXaxis()->SetRangeUser(xmin,xmax);
    // pipDeltap->GetYaxis()->SetRangeUser(ymin,ymax);
    pipDeltap->GetXaxis()->SetTitleSize(0.05);
    pipDeltap->GetYaxis()->SetTitleSize(0.05);
    // pipDeltap->GetXaxis()->SetTitle("#chi^{2}");
    pipDeltap->SetMarkerColor(4);
    pipDeltap->Draw();
    
    c2->cd(6);
    // gPad->SetLogy();
    ymin = 0;
    ymax = 4000;
    
    pimDeltap->SetTitle(filename);
    // pimDeltap->GetXaxis()->SetRangeUser(xmin,xmax);
    // pimDeltap->GetYaxis()->SetRangeUser(ymin,ymax);
    pimDeltap->GetXaxis()->SetTitleSize(0.05);
    pimDeltap->GetYaxis()->SetTitleSize(0.05);
    // pimDeltap->GetXaxis()->SetTitle("#chi^{2}");
    pimDeltap->SetMarkerColor(4);
    pimDeltap->Draw();
    
    
    c0->SaveAs(filename+".pdf(");
    c1->SaveAs(filename+".pdf");
    c2->SaveAs(filename+".pdf)");
}
