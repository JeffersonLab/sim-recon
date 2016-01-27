// File: ST_Monitoring_2D.C
// Created: 05/21/2014
// Creator: Mahmoud Kamel, mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes

{
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_tracking");
  //gDirectory->ls();
  if(dir) dir->cd();

  // Grab 2D histograms  
  TH2I *h_r_vs_z       = (TH2I*)gDirectory->FindObjectAny("r_vs_z");
  TH2I *h_y_vs_x       = (TH2I*)gDirectory->FindObjectAny("y_vs_x");
  TH2I *h_phi_vs_sector= (TH2I*)gDirectory->FindObjectAny("phi_vs_sector");
  TH2I *h2_tdct_adct_vs_de = (TH2I*)gDirectory->FindObjectAny("h2_tdct_adct_vs_de");
  TH2I *h2_tdctCorr_adct_vs_de = (TH2I*)gDirectory->FindObjectAny("h2_tdctCorr_adct_vs_de");
  TH2I *h2_hit_det_eff = (TH2I*)gDirectory->FindObjectAny("h2_hit_det_eff");
  TH1I *h_hit_det_eff= (TH1I*)gDirectory->FindObjectAny("h_hit_det_eff");
  TH2I *h2_dphi_sector = (TH2I*)gDirectory->FindObjectAny("h2_dphi_sector"); 
  TH2I *h2_adc_det_eff= (TH2I*)gDirectory->FindObjectAny("h2_adc_det_eff");
  //Create the canvas
  if(gPad == NULL)
    {
      TCanvas *c1 = new TCanvas("c1","Start Counter 2D Histograms", 200, 10, 600, 480);
      c1->cd(0);
      c1->Draw();
      c1->Update();
    }
  
  if(!gPad) return;
  TCanvas *c1 = gPad->GetCanvas();
  c1->Divide(2,2);

  // Pulse integral vs. pulse time histogram
  c1->cd(1);
  gPad->SetTicks();
  gPad->SetGrid();
  // pi_pt_dhit_2d->SetStats(0);
  // pi_pt_dhit_2d->GetXaxis()->CenterTitle();
  // pi_pt_dhit_2d->GetYaxis()->CenterTitle();
  if(h_r_vs_z) h_r_vs_z->Draw("colz");
  gPad->SetLogz();
  gStyle->SetOptStat(10);
  // Pulse integral vs. f1TDC DigiHit time
  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogz();
  // pt_tdc_dhit_2d->SetStats(0);
  // pt_tdc_dhit_2d->GetYaxis()->CenterTitle();
  // pt_tdc_dhit_2d->GetXaxis()->CenterTitle();
  if(h_y_vs_x) h_y_vs_x->Draw("colz");
  h_y_vs_x->GetZaxis()->SetRangeUser(0.0, 100.);
  // f1TDC DigiHit time vs. channel number
  c1->cd(3);
  gPad->SetTicks();
  gPad->SetGrid();
  // tdc_dhit_2d->SetStats(0);
  // tdc_dhit_2d->SetAxisRange(-1., 31.,"X");
  // tdc_dhit_2d->GetXaxis()->CenterTitle();
  // tdc_dhit_2d->GetYaxis()->CenterTitle();
  if(h_phi_vs_sector) h_phi_vs_sector->Draw("colz");
  // fADC250 pulse integral vs. channel histogram
  c1->cd(4);
  gPad->SetTicks();
  gPad->SetGrid();
  // pi_dhit_2d->SetStats(0);
  // pi_dhit_2d->SetAxisRange(-1., 31.,"X");
  // pi_dhit_2d->GetXaxis()->CenterTitle();
  // pi_dhit_2d->GetYaxis()->CenterTitle();
  if(h2_dphi_sector) h2_dphi_sector->Draw("colz");
 
}
