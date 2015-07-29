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
  gStyle->SetOptStat(10);
  // Pulse integral vs. f1TDC DigiHit time
  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
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
  // h2_tdct_adct_vs_de->SetAxisRange(0.0, 0.006.,"X");
  // h2_tdct_adct_vs_de->SetAxisRange(-2., 2.,"Y");
  // fADC250 pulse pedestal vs. channel histogram
  // c1->cd(5);
 //  gPad->SetTicks();
 //  gPad->SetGrid();
 //  // ped_dhit_2d->SetStats(0);
 //  // ped_dhit_2d->SetAxisRange(-1., 31.,"X");
 //  // ped_dhit_2d->GetXaxis()->CenterTitle();
 //  // ped_dhit_2d->GetYaxis()->CenterTitle();
 //  if(h2_tdctCorr_adct_vs_de) h2_tdctCorr_adct_vs_de->Draw("colz");
 //  h2_tdctCorr_adct_vs_de->SetAxisRange(0.0, 0.006.,"X");
 //  h2_tdctCorr_adct_vs_de->SetAxisRange(-2., 2.,"Y");
 //  // fADC250 pulse time vs. channel histogram
 //  c1->cd(6);
 //  gPad->SetTicks();
 //  gPad->SetGrid();
 //  // pt_dhit_2d->SetStats(0);
 //  // pt_dhit_2d->SetAxisRange(-1., 31.,"X");
 //  // pt_dhit_2d->GetXaxis()->CenterTitle();
 //  // pt_dhit_2d->GetYaxis()->CenterTitle();
 //  if(h2_hit_det_eff) h2_hit_det_eff->Draw("colz");
 //  // ST fADC250 multiplicity vs. channel histogram
 //  c1->cd(7);
 //  gPad->SetTicks();
 //  gPad->SetGrid();
 //  // adc_multi_2d->SetStats(0);
 //  // adc_multi_2d->SetAxisRange(-1., 31.,"X");
 //  // adc_multi_2d->GetXaxis()->CenterTitle();
 //  // adc_multi_2d->GetYaxis()->CenterTitle();
 //  if(h_hit_det_eff) h_hit_det_eff->Draw("");
 //  // ST f1TDC multiplicty vs. channel histogram
 //  c1->cd(8);
 //  gPad->SetTicks();
 //  gPad->SetGrid();
 //  // tdc_multi_2d->SetStats(0);
 //  // tdc_multi_2d->SetAxisRange(-1., 31.,"X");
 //  // tdc_multi_2d->GetXaxis()->CenterTitle();
 //  // tdc_multi_2d->GetYaxis()->CenterTitle();
 //  if(h2_dt_z) h2_dt_z->Draw("colz");
 // // ST ADC vs TDC Multiplicity histogram
 //  c1->cd(9);
 //  gPad->SetTicks();
 //  gPad->SetGrid();
 //  // adc_tdc_multi_2d->SetStats(0);
 //  // adc_tdc_multi_2d->GetXaxis()->CenterTitle();
 //  // adc_tdc_multi_2d->GetYaxis()->CenterTitle();
 //  if(h2_adc_det_eff) h2_adc_det_eff->Draw("colz");
}
