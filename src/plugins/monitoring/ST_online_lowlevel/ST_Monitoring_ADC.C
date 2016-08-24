// File: ST_Monitoring_1.C
// Created: 05/17/2015
// Creator: Mahmoud Kamel, mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes
{
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_lowlevel");
  if(dir) dir->cd();
  // Grab Raw ADC object histograms
  TH2I *h2_raw_pi_sector = (TH2I*)gDirectory->FindObjectAny("h2_raw_pi_sector"); 
  TH2I *h2_raw_ped_sector = (TH2I*)gDirectory->FindObjectAny("h2_raw_ped_sector");
  TH2I *h2_raw_pt_sector = (TH2I*)gDirectory->FindObjectAny("h2_raw_pt_sector");
    
  //Grab Offset applied ADC histograms 
  
  TH2I *h2_adc_pcpi_sector = (TH2I*)gDirectory->FindObjectAny("h2_adc_pcpi_sector");
  TH2I *h2_adc_pt_sector = (TH2I*)gDirectory->FindObjectAny("h2_adc_pt_sector");
  TH2I *h2_adc_ped_sector = (TH2I*)gDirectory->FindObjectAny("h2_adc_ped_sector");
  TH2I *h2_adc_pp_sector = (TH2I*)gDirectory->FindObjectAny("h2_adc_pp_sector");
  TH2I *h2_st_time_vs_pcpi= (TH2I*)gDirectory->FindObjectAny("h2_st_time_vs_pcpi");
  TH2I *h2_st_time_vs_pp= (TH2I*)gDirectory->FindObjectAny("h2_st_time_vs_pp");
  //Create the canvas
  if(gPad == NULL)
    {
      TCanvas *c1 = new TCanvas("c1","Start Counter ADC Histograms 1 (None offset)", 200, 10, 600, 480);
      c1->cd(0);
      c1->Draw();
      c1->Update();
    }
  if(!gPad) return;
  TCanvas *c1 = gPad->GetCanvas();
  c1->Divide(3,2);
  //TDC time when there is a hit in adc
  c1->cd(1);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_raw_pi_sector) h2_raw_pi_sector->Draw("colz");
  h2_raw_pi_sector->SetLabelSize(0.025,"Y");
  h2_raw_pi_sector->GetXaxis()->CenterTitle();
  h2_raw_pi_sector->GetYaxis()->CenterTitle();
  //t_tdc from hit object
  h2_raw_pi_sector->GetYaxis()->SetTitleOffset(1.5);
  gPad->SetLogz();
  c1->cd(2);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_adc_ped_sector) h2_adc_ped_sector->Draw("colz");
  h2_adc_ped_sector->GetXaxis()->CenterTitle();
  h2_adc_ped_sector->GetYaxis()->CenterTitle();
  h2_adc_ped_sector->GetYaxis()->SetTitleOffset(1.6);
  gPad->SetLogz();
  // t from hit object (tdc time walk corrected)
  c1->cd(3);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_raw_pt_sector) h2_raw_pt_sector->Draw("colz");
  h2_raw_pt_sector->GetXaxis()->CenterTitle();
  h2_raw_pt_sector->GetYaxis()->CenterTitle();
  h2_raw_pt_sector->GetYaxis()->SetTitleOffset(1.55);
  gPad->SetLogz();
  c1->cd(4);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_pi = h2_raw_pi_sector->ProjectionY();
  h_pi->Draw();
  h_pi->GetXaxis()->CenterTitle();
  h_pi->GetYaxis()->CenterTitle();
  c1->cd(5);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_adc_ped = h2_adc_ped_sector->ProjectionY();
  h_adc_ped->Draw();
  h_adc_ped->GetXaxis()->CenterTitle();
  h_adc_ped->GetYaxis()->CenterTitle();
 
  c1->cd(6);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_pt = h2_raw_pt_sector->ProjectionY();
  h_pt->Draw();
  h_pt->GetXaxis()->CenterTitle();
  h_pt->GetYaxis()->CenterTitle();
 
  //Create the canvas = 2
  TCanvas *c2 = new TCanvas("c2","Start Counter ADC Histograms 2(offsets applied)", 200, 10, 600, 480);
  c2->cd(1);
  c2->Draw();
  c2->Update();
  if(!gPad) return;
  c2 = gPad->GetCanvas();
  c2->Divide(3,2);
  // t_adc from hit object
  c2->cd(1);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_adc_pcpi_sector) h2_adc_pcpi_sector->Draw("colz");
  h2_adc_pcpi_sector->GetXaxis()->CenterTitle();
  h2_adc_pcpi_sector->GetYaxis()->CenterTitle(); 
  h2_adc_pcpi_sector->GetYaxis()->SetTitleOffset(1.5);
  h2_adc_pcpi_sector->SetLabelSize(0.03,"Y");
  gPad->SetLogz();
  // Energy loss from hit object
  c2->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_raw_ped_sector) h2_raw_ped_sector->Draw("colz");
  h2_raw_ped_sector->GetYaxis()->SetTitleOffset(1.5);
  h2_raw_ped_sector->GetXaxis()->CenterTitle();
  h2_raw_ped_sector->GetYaxis()->CenterTitle();
  gPad->SetLogz();
  c2->cd(3);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_adc_pt_sector) h2_adc_pt_sector->Draw("colz");
  h2_adc_pt_sector->GetXaxis()->CenterTitle();
  h2_adc_pt_sector->GetYaxis()->CenterTitle();
  gPad->SetLogz();
  c2->cd(4);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_adc_pcpi = h2_adc_pcpi_sector->ProjectionY();
  h_adc_pcpi->Draw();
  h_adc_pcpi->GetXaxis()->CenterTitle();
  h_adc_pcpi->GetYaxis()->CenterTitle();
  c2->cd(5);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_ped = h2_raw_ped_sector->ProjectionY();
  h_ped->Draw();
  h_ped->GetXaxis()->CenterTitle();
  h_ped->GetYaxis()->CenterTitle();
  c2->cd(6);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_adc_pt = h2_adc_pt_sector->ProjectionY();
  h_adc_pt->Draw();
  h_adc_pt->GetXaxis()->CenterTitle();
  h_adc_pt->GetYaxis()->CenterTitle();
  
  //Create the canvas = 3
  TCanvas *c3 = new TCanvas("c3","Start Counter ADC Histograms 3", 200, 10, 600, 480);
  c3->cd(1);
  c3->Draw();
  c3->Update();
  if(!gPad) return;
  c3 = gPad->GetCanvas();
  c3->Divide(3,2);
  // t_adc from hit object
  c3->cd(1);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_adc_pp_sector) h2_adc_pp_sector->Draw("colz");
  h2_adc_pp_sector->GetXaxis()->CenterTitle();
  h2_adc_pp_sector->GetYaxis()->CenterTitle();
  h2_adc_pp_sector->GetYaxis()->SetTitleOffset(1.5);
   gPad->SetLogz();
  c3->cd(2);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_st_time_vs_pp) h2_st_time_vs_pp->Draw("colz");
  h2_st_time_vs_pp->GetXaxis()->CenterTitle();
  h2_st_time_vs_pp->GetYaxis()->CenterTitle();
  h2_st_time_vs_pp->GetYaxis()->SetTitleOffset(1.4);
  gPad->SetLogz();
  c3->cd(3);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_st_time_vs_pcpi) h2_st_time_vs_pcpi->Draw("colz");
  h2_st_time_vs_pcpi->GetXaxis()->CenterTitle();
  h2_st_time_vs_pcpi->GetYaxis()->CenterTitle();
  gPad->SetLogz();
  c3->cd(4);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_adc_pp = h2_adc_pp_sector->ProjectionY();
  h_adc_pp->Draw();
  h_adc_pp->GetXaxis()->CenterTitle();
  h_adc_pp->GetYaxis()->CenterTitle();
  c3->cd(5);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_adc_delta_t = h2_st_time_vs_pp->ProjectionY();
  h_adc_delta_t->Draw();
  h_adc_delta_t->GetXaxis()->CenterTitle();
  h_adc_delta_t->GetYaxis()->CenterTitle();
  c3->cd(6);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_pcpi = h2_st_time_vs_pcpi->ProjectionX();
  h_pcpi->Draw();
  h_pcpi->GetXaxis()->CenterTitle();
  h_pcpi->GetYaxis()->CenterTitle();
}
