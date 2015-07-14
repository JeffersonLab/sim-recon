// File: ST_Monitoring_1.C
// Created: 05/17/2015
// Creator: Mahmoud Kamel, mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes
{
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_lowlevel");
  if(dir) dir->cd();
  // Grab Hit object histograms
  TH2I *h2_t_sec = (TH2I*)gDirectory->FindObjectAny("h2_t_sec"); 
  TH2I *h2_tTDC_sec = (TH2I*)gDirectory->FindObjectAny("h2_tTDC_sec");
  TH2I *h2_tfADC_sec = (TH2I*)gDirectory->FindObjectAny("h2_tfADC_sec");
  TH2I *h2_dE_sec = (TH2I*)gDirectory->FindObjectAny("h2_dE_sec");
  
  //Grab TDC Digihitobject histograms 
  TH2I *h2_raw_tdcTime_sec = (TH2I*)gDirectory->FindObjectAny("h2_raw_tdcTime_sec");
  //if adc_sec = tdc_sec grab the TDC vs Sector histo
  TH2I *h2_tdcTime_sec = (TH2I*)gDirectory->FindObjectAny("h2_tdcTime_sec");
  //Create the canvas
  if(gPad == NULL)
    {
      TCanvas *c1 = new TCanvas("c1","Start Counter Raw Histograms 1", 200, 10, 600, 480);
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
  if (h2_tdcTime_sec) h2_tdcTime_sec->Draw("colz");
  h2_tdcTime_sec->SetLabelSize(0.025,"Y");
  h2_tdcTime_sec->GetXaxis()->CenterTitle();
  h2_tdcTime_sec->GetYaxis()->CenterTitle();
  h2_tdcTime_sec->GetYaxis()->SetTitleOffset(1.5);
  //t_tdc from hit object
  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_tTDC_sec) h2_tTDC_sec->Draw("colz");
  h2_tTDC_sec->GetXaxis()->CenterTitle();
  h2_tTDC_sec->GetYaxis()->CenterTitle();
  // t from hit object (tdc time walk corrected)
  c1->cd(3);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_t_sec) h2_t_sec->Draw("colz");
  h2_t_sec->GetXaxis()->CenterTitle();
  h2_t_sec->GetYaxis()->CenterTitle();
  c1->cd(4);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_tdcTime_sec = h2_tdcTime_sec->ProjectionY();
  h_tdcTime_sec->Draw();
  h_tdcTime_sec->GetXaxis()->CenterTitle();
  h_tdcTime_sec->GetYaxis()->CenterTitle();
  c1->cd(5);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_tTDC = h2_tTDC_sec->ProjectionY();
  h_tTDC->Draw();
  h_tTDC->GetXaxis()->CenterTitle();
  h_tTDC->GetYaxis()->CenterTitle();
  c1->cd(6);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_t = h2_t_sec->ProjectionY();
  h_t->Draw();
  h_t->GetXaxis()->CenterTitle();
  h_t->GetYaxis()->CenterTitle();
 
  //Create the canvas = 2
  TCanvas *c2 = new TCanvas("c2","Start Counter Raw Histograms 2", 200, 10, 600, 480);
  c2->cd(1);
  c2->Draw();
  c2->Update();
  if(!gPad) return;
  TCanvas *c2 = gPad->GetCanvas();
  c2->Divide(2,2);
  // t_adc from hit object
  c2->cd(1);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_tfADC_sec) h2_tfADC_sec->Draw("colz");
  h2_tfADC_sec->GetXaxis()->CenterTitle();
  h2_tfADC_sec->GetYaxis()->CenterTitle();
  // Energy loss from hit object
  c2->cd(2);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_dE_sec) h2_dE_sec->Draw("colz");
  h2_dE_sec->GetXaxis()->CenterTitle();
  h2_dE_sec->GetYaxis()->CenterTitle();
  c2->cd(3);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_tfADC = h2_tfADC_sec->ProjectionY();
  h_tfADC->Draw();
  //h1_adc_sec->SetMinimum(0);
  h_tfADC->GetXaxis()->CenterTitle();
  h_tfADC->GetYaxis()->CenterTitle();
  c2->cd(4);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  TH1D *h_dE = h2_dE_sec->ProjectionY();
  h_dE->Draw();
  h_dE->GetXaxis()->CenterTitle();
  h_dE->GetYaxis()->CenterTitle();
}
