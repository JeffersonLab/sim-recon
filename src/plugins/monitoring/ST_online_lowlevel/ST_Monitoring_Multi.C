// File: ST_Monitoring_1.C
// Created: 05/17/2015
// Creator: Mahmoud Kamel, mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes
{
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_lowlevel");
  if(dir) dir->cd();
  // Grab Multiplicity histograms
  TH2I *h2_st_adc_tdc_multi = (TH2I*)gDirectory->FindObjectAny("h2_st_adc_tdc_multi"); 
  TH2I *h2_st_adc_hit_multi = (TH2I*)gDirectory->FindObjectAny("h2_st_adc_hit_multi");
  TH1I *h1_adc_sec = (TH1I*)gDirectory->FindObjectAny("h1_adc_sec");
  TH1I *h1_tdc_sec = (TH1I*)gDirectory->FindObjectAny("h1_tdc_sec");
  TH1I *h1_hit_sec = (TH1I*)gDirectory->FindObjectAny("h1_hit_sec");
  //Create the canvas
  if(gPad == NULL)
    {
      TCanvas *c1 = new TCanvas("c1","Start Counter Multiplicity Histograms", 200, 10, 600, 480);
      c1->cd(0);
      c1->Draw();
      c1->Update();
    }
  
  if(!gPad) return;
  TCanvas *c1 = gPad->GetCanvas();
  c1->Divide(3,2);
  // ADC vs TDC
  c1->cd(1);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_st_adc_tdc_multi) h2_st_adc_tdc_multi->Draw("colz");
  h2_st_adc_tdc_multi->SetAxisRange(0., 40.,"X");
  h2_st_adc_tdc_multi->SetTitle("ST ADC Vs TDC multiplicities");
  h2_st_adc_tdc_multi->GetZaxis()->SetLabelFont(10);
  h2_st_adc_tdc_multi->GetZaxis()->SetLabelOffset(0.0);

  //ADC vs Hit
  c1->cd(4);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h2_st_adc_hit_multi) h2_st_adc_hit_multi->Draw("colz");
  h2_st_adc_hit_multi->SetAxisRange(0., 40.,"X");
  h2_st_adc_hit_multi->SetTitle("ST ADC Vs Hit multiplicities");
  h2_st_adc_hit_multi->GetZaxis()->SetLabelFont(10);
  h2_st_adc_hit_multi->GetZaxis()->SetLabelOffset(0.0);
  // TDC multiplicity
  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogy();
  TH1D *h_tdc_multi = h2_st_adc_tdc_multi->ProjectionX();
  h_tdc_multi->Draw();
  h_tdc_multi->SetAxisRange(0., 40.,"X");
  h_tdc_multi->GetYaxis()->SetTitle("Counts");
  h_tdc_multi->GetXaxis()->CenterTitle();
  h_tdc_multi->GetYaxis()->CenterTitle();
  // Hit multiplicity
  c1->cd(5);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogy();
  TH1D *h_hit_multi = h2_st_adc_hit_multi->ProjectionX();
  h_hit_multi->Draw();
  h_hit_multi->SetAxisRange(0., 40.,"X");
  h_hit_multi->GetYaxis()->SetTitle("Counts");
  h_hit_multi->GetXaxis()->CenterTitle();
  h_hit_multi->GetYaxis()->CenterTitle();
  // 3, and 6 are the same quantity from two different 2D histograms
  c1->cd(3);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogy();
  TH1D *h_adc_multi = h2_st_adc_tdc_multi->ProjectionY();
  h_adc_multi->Draw();
  h_adc_multi->SetAxisRange(0., 40.,"X");
  h_adc_multi->GetYaxis()->SetTitle("Counts");
  h_adc_multi->GetXaxis()->CenterTitle();
  h_adc_multi->GetYaxis()->CenterTitle();
 
  //***********  Create the OCCUPANCY canvas *****************************
  //**********************************************************************
  TCanvas *c2 = new TCanvas("c2","Start Counter Occupancy Histograms", 200, 10, 600, 480);
  c2->cd(0);
  c2->Draw();
  c2->Update();
  if(!gPad) return;
  TCanvas *c2 = gPad->GetCanvas();
  c2->Divide(3,1);
  
  // Occupancy Histos 
  c2->cd(1);
  gStyle->SetOptStat(10);
  gStyle->SetErrorX(0); 
  gPad->SetTicks();
  gPad->SetGrid();
  if (h1_tdc_sec) h1_tdc_sec->Draw("err");
  h1_tdc_sec->SetMinimum(0);
 
  h1_tdc_sec->GetXaxis()->CenterTitle();
  h1_tdc_sec->GetYaxis()->CenterTitle();
  h1_tdc_sec->GetYaxis()->SetTitleOffset(1.5);
  h1_tdc_sec->SetMarkerStyle(21);
  h1_tdc_sec->SetMarkerSize(1.2);
  h1_tdc_sec->SetMarkerColor(4);
  c2->cd(3);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h1_hit_sec) h1_hit_sec->Draw("err");
  h1_hit_sec->SetMinimum(0);
  h1_hit_sec->GetXaxis()->CenterTitle();
  h1_hit_sec->GetYaxis()->CenterTitle();
  h1_hit_sec->GetYaxis()->SetTitleOffset(1.5);
  h1_hit_sec->SetMarkerStyle(21);
  h1_hit_sec->SetMarkerSize(1.2);
  h1_hit_sec->SetMarkerColor(4);
  c2->cd(2);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h1_adc_sec) h1_adc_sec->Draw("err");
  h1_adc_sec->SetMinimum(0);
  h1_adc_sec->GetXaxis()->CenterTitle();
  h1_adc_sec->GetYaxis()->CenterTitle();
  h1_adc_sec->GetYaxis()->SetTitleOffset(1.5);
  h1_adc_sec->SetMarkerStyle(21);
  h1_adc_sec->SetMarkerSize(1.2);
  h1_adc_sec->SetMarkerColor(4);
}
