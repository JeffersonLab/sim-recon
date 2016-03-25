// File: ST_Monitoring_1.C
// Created: 05/17/2015
// Creator: Mahmoud Kamel, mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes
{
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_multiplicity");
  if(dir) dir->cd();
  // Grab Multiplicity histograms
  TH2I *h2_adc_multip = (TH2I*)gDirectory->FindObjectAny("h2_sector_adc_multip"); 
  TH2I *h2_tdc_multip = (TH2I*)gDirectory->FindObjectAny("h2_sector_tdc_multip");
  TH2I *h2_hit_multip = (TH2I*)gDirectory->FindObjectAny("h2_sector_hit_multip");
  TH2I *h2_adc_unmatched = (TH2I*)gDirectory->FindObjectAny("h2_adc_unmatched"); 
  TH2I *h2_tdc_unmatched = (TH2I*)gDirectory->FindObjectAny("h2_tdc_unmatched"); 


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
  c1->Divide(3,1);
  // ADC 
  c1->cd(1);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogz();
  if (h2_adc_multip) h2_adc_multip->Draw("colz");
  //TDC 
  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogz();
  if (h2_tdc_multip) h2_tdc_multip->Draw("colz");
  // Hit multiplicity
  c1->cd(3);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogz();
  if (h2_hit_multip) h2_hit_multip->Draw("colz");
  

  TCanvas *c2 = new TCanvas("c2","Start Counter unmatched ADC & TDC hits", 200, 10, 600, 480);
  c2->cd(0);
  c2->Draw();
  c2->Update();
  TCanvas *c2 = gPad->GetCanvas();
  c2->Divide(2,1);
  // ADC 
  c2->cd(1);
  gPad->SetTicks();
  gPad->SetGrid();
  // gPad->SetLogz();
  if (h2_adc_unmatched) h2_adc_unmatched->Draw("colz");
  h2_adc_unmatched->GetZaxis()->SetRangeUser(100,10000.);
  // TDC 
  c2->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  // gPad->SetLogz();
  if (h2_tdc_unmatched) h2_tdc_unmatched->Draw("colz");
  h2_tdc_unmatched->GetZaxis()->SetRangeUser(100,10000.);


}
