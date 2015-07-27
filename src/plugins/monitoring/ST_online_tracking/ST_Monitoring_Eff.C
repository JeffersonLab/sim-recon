// File: ST_Monitoring_Efficiency.C
// Last Modified: 05/27/2015
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes
{
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_tracking");
  //gDirectory->ls();
  if(dir) dir->cd();

  // Grab 1D histograms 
  TH1D *pEff     = (TH1D*)dir->FindObjectAny("pEff");
  TH1D *pEff_adc     = (TH1D*)dir->FindObjectAny("pEff_adc");
  

  //Create the canvas
  if(gPad == NULL)
    {
      TCanvas *c1 = new TCanvas("c1","Start Counter 1D Histograms", 200, 10, 600, 480);
      c1->cd(0);
      c1->Draw();
      c1->Update();
    }
  
  if(!gPad) return;
  TCanvas *c1 = gPad->GetCanvas();
  c1->Divide(2,1);
  // ST ADC Efficiency histogram
  c1->cd(1);

  gStyle->SetStatW(0.13);
  gStyle->SetStatH(0.05);
  gStyle->SetStatY(0.97);
  gStyle->SetStatX(0.9);
  gStyle->SetOptStat(10);
  gStyle->SetOptFit(0100); 
  gStyle->SetErrorX(0); 
  gPad->SetTicks();
  gPad->SetGrid();
  if (pEff_adc) {
    pEff_adc->Draw("E1");
    pEff_adc->SetMarkerStyle(4);
    pEff_adc->SetMarkerSize(0.7);
    pEff_adc->SetAxisRange(0.7, 1.,"Y");
    pEff_adc->GetYaxis()->SetTitleOffset(1.26);
    pEff_adc->Fit("pol0");
  }
   
  // ST Hit Efficiency histogram
  c1->cd(2);
  gStyle->SetOptStat(10);
  gPad->SetTicks();
  gPad->SetGrid();
  if (pEff) {
    pEff->Draw("E1");
    pEff->SetMarkerStyle(4);
    pEff->SetMarkerSize(0.7);
    pEff->SetAxisRange(0.7, 1.,"Y");
    pEff->GetYaxis()->SetTitleOffset(1.26);
    pEff->Fit("pol0");
  }

}
