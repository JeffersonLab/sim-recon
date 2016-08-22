// File: ST_Monitoring_2D.C
// Created: 05/21/2015
// Creator: Mahmoud Kamel, mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes
//void ST_Monitoring_2D()
void ST_Monitoring_Pid() {
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_tracking");
  //gDirectory->ls();
  if(dir) dir->cd();

  // Grab 2D histograms  
  TH2I *h2_dedx_P_mag       = (TH2I*)gDirectory->FindObjectAny("h2_dedx_P_mag");
  TH2I *h2_dedx_P_mag_postv = (TH2I*)gDirectory->FindObjectAny("h2_dedx_P_mag_postv");
  TH2I *h2_dedx_P_mag_negtv= (TH2I*)gDirectory->FindObjectAny("h2_dedx_P_mag_negtv");
  TH2I *h2_dedx__theta = (TH2I*)gDirectory->FindObjectAny("h2_dedx__theta");
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
  c1->Divide(2,1);
  gStyle->SetOptStat(10);
  gStyle->SetTitleFontSize(0.04);
  c1->cd(1);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogz();
  
  if(h2_dedx_P_mag_postv) h2_dedx_P_mag_postv->Draw("colz");
  
  h2_dedx_P_mag_postv->GetYaxis()->SetLabelSize(.02);
  h2_dedx_P_mag_postv->GetYaxis()->SetTitle("#frac{dE}{dx} (au)");
  h2_dedx_P_mag_postv->GetYaxis()->CenterTitle();
  h2_dedx_P_mag_postv->GetYaxis()->SetTitleOffset(1.25);

  h2_dedx_P_mag_postv->GetXaxis()->CenterTitle();
  h2_dedx_P_mag_postv->GetZaxis()->SetLabelFont(20);
  h2_dedx_P_mag_postv->GetZaxis()->SetLabelSize(0.02);
 

  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  gPad->SetLogz(); 
  if(h2_dedx_P_mag_negtv) h2_dedx_P_mag_negtv->Draw("colz");
  h2_dedx_P_mag_negtv->GetYaxis()->SetLabelSize(.02);
  h2_dedx_P_mag_negtv->GetYaxis()->SetTitle("#frac{dE}{dx} (au)");
  h2_dedx_P_mag_negtv->GetYaxis()->CenterTitle();
  h2_dedx_P_mag_negtv->GetXaxis()->CenterTitle();
  h2_dedx_P_mag_negtv->GetYaxis()->SetTitleOffset(1.25);
  h2_dedx_P_mag_negtv->GetZaxis()->SetLabelFont(20);
  h2_dedx_P_mag_negtv->GetZaxis()->SetLabelSize(0.02);
   
}
