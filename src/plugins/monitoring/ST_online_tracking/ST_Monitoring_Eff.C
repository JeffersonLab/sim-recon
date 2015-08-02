// File: ST_Monitoring_Efficiency.C
// Last Modified: 05/27/2015
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms for online monitoring purposes
{
  // Define the directory that contains the histograms
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_tracking");
  if(dir) dir->cd();
  // Grab 1D histograms 
  TH1D *MacropEff     = (TH1D*)dir->FindObjectAny("MacropEff");
  TH1D *MacropEff_adc     = (TH1D*)dir->FindObjectAny("MacropEff_adc");
  TH1D *h_phi_sec_pred_hit_cntr = (TH1D*)dir->FindObjectAny("h_phi_sec_pred_hit_cntr");
  TH1D *h_phi_sec_hit_cntr  = (TH1D*)dir->FindObjectAny("h_phi_sec_hit_cntr");
  TH1D *h_phi_sec_adc_cntr = (TH1D*)dir->FindObjectAny("h_phi_sec_adc_cntr");

  // get Binomial errors in an efficiency plot
  // hit object efficiency
  h_phi_sec_hit_cntr->Sumw2();
  h_phi_sec_pred_hit_cntr->Sumw2();
  MacropEff->Sumw2();
  MacropEff->Divide(h_phi_sec_hit_cntr,h_phi_sec_pred_hit_cntr,1,1,"B");

  // adc efficiency
  h_phi_sec_adc_cntr->Sumw2();
  MacropEff_adc->Sumw2();
  MacropEff_adc->Divide(h_phi_sec_adc_cntr,h_phi_sec_pred_hit_cntr,1,1,"B");

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
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0); 
  gPad->SetTicks();
  gPad->SetGrid();
  if (MacropEff_adc) {
   MacropEff_adc->Draw("E1");
   MacropEff_adc->SetMarkerStyle(21);
   MacropEff_adc->SetMarkerSize(1.5);
   MacropEff_adc->SetMarkerColor(4.0);
   MacropEff_adc->SetAxisRange(0., 1.,"Y");
   MacropEff_adc->GetYaxis()->SetTitleOffset(1.26);
   Double_t TWA_adc=0;
   Double_t sumE_adc=0;
   Double_t Final_adc=0;
   for (int i = 0; i < 30; i++)
     {
       Double_t error_adc =   MacropEff_adc->GetBinError(i+2);
       sumE_adc = sumE_adc + error_adc;
       Double_t BC_adc =   MacropEff_adc->GetBinContent(i+2);
       Double_t WA_adc = BC_adc*error_adc;
       TWA_adc = TWA_adc + WA_adc ;
       Double_t Final_adc = TWA_adc/sumE_adc;
       // cout << "error_adc =  "<< error_adc << endl;
       // cout << "BC_adc =  "<< BC_adc << endl;
       // cout << "Final_adc =  "<< Final_adc << endl;
     } 
   TLine *line = new TLine(1,Final_adc,30,Final_adc);
   line->SetLineWidth(3);
   line->SetLineColor(2);
   //Write the eff value on the histogram 
   char tFinal_adc[40];
   char terror_adc[40];
   sprintf(tFinal_adc,"ADC Efficiency (%) = %g",Final_adc*100);
   sprintf(terror_adc,"ADC Efficiency error (%) = %g",error_adc*100);
   line->Draw();
   TPaveLabel *p = new TPaveLabel(0.3,0.6,0.7,0.7,tFinal_adc,"brNDC");
   p->Draw();
   TPaveLabel *p1 = new TPaveLabel(0.3,0.5,0.7,0.6,terror_adc,"brNDC");
   p1->Draw();
  }
  // ST Hit Efficiency histogram
  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  if (MacropEff) {
    MacropEff->Draw("E1");
    MacropEff->SetMarkerStyle(21);
    MacropEff->SetMarkerSize(1.5);
    MacropEff->SetMarkerColor(4.0);
    MacropEff->SetAxisRange(0., 1.,"Y");
    MacropEff->GetYaxis()->SetTitleOffset(1.26);
    Double_t TWA=0;
    Double_t sumE=0;
    Double_t Final=0; 
    for (int i = 0; i < 30; i++)
      {
	Double_t error =   MacropEff->GetBinError(i+2);
	sumE = sumE+error;
	Double_t BC =   MacropEff->GetBinContent(i+2);
	Double_t WA = BC*error;
	TWA = TWA + WA ;
	Double_t Final = TWA/sumE;
	//	cout << "error =  "<< error << endl;
	//	cout << "BC =  "<< BC << endl;
	//	cout << "Final =  "<< Final << endl;
      } 
    TLine *line = new TLine(1,Final,30,Final);
    line->SetLineWidth(3);
    line->SetLineColor(2);
    //Write the eff value on the histogram 
   char tFinal[40];
   char terror[40];
   sprintf(tFinal,"Hit Efficiency (%) = %g",Final*100);
   sprintf(terror,"Hit Efficiency error (%) = %g",error*100);
   line->Draw();
   TPaveLabel *p = new TPaveLabel(0.3,0.6,0.7,0.7,tFinal,"brNDC");
   p->Draw();
   TPaveLabel *p1 = new TPaveLabel(0.3,0.5,0.7,0.6,terror,"brNDC");
   p1->Draw();
  }
}
