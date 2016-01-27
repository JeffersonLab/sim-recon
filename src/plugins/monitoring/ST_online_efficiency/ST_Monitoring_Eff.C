// File: ST_Monitoring_Efficiency.C
// Last Modified: 08/27/2015
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms for efficiency studies
{
  // Define the directory that contains the histograms
   TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_efficiency");
   //TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("st_tracking");
  if(dir) dir->cd();
  // Grab 1D histograms 
  TH1D *h_N_trck_prd_All     = (TH1D*)dir->FindObjectAny("h_N_trck_prd_All");
  TH1D *h_N_recd_hit_All     = (TH1D*)dir->FindObjectAny("h_N_recd_hit_All");
  TH1D *h_ST_Eff_All = (TH1D*)dir->FindObjectAny("h_ST_Eff_All");
  TH1D *h_N_trck_prd_ss  = (TH1D*)dir->FindObjectAny("h_N_trck_prd_ss");
  TH1D *h_N_recd_hit_ss  = (TH1D*)dir->FindObjectAny("h_N_recd_hit_ss");
  TH1D *h_ST_Eff_ss  = (TH1D*)dir->FindObjectAny("h_ST_Eff_ss");
  TH1D *h_N_trck_prd_bs  = (TH1D*)dir->FindObjectAny("h_N_trck_prd_bs");
  TH1D *h_N_recd_hit_bs  = (TH1D*)dir->FindObjectAny("h_N_recd_hit_bs");
  TH1D *h_ST_Eff_bs  = (TH1D*)dir->FindObjectAny("h_ST_Eff_bs");
  TH1D *h_N_trck_prd_ns  = (TH1D*)dir->FindObjectAny("h_N_trck_prd_ns");
  TH1D *h_N_recd_hit_ns  = (TH1D*)dir->FindObjectAny("h_N_recd_hit_ns");
  TH1D *h_ST_Eff_ns  = (TH1D*)dir->FindObjectAny("h_ST_Eff_ns");

  // get Binomial errors in an efficiency plot
  // Whole ST efficiency
  h_N_recd_hit_All->Sumw2();
  h_N_trck_prd_All->Sumw2();
  h_ST_Eff_All->Sumw2();
  h_ST_Eff_All->Divide(h_N_recd_hit_All,h_N_trck_prd_All,1,1,"B");

  // SS ST efficiency
  h_N_recd_hit_ss->Sumw2();
  h_N_trck_prd_ss->Sumw2();
  h_ST_Eff_ss->Sumw2();
  h_ST_Eff_ss->Divide(h_N_recd_hit_ss,h_N_trck_prd_ss,1,1,"B");

  // BS ST efficiency
  h_N_recd_hit_bs->Sumw2();
  h_N_trck_prd_bs->Sumw2();
  h_ST_Eff_bs->Sumw2();
  h_ST_Eff_bs->Divide(h_N_recd_hit_bs,h_N_trck_prd_bs,1,1,"B");
  // NS ST efficiency
  h_N_recd_hit_ns->Sumw2();
  h_N_trck_prd_ns->Sumw2();
  h_ST_Eff_ns->Sumw2();
  h_ST_Eff_ns->Divide(h_N_recd_hit_ns,h_N_trck_prd_ns,1,1,"B");


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
  c1->Divide(2,2);
  // ST All Efficiency histogram
  c1->cd(1);
  gStyle->SetOptStat(0);
  gStyle->SetErrorX(0); 
  gPad->SetTicks();
  gPad->SetGrid();
  if (h_ST_Eff_All) {
   h_ST_Eff_All->Draw("E1");
   h_ST_Eff_All->SetMarkerStyle(21);
   h_ST_Eff_All->SetMarkerSize(1.5);
   h_ST_Eff_All->SetMarkerColor(4.0);
   h_ST_Eff_All->SetAxisRange(0.5, 1.1,"Y");
   h_ST_Eff_All->GetYaxis()->SetTitleOffset(1.26);
   h_ST_Eff_All->GetYaxis()->SetTitle("Efficiency");
   Double_t TWA_All=0;
   Double_t sumE_All=0;
   Double_t Final_All=0;
   for (int i = 0; i < 30; i++)
     {
       Double_t error_All =   h_ST_Eff_All->GetBinError(i+2);
       sumE_All = sumE_All + error_All;
       Double_t BC_All =   h_ST_Eff_All->GetBinContent(i+2);
       Double_t WA_All = BC_All*error_All;
       TWA_All = TWA_All + WA_All ;
       Double_t Final_All = TWA_All/sumE_All;
     } 
   TLine *lineAll = new TLine(1,Final_All,30,Final_All);
   lineAll->SetLineWidth(3);
   lineAll->SetLineColor(2);
   //Write the eff value on the histogram 
   char tFinal_All[20];
   char terror_All[20];
   sprintf(tFinal_All,"All Efficiency = %0.2f #pm %0.2f",Final_All*100,sumE_All*100/30);
   //sprintf(terror_All,"All Efficiency error (%) = %g",error_All*100);
   lineAll->Draw();
   TPaveLabel *pAll = new TPaveLabel(0.2,0.2,0.7,0.4,tFinal_All,"brNDC");
   pAll->Draw();
  }
  //SS Efficiency histogram
  c1->cd(2);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h_ST_Eff_ss) {
    h_ST_Eff_ss->Draw("E1");
    h_ST_Eff_ss->SetMarkerStyle(21);
    h_ST_Eff_ss->SetMarkerSize(1.5);
    h_ST_Eff_ss->SetMarkerColor(4.0);
    h_ST_Eff_ss->SetAxisRange(0.5, 1.1,"Y");
    h_ST_Eff_ss->GetYaxis()->SetTitleOffset(1.26);
    h_ST_Eff_ss->GetYaxis()->SetTitle("Efficiency");
    Double_t TWA_ss=0;
    Double_t sumE_ss=0;
    Double_t Final_ss=0; 
    for (int i = 0; i < 30; i++)
      {
	Double_t error_ss =    h_ST_Eff_ss->GetBinError(i+2);
	sumE_ss = sumE_ss+error_ss;
	Double_t BC_ss =    h_ST_Eff_ss->GetBinContent(i+2);
	Double_t WA_ss = BC_ss*error_ss;
	TWA_ss = TWA_ss + WA_ss ;
	Double_t Final_ss = TWA_ss/sumE_ss;
      } 
    TLine *liness = new TLine(1,Final_ss,30,Final_ss);
    liness->SetLineWidth(3);
    liness->SetLineColor(2);
    //Write the eff value on the histogram 
   char tFinal_ss[20];
   char terror_ss[20];
   sprintf(tFinal_ss,"SS Hit Efficiency  = %0.2f #pm %0.2f",Final_ss*100,sumE_ss*100/30);
   //   sprintf(terror_ss,"SS  Hit Efficiency error (%) = %g",error_ss*100);
   liness->Draw();
   TPaveLabel *pss = new TPaveLabel(0.2,0.2,0.7,0.4,tFinal_ss,"brNDC");
   pss->Draw();
  }
//BS Efficiency histogram
  c1->cd(3);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h_ST_Eff_bs) {
    h_ST_Eff_bs->Draw("E1");
    h_ST_Eff_bs->SetMarkerStyle(21);
    h_ST_Eff_bs->SetMarkerSize(1.5);
    h_ST_Eff_bs->SetMarkerColor(4.0);
    h_ST_Eff_bs->SetAxisRange(0.5, 1.1,"Y");
    h_ST_Eff_bs->GetYaxis()->SetTitleOffset(1.26);
    h_ST_Eff_bs->GetYaxis()->SetTitle("Efficiency");
    Double_t TWA_bs=0;
    Double_t sumE_bs=0;
    Double_t Final_bs=0; 
    for (int i = 0; i < 30; i++)
      {
	Double_t error_bs =    h_ST_Eff_bs->GetBinError(i+2);
	sumE_bs = sumE_bs+error_bs;
	Double_t BC_bs =    h_ST_Eff_bs->GetBinContent(i+2);
	Double_t WA_bs = BC_bs*error_bs;
	TWA_bs = TWA_bs + WA_bs ;
	Double_t Final_bs = TWA_bs/sumE_bs;
      } 
    TLine *linebs = new TLine(1,Final_bs,30,Final_bs);
    linebs->SetLineWidth(3);
    linebs->SetLineColor(2);
    //Write the eff value on the histogram 
   char tFinal_bs[20];
   char terror_bs[20];
   sprintf(tFinal_bs,"BS Hit Efficiency = %0.2f #pm %0.2f ",Final_bs*100,sumE_bs*100/30);
   linebs->Draw();
   TPaveLabel *pbs = new TPaveLabel(0.2,0.2,0.7,0.4,tFinal_bs,"brNDC");
   pbs->Draw();
  }
//NS Efficiency histogram
  c1->cd(4);
  gPad->SetTicks();
  gPad->SetGrid();
  if (h_ST_Eff_ns) {
    h_ST_Eff_ns->Draw("E1");
    h_ST_Eff_ns->SetMarkerStyle(21);
    h_ST_Eff_ns->SetMarkerSize(1.5);
    h_ST_Eff_ns->SetMarkerColor(4.0);
    h_ST_Eff_ns->SetAxisRange(0.5, 1.1,"Y");
    h_ST_Eff_ns->GetYaxis()->SetTitleOffset(1.26);
    h_ST_Eff_ns->GetYaxis()->SetTitle("Efficiency");
    Double_t TWA_ns=0;
    Double_t sumE_ns=0;
    Double_t Final_ns=0; 
    for (int i = 0; i < 30; i++)
      {
	Double_t error_ns =    h_ST_Eff_ns->GetBinError(i+2);
	sumE_ns = sumE_ns+error_ns;
	Double_t BC_ns =    h_ST_Eff_ns->GetBinContent(i+2);
	Double_t WA_ns = BC_ns*error_ns;
	TWA_ns = TWA_ns + WA_ns ;
	Double_t Final_ns = TWA_ns/sumE_ns;
      } 
    TLine *linens = new TLine(1,Final_ns,30,Final_ns);
    linens->SetLineWidth(3);
    linens->SetLineColor(2);
    //Write the eff value on the histogram 
   char tFinal_ns[20];
   char terror_ns[20];
   sprintf(tFinal_ns,"NS Hit Efficiency = %0.2f #pm %0.2f",Final_ns*100,sumE_ns*100/30);
   linens->Draw();
   TPaveLabel *pns = new TPaveLabel(0.2,0.2,0.7,0.4,tFinal_ns,"brNDC");
   pns->Draw();
   }
}
