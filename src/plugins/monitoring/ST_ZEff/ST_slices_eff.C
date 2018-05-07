// File: ST_Monitoring_Efficiency.C
// Last Modified: 09/15/2015
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms for efficiency studies
#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <TDirectory.h>
#include <TLine.h>
#include <TPaveLabel.h>
// Declare data files
TFile *df;
// Define some constants
const int Nof_ss_intervals = 8;
const int Nof_bs_intervals = 4;
const int Nof_ns_intervals = 8;
const int NChanel = 30;
// Declare canvas
//TCanvas *SS_can;
void ST_slices_eff()
{
  double error_ss[Nof_ss_intervals];
  double sumE_ss[Nof_ss_intervals];
  double BC_ss[Nof_ss_intervals];
  double WA_ss[Nof_ss_intervals];
  double TWA_ss[Nof_ss_intervals];
  double Final_ss[Nof_ss_intervals];

  double error_bs[Nof_bs_intervals];
  double sumE_bs[Nof_bs_intervals];
  double BC_bs[Nof_bs_intervals];
  double WA_bs[Nof_bs_intervals];
  double TWA_bs[Nof_bs_intervals];
  double Final_bs[Nof_bs_intervals];

  double error_ns[Nof_ns_intervals];
  double sumE_ns[Nof_ns_intervals];
  double BC_ns[Nof_ns_intervals];
  double WA_ns[Nof_ns_intervals];
  double TWA_ns[Nof_ns_intervals];
  double Final_ns[Nof_ns_intervals];

  //  df = new TFile("./OUTPUT_JOBS/all_runs.root");
  //df = new TFile("hd_root_SIM600K.root");
  //   df = new TFile("./hists/hd_root_10runs.root");
  // df = new TFile("Five_runs.root");
  df = new TFile("hd_root.root");
 
  //df = new TFile("./hd_root_Sim1_1Mevent.root");
  TDirectory *dir = (TDirectory*)gDirectory->FindObjectAny("./st_efficiency");
  if(dir) dir->cd();
  //TDirectory *dir1 = (TDirectory*)gDirectory->FindObjectAny("st_eff_nebrg");
  //if(dir) dir1->cd();
  //*****************************************************************
  //***********          SS slices efficiency      ******************
  //*****************************************************************
  //Create the canvas
  SS_can = new TCanvas("SS_can", "SS_can", 800, 800);
  SS_can->Divide(3, 3);
  for (unsigned int i = 1; i < Nof_ss_intervals+1; i++)
    {
      char* prd_z_ss = Form("h_N_trck_prd_z_ss_%i",i);
      char* hit_z_ss = Form("h_N_recd_hit_z_ss_%i",i);
      char* eff_z_ss = Form("h_N_trck_prd_z_ss_eff_%i",i);
      char* ACC_z_ss = Form("h_N_recd_hit_z_ss_ACC_%i",i);

      TH1D *hp = (TH1D*) gDirectory->FindObjectAny(prd_z_ss);
      TH1D *hh = (TH1D*) gDirectory->FindObjectAny(hit_z_ss);
      TH1D *he = (TH1D*) gDirectory->FindObjectAny(eff_z_ss);
      TH1D *hA = (TH1D*) gDirectory->FindObjectAny(ACC_z_ss);
    
    
      hh->Sumw2();
      hp->Sumw2();
      he->Sumw2();
      hA->Sumw2();
      TH1D * hh_final = (TH1D*) hh->Clone();     
      hh_final->Add(hA,-1.0);    
      he->Divide(hh_final,hp,1.,1.,"B");       
      
      SS_can->cd(i);
      gStyle->SetOptStat(0);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      
      he->Draw("E1");
      he->SetMarkerStyle(21);
      he->SetMarkerSize(1.1);
      he->SetMarkerColor(4.0);
      he->SetAxisRange(0.5, 1.0,"Y");
      he->GetYaxis()->SetTitleOffset(1.26);
      
      for (int j = 0; j < 30; j++)
	{
	  error_ss[i-1]= he->GetBinError(j+2);
	  sumE_ss[i-1] = sumE_ss[i-1] + error_ss[i-1];
	  BC_ss[i-1]   = he->GetBinContent(j+2);
	  WA_ss[i-1]   = BC_ss[i-1]*error_ss[i-1];
	  TWA_ss[i-1]  = TWA_ss[i-1]+ WA_ss[i-1];
	  Final_ss[i-1]= TWA_ss[i-1]/sumE_ss[i-1];
	}  
      TLine *line = new TLine(1,Final_ss[i-1],30,Final_ss[i-1]);
      line->SetLineWidth(3);
      line->SetLineColor(2);
      //Write the eff value on the histogram 
      char tFinal_ss[20][Nof_ss_intervals];
      char terror_ss[20][Nof_ss_intervals];
      sprintf(tFinal_ss[i-1],"Straight Section Efficiency = %0.2f #pm %0.2f",Final_ss[i-1]*100,sumE_ss[i-1]*100/30);
      line->Draw();
      TPaveLabel *p[Nof_ss_intervals];
      p[i-1] = new TPaveLabel(0.2,0.2,0.7,0.4,tFinal_ss[i-1],"brNDC");
      p[i-1]->Draw();
      
      memset(error_ss, 0, sizeof(error_ss));
      memset(sumE_ss, 0, sizeof(sumE_ss));
      memset(BC_ss, 0, sizeof(BC_ss));
      memset(WA_ss, 0, sizeof(WA_ss));
      memset(TWA_ss, 0, sizeof(TWA_ss));
      memset(Final_ss, 0, sizeof(Final_ss));
      
    }
  //*****************************************************************
  //***********          BS slices efficiency      ******************
  //*****************************************************************
  //Create the canvas
  BS_can = new TCanvas("BS_can", "BS_can", 800, 800);
  BS_can->Divide(2, 2);
  for (unsigned int i = 1; i < Nof_bs_intervals+1; i++)
    {
      char* prd_z_bs = Form("h_N_trck_prd_z_bs_%i",i);
      char* hit_z_bs = Form("h_N_recd_hit_z_bs_%i",i);
      char* eff_z_bs = Form("h_N_trck_prd_z_bs_eff_%i",i);
      char* ACC_z_bs = Form("h_N_recd_hit_z_bs_ACC_%i",i);

      TH1D *hp_bs = (TH1D*) gDirectory->FindObjectAny(prd_z_bs);
      TH1D *hh_bs = (TH1D*) gDirectory->FindObjectAny(hit_z_bs);
      TH1D *he_bs = (TH1D*) gDirectory->FindObjectAny(eff_z_bs);
      TH1D *hA_bs = (TH1D*) gDirectory->FindObjectAny(ACC_z_bs);
    
      hh_bs->Sumw2();
      hp_bs->Sumw2();
      he_bs->Sumw2();
      hA_bs->Sumw2();
      TH1D * hh_bs_final = (TH1D*) hh_bs->Clone();     
      hh_bs_final->Add(hA_bs,-1.0);    

      he_bs->Divide(hh_bs_final,hp_bs,1.,1.,"B");       
      
      BS_can->cd(i);
      gStyle->SetOptStat(0);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      
      he_bs->Draw("E1");
      he_bs->SetMarkerStyle(21);
      he_bs->SetMarkerSize(1.1);
      he_bs->SetMarkerColor(4.0);
      he_bs->SetAxisRange(0.9, 1.0,"Y");
      he_bs->GetYaxis()->SetTitleOffset(1.26);
      
      for (int j = 0; j < 30; j++)
	{
	  error_bs[i-1]= he_bs->GetBinError(j+2);
	  sumE_bs[i-1] = sumE_bs[i-1] + error_bs[i-1];
	  BC_bs[i-1]   = he_bs->GetBinContent(j+2);
	  WA_bs[i-1]   = BC_bs[i-1]*error_bs[i-1];
	  TWA_bs[i-1]  = TWA_bs[i-1]+ WA_bs[i-1];
	  Final_bs[i-1]= TWA_bs[i-1]/sumE_bs[i-1];
	}  
      TLine *line = new TLine(1,Final_bs[i-1],30,Final_bs[i-1]);
      line->SetLineWidth(3);
      line->SetLineColor(2);
      //Write the eff value on the histogram 
      char tFinal_bs[20][Nof_bs_intervals];
      char terror_bs[20][Nof_bs_intervals];
      sprintf(tFinal_bs[i-1],"Bend Section Efficiency = %0.2f #pm %0.2f",Final_bs[i-1]*100,sumE_bs[i-1]*100/30);
      line->Draw();
      TPaveLabel *p[Nof_bs_intervals];
      p[i-1] = new TPaveLabel(0.2,0.2,0.7,0.4,tFinal_bs[i-1],"brNDC");
      p[i-1]->Draw();
      
      memset(error_bs, 0, sizeof(error_bs));
      memset(sumE_bs, 0, sizeof(sumE_bs));
      memset(BC_bs, 0, sizeof(BC_bs));
      memset(WA_bs, 0, sizeof(WA_bs));
      memset(TWA_bs, 0, sizeof(TWA_bs));
      memset(Final_bs, 0, sizeof(Final_bs));
    }
  //*****************************************************************
  //***********          NS slices efficiency      ******************
  //*****************************************************************
  //Create the canvas
  NS_can = new TCanvas("NS_can", "NS_can", 800, 800);
  NS_can->Divide(3, 3);
  for (unsigned int i = 1; i < Nof_ns_intervals+1; i++)
    {
      char* prd_z_ns = Form("h_N_trck_prd_z_ns_%i",i);
      char* hit_z_ns = Form("h_N_recd_hit_z_ns_%i",i);
      char* eff_z_ns = Form("h_N_trck_prd_z_ns_eff_%i",i);
      char* ACC_z_ns = Form("h_N_recd_hit_z_ns_ACC_%i",i);

      TH1D *hp_ns = (TH1D*) gDirectory->FindObjectAny(prd_z_ns);
      TH1D *hh_ns = (TH1D*) gDirectory->FindObjectAny(hit_z_ns);
      TH1D *he_ns = (TH1D*) gDirectory->FindObjectAny(eff_z_ns);
      TH1D *hA_ns = (TH1D*) gDirectory->FindObjectAny(ACC_z_ns);
    
      hh_ns->Sumw2();
      hp_ns->Sumw2();
      he_ns->Sumw2();
      hA_ns->Sumw2();
      TH1D * hh_ns_final = (TH1D*) hh_ns->Clone();     
      hh_ns_final->Add(hA_ns,-1.0);    


      he_ns->Divide(hh_ns_final,hp_ns,1.,1.,"B");       
      
      NS_can->cd(i);
      gStyle->SetOptStat(0);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      
      he_ns->Draw("E1");
      he_ns->SetMarkerStyle(21);
      he_ns->SetMarkerSize(1.1);
      he_ns->SetMarkerColor(4.0);
      he_ns->SetAxisRange(0.9, 1.0,"Y");
      he_ns->GetYaxis()->SetTitleOffset(1.26);
      
      for (int j = 0; j < 30; j++)
	{
	  error_ns[i-1]= he_ns->GetBinError(j+2);
	  sumE_ns[i-1] = sumE_ns[i-1] + error_ns[i-1];
	  BC_ns[i-1]   = he_ns->GetBinContent(j+2);
	  WA_ns[i-1]   = BC_ns[i-1]*error_ns[i-1];
	  TWA_ns[i-1]  = TWA_ns[i-1]+ WA_ns[i-1];
	  Final_ns[i-1]= TWA_ns[i-1]/sumE_ns[i-1];
	}  
      TLine *line = new TLine(1,Final_ns[i-1],30,Final_ns[i-1]);
      line->SetLineWidth(3);
      line->SetLineColor(2);
      //Write the eff value on the histogram 
      char tFinal_ns[20][Nof_ns_intervals];
      char terror_ns[20][Nof_ns_intervals];
      sprintf(tFinal_ns[i-1],"Nose Section Efficiency = %0.2f #pm %0.2f",Final_ns[i-1]*100,sumE_ns[i-1]*100/30);
      line->Draw();
      TPaveLabel *p[Nof_ns_intervals];
      p[i-1] = new TPaveLabel(0.2,0.2,0.7,0.4,tFinal_ns[i-1],"brNDC");
      p[i-1]->Draw();
      
      memset(error_ns, 0, sizeof(error_ns));
      memset(sumE_ns, 0, sizeof(sumE_ns));
      memset(BC_ns, 0, sizeof(BC_ns));
      memset(WA_ns, 0, sizeof(WA_ns));
      memset(TWA_ns, 0, sizeof(TWA_ns));
      memset(Final_ns, 0, sizeof(Final_ns));
      
    }
}

