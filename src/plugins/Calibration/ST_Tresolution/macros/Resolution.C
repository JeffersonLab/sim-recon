// File: st_tw_fits.C
// Last Modified: 11/10/2015
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms and tw automatic process.
#include "TH1.h"
#include <TH2I.h>
#include <TH1D.h>
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
#include <TPaletteAxis.h>
#include <stdio.h>
#include <stdint.h>
#include <fstream>
#include <TNtuple.h>
#include <TGraphErrors.h>
using namespace std;
using std::string;
// ***************** define constants and varibles*************************
const Int_t NCHANNELS          = 30;

// Double_t sudo_mpv_chan[NCHANNELS];
Double_t t_all_fit[NCHANNELS][3];
Double_t t_all_fit_err[NCHANNELS][3];

Double_t t_ss_fit[NCHANNELS][3];
Double_t t_ss_fit_err[NCHANNELS][3];

Double_t t_bs_fit[NCHANNELS][3];
Double_t t_bs_fit_err[NCHANNELS][3];

Double_t t_ns_fit[NCHANNELS][3];
Double_t t_ns_fit_err[NCHANNELS][3];

Double_t t_bn_fit[NCHANNELS][3];
Double_t t_bn_fit_err[NCHANNELS][3];

Double_t t_total_fit[NCHANNELS][3];
Double_t t_total_fit_err[NCHANNELS][3];

Double_t ex[NCHANNELS];
Float_t t_bin[NCHANNELS];
Double_t t_max[NCHANNELS];
Double_t t_min[NCHANNELS];
Double_t tbn_max[NCHANNELS];
Double_t tbn_min[NCHANNELS];
Double_t tns_max[NCHANNELS];
Double_t tns_min[NCHANNELS];
Double_t tbs_max[NCHANNELS];
Double_t tbs_min[NCHANNELS];
Double_t tss_max[NCHANNELS];
Double_t tss_min[NCHANNELS];
Double_t tall_max[NCHANNELS];
Double_t tall_min[NCHANNELS];
Float_t bins[210];
Float_t binsall[210];
Float_t binsss[210];
Float_t binsbs[210];
Float_t binsns[210];
Float_t binsbn[210];
//TNtuple *ntuple;
//TGraphErrors *gr;
//****   Declare fits *****************************
TF1 *t_vs_z_fit_chan[NCHANNELS];
TH1D *py[NCHANNELS];
// Declare canvas
TCanvas *PT_can[30];
TDirectory* TopDirectory;
TH2I* h2_total;
Double_t fitf_pp(Double_t *x, Double_t *par)
{
  Double_t fitval_pp = par[0] + par[1]*x[0];//(TMath::Power(x[0]/adc_thresh_calc, par[2]));
  return fitval_pp;
}

void Resolution(char const*input_filename)
//void st_tw_fits()
{
  TFile *df = new TFile(input_filename);
  std::ofstream st_time_res;
  st_time_res.open ("st_time_res.txt", std::ofstream::out);

  //*****************************************************************
  //*********** Grab the histograms from the root file **************
  //*****************************************************************
 
  for (unsigned int j = 0; j < 30; j++)
    {
      //Create the canvas
      PT_can[j] = new TCanvas( Form("PT_can_%i",j+1), Form("PT_can_%i",j+1), 800, 450);
      PT_can[j]->Divide(2, 1);
	 // The top directory 
      TopDirectory = (TDirectory*) df->FindObjectAny("ST_Tresolution"); 
      TopDirectory->cd(); 
      // Grab the histograms
      //Corrected
      char* total = Form("h2_CorrectedTime_z_%i",j+1);
      h2_total = (TH2I*) TopDirectory->FindObjectAny(total);
      char* pj = Form("py_%i",j+1);
      TH1D* pytotal = (TH1D*)pj;
      
      cout << "==================================================" << endl;
      cout << "Processing Channel " << j+1 << endl;
      cout << "==================================================" << endl;
      ex[j]=0.0;
      
      
      //*****************************************************************
      //***********  Plot Total corrected time from each division *******
      //*****************************************************************
      PT_can[j]->cd(1);
      gStyle->SetOptStat(0);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h2_total->Draw("colz");
      h2_total->GetZaxis()->SetRangeUser(0,10000);
      h2_total->GetYaxis()->SetTitle("ST Corrected Time (ns)");
      h2_total->GetYaxis()->SetRangeUser(-5.0,5.0);
      gPad->SetLogz();
      PT_can[j]->cd(2);
      gStyle->SetOptFit(11);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      
      pytotal =  h2_total->ProjectionY();    
      t_max[j]=  pytotal->GetMaximum();
      t_min[j]=  pytotal->GetMinimum();
      // for (unsigned int i = 0; i < 200; i++)
      // 	{bins[i]= pytotal->GetBinContent(i);
      // 	  if (bins[i] < (0.15*t_max[j]))
      // 	    {	
      // 	      pytotal->SetBinContent(i,0);
      // 	      pytotal->Draw("][");
      // 	      pytotal->GetXaxis()->SetRangeUser(-1.5,1.5);
      // 	      pytotal->GetXaxis()->SetTitle("ST Corrected Time (ns)");
      // 	    } 
      // 	}
      pytotal->GetXaxis()->SetRangeUser(-5.0,5.0);
      pytotal->Fit("gaus","","",-0.6,0.6);
      TF1 *gaus = pytotal->GetFunction("gaus");
      gPad->Update();
       PT_can[j]->Print(Form("ST_Sector_%i.png",j+1));
      t_total_fit[j][2] = gaus->GetParameter(2)*1000;
      t_total_fit_err[j][2] = gaus->GetParError(2)*1000;
      
      st_time_res << j+1 <<" " << ex[j] << "   "<< t_total_fit[j][2] << "  "<< t_total_fit_err[j][2]<< endl;
      
    }
  st_time_res.close();
  ifstream in;
  in.open(Form("st_time_res.txt"));
  float sector,ser,tall,taller;
  TNtuple *ntuple = new TNtuple("ntuple","data from text file","sector:ser:tall:taller");
  ntuple->ReadFile("st_time_res.txt");
  ntuple->Write();
  in.close();
  //Create the canvas
  TCanvas *Time_can = new TCanvas( "Time_can", "Time_can", 800, 600);
  Time_can->SetFillColor(41);
  gStyle->SetOptFit(011);
  ntuple->Draw("sector:tall:taller:ser");
  TGraphErrors *gr=new TGraphErrors(ntuple->GetSelectedRows(),ntuple->GetV1(),ntuple->GetV2(),ntuple->GetV4(),ntuple->GetV3());
  gr->Draw("AP");
  gr->SetTitle("Time Resolution After Applying Propagation Time Corrections");
  gr->GetXaxis()->SetTitle("Sector Number");
  gr->GetYaxis()->SetTitle("Time Resolution (ps)");
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(3);
  gr->SetMarkerColor(4);
  gr->Fit("pol0");
  gr->SetMinimum(150);
  gr->SetMaximum(450);
  TLine *line = new TLine(0,350,33,350);
  line->SetLineWidth(4);
  line->SetLineColor(28);
  line->SetLineStyle(9);

  line->Draw();
  Time_can->Print(Form("ST_Tresolution.png"));
}

