// File: st_tw_fits.C
// Last Modified: 11/10/2015
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms and tw automatic process.
#include "TH1.h"
#include <TH2I.h>
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TAxis.h"
#include <TDirectory.h>
#include <TLine.h>
#include <TGraph.h>
#include <TPaveLabel.h>
#include <TPaletteAxis.h>
#include <stdio.h>
#include <stdint.h>
#include <fstream>

using namespace std;
using std::string;
// ***************** define constants and varibles*************************
const Int_t NCHANNELS          = 30;
const Double_t tdc_thresh_mV   = 50.0;
const Double_t tdc_gain_factor = 5.0;
const Double_t adc_max_chan    = 4096.0;
const Double_t adc_max_mV      = 2000.0;
const Double_t adc_thresh_calc = (tdc_thresh_mV/tdc_gain_factor)*(adc_max_chan/adc_max_mV);
Double_t sudo_mpv_chan[NCHANNELS];
Double_t t_pp_fit_params[NCHANNELS][3];
Double_t t_pp_fit_params_err[NCHANNELS][3];
Double_t stt_fit_params[NCHANNELS][3];
Double_t stt_fit_params_err[NCHANNELS][3];
Double_t sigma[NCHANNELS];
Double_t Err[NCHANNELS];
Double_t x[NCHANNELS];
Double_t ex[NCHANNELS];
//****   Declare fits *****************************
TF1 *t_fit_chan[NCHANNELS];
TGraph *gr[NCHANNELS];
// Declare canvas
TCanvas *TW_can[30];
TCanvas *Time_can;
Double_t fitf_pp(Double_t *x, Double_t *par)
{
  Double_t fitval_pp = par[0] + par[1]*(TMath::Power(x[0]/adc_thresh_calc, par[2]));
  return fitval_pp;
}

void st_tw_resols(char const*input_filename)

{
  TFile *df = new TFile(input_filename);
  std::ofstream Time_resol;
  Time_resol.open ("SC_time_resol.txt", std::ofstream::out);

  //*****************************************************************
  //*********** Grab the histograms from the root file **************
  //*****************************************************************
  for (unsigned int j = 0; j < 30; j++)
    {
        t_fit_chan[j] = new TF1(Form("t_fit_chan_%i", j+1), "gaus", -2.0,2.0);
        //t_fit_chan[j] = new TF1(Form("t_fit_chan_%i", j), "[0]*exp(-0.5*((x-[1])/[2])**2)", -2.0,2.0);
        t_fit_chan[j]->SetParameters(1000,0.,0.3);
	x[j]=j+1;
	ex[j]=0.0;
      //Create the canvas
      TW_can[j] = new TCanvas( Form("TW_can_%i",j+1), "TW_can", 800, 600);
      TW_can[j]->Divide(2, 2);
	  
      // Grab the histograms 
      char* stt = Form("h_stt_chan_%i",j+1);
      TH1I* h_st = (TH1I*) df->Get(stt);
      char* stt_corr = Form("h1_st_corr_time_%i",j+1);
      TH1I* h_st_cr = (TH1I*) df->Get(stt_corr);
      char* h2_st_corr_vs_pp = Form("h2_st_corr_vs_pp_%i",j+1);
      TH2I* h2_st_cr_pp = (TH2I*) df->Get(h2_st_corr_vs_pp);

      //*****************************************************************
      //*********** Plot 1D stt histos                 ******************
      //*****************************************************************
      TW_can[j]->cd(1);
      gStyle->SetOptStat(10);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h_st->Draw();
      h_st->GetXaxis()->SetRangeUser(-2.,2.);
      //*****************************************************************
      //*********** Plot 1D  histos for the Corrected stt ***************
      //*****************************************************************
      TW_can[j]->cd(2);
      gStyle->SetOptStat(10);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h_st_cr->Draw();
      h_st_cr->GetXaxis()->SetRangeUser(-2.,2.);
      h_st_cr->Fit(t_fit_chan[j]);
      sigma[j] =  t_fit_chan[j]->GetParameter(2);
      Err[j] =  t_fit_chan[j]->GetParError(2);
      Time_resol << x[j] << "\t" <<sigma[j]*1000 << "\t" << ex[j] << "\t" << Err[j]*1000 << endl;

      //*****************************************************************
      //*********** Plot 2D  histos the Corrected stt vs pp *************
      //*****************************************************************
      TW_can[j]->cd(3);
      gStyle->SetOptStat(10);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h2_st_cr_pp->Draw("colz");
      h2_st_cr_pp->GetXaxis()->SetRangeUser(120.,3500.);
      h2_st_cr_pp->GetYaxis()->SetRangeUser(-2.0,2.0);
    }
  //*****************************************************************
  //*********** Plot 1D  histos of the time resolution  *************
  //*****************************************************************
  Time_resol.close();
  ifstream in;
  in.open(Form("SC_time_resol.txt"));
  float sector, t, es ,e;
  TNtuple *ntuple = new TNtuple("ntuple","data from text file","sector:t:es:e");
  ntuple->ReadFile("SC_time_resol.txt");
  ntuple->Write();
  in.close();
  //Create the canvas
  Time_can = new TCanvas( "Time_can", "Time_can", 800, 600);
  Time_can->SetFillColor(41);
  gStyle->SetOptFit(0111);
  ntuple->Draw("sector:t:es:e");
  TGraphErrors *gr=new TGraphErrors(ntuple->GetSelectedRows(),ntuple->GetV1(),ntuple->GetV2(),ntuple->GetV3(),ntuple->GetV4());
  gr->Draw("AP");
  gr->SetTitle("ST Self Time Resolution Vs Sector Number");
  gr->GetXaxis()->SetTitle("Sector Number");
  gr->GetYaxis()->SetTitle("ST Self Time Resolution (ps)");
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(3);
  gr->SetMarkerColor(4);
  gr->Fit("pol0");
 
  
}

