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
//****   Declare fits *****************************
TF1 *t_vs_pp_fit_chan[NCHANNELS];
// Declare canvas
TCanvas *TW_can[30];

Double_t fitf_pp(Double_t *x, Double_t *par)
{
  Double_t fitval_pp = par[0] + par[1]*(TMath::Power(x[0]/adc_thresh_calc, par[2]));
  return fitval_pp;
}

void st_tw_fits(char const*input_filename)
//void st_tw_fits()
{
  // df = new TFile("/lustre/expphy/work/halld/home/mkamel/st_tw_corr_auto/hd_root_pass0.root");
  // df = new TFile("/lustre/expphy/work/halld/home/mkamel/st_tw_corr_auto/hd_root.root");
   TFile *df = new TFile(input_filename);
  std::ofstream results;
  results.open ("st_timewalks.txt", std::ofstream::out);
  //*****************************************************************
  //*********** Grab the histograms from the root file **************
  //*****************************************************************
  for (unsigned int j = 0; j < 30; j++)
    {
      //Create the canvas
      TW_can[j] = new TCanvas( Form("TW_can_%i",j+1), "TW_can", 800, 450);
      TW_can[j]->Divide(2, 1);
	  
      // Grab the histograms 
      char* sttpp = Form("stt_vs_pp_chan_%i",j+1);
      TH2I* h2_tpp = (TH2I*) df->Get(sttpp);
      char* pp = Form("pp_chan_%i",j+1);
      TH1I* h_pp = (TH1I*) df->Get(pp);
      //*****************************************************************
      //***********      Plot stt vs pp                ******************
      //*****************************************************************
      TW_can[j]->cd(1);
      gStyle->SetOptStat(10);
      gPad->SetTicks();
      gPad->SetGrid();
    
      h2_tpp->Draw("colz"); 
      h2_tpp->GetZaxis()->SetRangeUser(0,300);
      h2_tpp->GetXaxis()->SetRangeUser(120.,3500.);
      h2_tpp->GetYaxis()->SetRangeUser(-2.0,2.0);
      TPaletteAxis *palette = new TPaletteAxis(3501.0,-2.0,3599.0,2.0,h2_tpp);
      palette->SetLabelOffset(0.002);
      palette->SetLabelSize(0.05);
      h2_tpp->GetListOfFunctions()->Add(palette,"br");
      //*****************************************************************
      //***********          Do the fit                ******************
      //*****************************************************************
      cout << "==================================================" << endl;
      cout << "Processing Channel " << j+1 << endl;
      cout << "==================================================" << endl;
      t_vs_pp_fit_chan[j] = new TF1(Form("t_vs_pp_fit_chan_%i", j+1), fitf_pp, 120.0,3500.0, 3);   
      t_vs_pp_fit_chan[j]->SetParameters(-1.0, 7.0, -0.5);
     
      h2_tpp->Fit(Form("t_vs_pp_fit_chan_%i", j+1));
      
      t_pp_fit_params[j][0] =  t_vs_pp_fit_chan[j]->GetParameter(0);
      t_pp_fit_params_err[j][0] =  t_vs_pp_fit_chan[j]->GetParError(0);
      t_pp_fit_params[j][1] =  t_vs_pp_fit_chan[j]->GetParameter(1);
      t_pp_fit_params_err[j][1] =  t_vs_pp_fit_chan[j]->GetParError(1);
      t_pp_fit_params[j][2] =  t_vs_pp_fit_chan[j]->GetParameter(2);
      t_pp_fit_params_err[j][2] =  t_vs_pp_fit_chan[j]->GetParError(2);
      gPad->Update();
      //*****************************************************************
      //*********** Plot 1D pp histos                  ******************
      //*****************************************************************
      TW_can[j]->cd(2);
      gStyle->SetOptStat(10);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h_pp->Draw(); 
      h_pp->GetXaxis()->SetRangeUser(120.0,3500.);
      sudo_mpv_chan[j] = (Double_t) h_pp->GetMaximumBin();
      //*****************************************************************
      //*********** Write the Results                  ******************
      //*****************************************************************
      results << t_pp_fit_params[j][0] << "\t" << t_pp_fit_params[j][1] << "\t" << "\t" << t_pp_fit_params[j][2] << "\t" << adc_thresh_calc  << "000" << "\t"  <<  sudo_mpv_chan[j] << ".0000"<< endl;
      TW_can[j]->Print(Form("stt_tw_plot_%i.png",j+1));
    }
  results.close();
}

