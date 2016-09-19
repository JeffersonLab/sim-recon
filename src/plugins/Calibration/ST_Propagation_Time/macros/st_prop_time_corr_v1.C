// File: st_tw_fits.C
// Last Modified: 01/12/2016
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms for propagation time corrections
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
using namespace std;
using std::string;
// ***************** define constants and varibles*************************
const Int_t NCHANNELS          = 30;
// Double_t sudo_mpv_chan[NCHANNELS];
Double_t t_ss_fit[NCHANNELS][3];
Double_t t_ss_fit_err[NCHANNELS][3];
Double_t t_bs_fit[NCHANNELS][3];
Double_t t_bs_fit_err[NCHANNELS][3];
Double_t t_ns_fit[NCHANNELS][3];
Double_t t_ns_fit_err[NCHANNELS][3];
Double_t t_total_fit[NCHANNELS][3];
Double_t t_total_fit_err[NCHANNELS][3];
Double_t ex[NCHANNELS];

Double_t t_bin[NCHANNELS];
Double_t t_max[NCHANNELS];
Double_t t_min[NCHANNELS];
Double_t tns_max[NCHANNELS];
Double_t tns_min[NCHANNELS];
Double_t tbs_max[NCHANNELS];
Double_t tbs_min[NCHANNELS];
Double_t tss_max[NCHANNELS];
Double_t tss_min[NCHANNELS];

Float_t bins[210];
Float_t binsss[210];
Float_t binsbs[210];
Float_t binsns[210];
//******************************** Cut Parameters **********************
Double_t low_cut = 0.00 ; // set bin content to zero if bin content is below low_cut
//****   Declare 1D Histogram ******************************************
TH1D *py[NCHANNELS];
// Declare canvas
TCanvas *PT_can[30];
TDirectory* TopDirectory;
TH2I* h2_ss;
TH2I* h2_bs;
TH2I* h2_ns;
TH2I* h2_total;

void st_prop_time_corr_v1(char const* input_filename)
//void st_tw_fits()
{
  TFile *df = new TFile(input_filename);
  std::ofstream st_prop_time_constants;
  st_prop_time_constants.open ("st_prop_timeCorr.txt", std::ofstream::out);
  //*****************************************************************
  //*********** Grab the histograms from the root file **************
  //*****************************************************************
  for (unsigned int j = 0; j < 30; j++)
    {
      //Create the canvas
      PT_can[j] = new TCanvas( Form("PT_can_%i",j+1), Form("PT_can_%i",j+1), 800, 450);
      PT_can[j]->Divide(2, 2);
      // The top directory 
      TopDirectory = (TDirectory*) df->FindObjectAny("ST_Propagation_Time"); 
      TopDirectory->cd();
      // Grab the histograms 
      char* ss = Form("h2_PropTimeCorr_z_SS_chan_%i",j+1);
      h2_ss = (TH2I*) TopDirectory->FindObjectAny(ss);
      char* bs = Form("h2_PropTimeCorr_z_BS_chan_%i",j+1);
      h2_bs = (TH2I*) TopDirectory->FindObjectAny(bs);
      char* ns = Form("h2_PropTimeCorr_z_NS_chan_%i",j+1);
      h2_ns = (TH2I*) TopDirectory->FindObjectAny(ns);
      char* total = Form("h2_CorrectedTime_z_%i",j+1);
      h2_total = (TH2I*) TopDirectory->FindObjectAny(total);      
      //*****************************************************************
      //***********  Plot stt vs z SS                  ******************
      //*****************************************************************
      PT_can[j]->cd(1);
      gStyle->SetOptStat(10);
      gStyle->SetOptFit(0111);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      h2_ss->Draw("colz");
      h2_ss->SetTitle("Straight Section Corrected Time");
      h2_ss->SetXTitle("Z (cm)");
      h2_ss->SetYTitle("SC Time (ns)");
      h2_ss->GetZaxis()->SetRangeUser(0,4);
      h2_ss->GetXaxis()->SetRangeUser(0,40);
      h2_ss->GetYaxis()->SetRangeUser(-5.0,5.0);
      gPad->Update();
      //*****************************************************************
      //***********  Plot stt vs z BS                  ******************
      //*****************************************************************
      PT_can[j]->cd(2);
      gStyle->SetOptStat(10);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
            
      h2_bs->Draw("colz");
      h2_bs->SetTitle("Bend Section Corrected Time");
      h2_bs->SetXTitle("Z (cm)");
      h2_bs->SetYTitle("SC Time (ns)");
      h2_bs->GetZaxis()->SetRangeUser(0,4);
      h2_bs->GetXaxis()->SetRangeUser(39,43);
      h2_bs->GetYaxis()->SetRangeUser(-5.0,5.0);
      gPad->Update();
      //*****************************************************************
      //***********  Plot stt vs z NS                  ******************
      //*****************************************************************
      PT_can[j]->cd(3);
      gStyle->SetOptStat(10);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
      
      h2_ns->Draw("colz");
      h2_ns->SetTitle("Nose Section Corrected Time");
      h2_ns->SetXTitle("Z (cm)");
      h2_ns->SetYTitle("SC Time (ns)");
      h2_ns->GetZaxis()->SetRangeUser(0,4);
      h2_ns->GetXaxis()->SetRangeUser(43,60);
      h2_ns->GetYaxis()->SetRangeUser(-5.0,5.0);
      gPad->Update();
      //*****************************************************************
      //***********  Plot Total corrected time from each division *******
      //*****************************************************************
      PT_can[j]->cd(4);
      gStyle->SetOptStat(10);
      gStyle->SetErrorX(0); 
      gPad->SetTicks();
      gPad->SetGrid();
            
      h2_total->Draw("colz");
      h2_total->SetTitle("Paddle Corrected Time");
      h2_total->SetXTitle("Z (cm)");
      h2_total->SetYTitle("SC Time (ns)");
      h2_total->GetZaxis()->SetRangeUser(0,4);
      h2_total->GetXaxis()->SetRangeUser(0,60);
      h2_total->GetYaxis()->SetRangeUser(-5.0,5.0);
      gPad->Update();
    }
  // Close the output file
  st_prop_time_constants.close();
}

