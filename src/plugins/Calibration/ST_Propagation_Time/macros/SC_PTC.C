// File: st_tw_fits.C
// Last Modified: 01/12/2016
// Creator: Mahmoud Kamel mkame006@fiu.edu
// Purpose: Displaying histograms and propagation time correction.
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
Double_t t_ns_fit_mean[NCHANNELS][300][3];
Double_t t_ns_fit_mean_err[NCHANNELS][300][3];
Double_t t_ss_fit_mean[NCHANNELS][300][3];
Double_t t_ss_fit_mean_err[NCHANNELS][300][3];
Double_t t_bs_fit_mean[NCHANNELS][300][3];
Double_t t_bs_fit_mean_err[NCHANNELS][300][3];
Double_t t_ss_fit_sigma[NCHANNELS][300][3];
Double_t t_ss_fit_sigma_err[NCHANNELS][300][3];
Double_t t_bs_fit_sigma[NCHANNELS][300][3];
Double_t t_bs_fit_sigma_err[NCHANNELS][300][3];
Double_t t_ns_fit_sigma[NCHANNELS][300][3];
Double_t t_ns_fit_sigma_err[NCHANNELS][300][3];
// Declare canvas
TCanvas *PT_can[30];
TCanvas *PTC[30];
TCanvas *PTC_BS[30];
TCanvas *PTC_NS[30];
TDirectory* TopDirectory;
TH2I* h2_ss;
TH2I* h2_bs;
TH2I* h2_ns;
TH2I* h2All;
double bin1[30][100];
double bin2[30][100];
TH1D* h1_ns_Bin[30][100];
TH1D* h1_ss_Bin[30][100];
TH1D* h1_bs_Bin[30][100];
Double_t t_max[NCHANNELS][300];
Double_t t_max_BS[NCHANNELS][300];
Double_t t_max_NS[NCHANNELS][300];
double bin1bs[30][100];
double bin2bs[30][100]; 
double bin1ns[30][100];
double bin2ns[30][100]; 


Double_t t_z_ss_fit_slope[NCHANNELS][3];
Double_t t_z_ss_fit_slope_err[NCHANNELS][3];
Double_t t_z_ss_fit_intercept[NCHANNELS][3];
Double_t t_z_ss_fit_intercept_err[NCHANNELS][3];
Double_t t_z_bs_fit_slope[NCHANNELS][3];
Double_t t_z_bs_fit_slope_err[NCHANNELS][3];
Double_t t_z_bs_fit_intercept[NCHANNELS][3];
Double_t t_z_bs_fit_intercept_err[NCHANNELS][3];
Double_t t_z_ns_fit_slope[NCHANNELS][3];
Double_t t_z_ns_fit_slope_err[NCHANNELS][3];
Double_t t_z_ns_fit_intercept[NCHANNELS][3];
Double_t t_z_ns_fit_intercept_err[NCHANNELS][3];
const int NOXbins = 300;
const int NOYbins = 200;
int sss[NCHANNELS][NOXbins][NOYbins];
Double_t binsss[NCHANNELS][NOXbins][NOYbins];
Double_t ss_max[NCHANNELS];
Double_t low_cut = 0.10;
Double_t global[NCHANNELS];

Double_t N_x_bin[NCHANNELS];
Double_t N_y_bin[NCHANNELS];
Double_t max_y_axis_bin[NCHANNELS];
Double_t max_y_axis_value[NCHANNELS];
Double_t bins[500][500];

//h1_ns_Bin = new TH1D *[300];
//TH1D** h1_ns_Bin = new TH1D *[300];
void SC_PTC(char const* input_filename)
//void st_tw_fits()
{
  TFile *df = new TFile(input_filename);
  std::ofstream SC_PTC;
  SC_PTC.open ("SC_PTC.txt", std::ofstream::out);
  //*****************************************************************
  //*********** Grab the histograms from the root file **************
  //*****************************************************************
    for (unsigned int j = 0; j < NCHANNELS; j++)
      //  for (unsigned int j = 0; j < 1; j++)
    {
      PTC[j] = new TCanvas( Form("PTC_%i",j+1),  Form("PTC_%i",j+1), 800, 450);
      PTC[j]->Divide(3,2);
      PTC_BS[j] = new TCanvas( Form("PTC_BS_%i",j+1),  Form("PTC_BS_%i",j+1), 800, 450);
      PTC_BS[j]->Divide(3,1);
      PTC_NS[j] = new TCanvas( Form("PTC_NS_%i",j+1),  Form("PTC_NS_%i",j+1), 800, 450);
      PTC_NS[j]->Divide(3,3);
      // The top directory 
      TopDirectory = (TDirectory*) df->FindObjectAny("ST_Propagation_Time"); 
      TopDirectory->cd();
      // Grab the uncorrected 2D histogram 
      char* All = Form("h2_PpropTime_z_%i",j+1);
      h2All = (TH2I*) TopDirectory->FindObjectAny(All);
      
      cout << "==================================================" << endl;
      cout << "Processing Channel " << j+1 << endl;
      cout << "==================================================" << endl;
      //*****************************************************************
      //***********  SS intervals                      ******************
      //*****************************************************************
      
      for (Int_t i = 10; i < 40; i+=5)
	{ 
	  Int_t k =(i/5)-1;
	  cout << " k = " << k << endl;
	  PTC[j]->cd(k);
	  bin1[j][i] = h2All->GetXaxis()->FindBin(i);
	  bin2[j][i] = h2All->GetXaxis()->FindBin(i+5);
	  h1_ss_Bin[j][i] = h2All->ProjectionY(Form("h1_ss_Bin_%i_%i",i,j+1),bin1[j][i],bin2[j][i],"");

	  t_max[j][i]=  h1_ss_Bin[j][i]->GetMaximumBin();
	  double LL_fit = (t_max[j][i] -100)*0.1 - 0.5 ;
	  double UL_fit = (t_max[j][i] -100)*0.1 + 0.5 ;
	  h1_ss_Bin[j][i]->Draw();
	  h1_ss_Bin[j][i]->Fit("gaus","","",LL_fit,UL_fit);
	  TF1 *gaus = h1_ss_Bin[j][i]->GetFunction("gaus");
	  t_ss_fit_mean[j][i][1] = gaus->GetParameter(1);
	  t_ss_fit_mean_err[j][i][1] = gaus->GetParError(1);
	  t_ss_fit_sigma[j][i][2] = gaus->GetParameter(2);
	  t_ss_fit_sigma_err[j][i][2] = gaus->GetParError(2);
	  
	  SC_PTC <<"\t"<< j+1 <<"\t"<< i+2.5 <<"\t"<< 0.0 <<"\t"<< t_ss_fit_mean[j][i][1] <<"\t"<< t_ss_fit_mean_err[j][i][1]<<"\t"<< t_ss_fit_sigma[j][i][2] << "\t"<< t_ss_fit_sigma_err[j][i][2] <<endl;
	}
      //*****************************************************************
      //***********  BS intervals                      ******************
      //*****************************************************************
      for (Int_t i = 40; i < 43; i+=1)
      	{ 
	  
	  Int_t k = i - 39;
      	  PTC_BS[j]->cd(k);
      	  bin1bs[j][i] = h2All->GetXaxis()->FindBin(i);
      	  bin2bs[j][i] = h2All->GetXaxis()->FindBin(i+1);
      	  h1_bs_Bin[j][i] = h2All->ProjectionY(Form("h1_bs_Bin_%i_%i",i,j+1),bin1bs[j][i],bin2bs[j][i],"");

      	  t_max_BS[j][i]=  h1_bs_Bin[j][i]->GetMaximumBin();
      	  double LL_fit = (t_max_BS[j][i] -100)*0.1 - 0.5 ;
      	  double UL_fit = (t_max_BS[j][i] -100)*0.1 + 0.5 ;
      	  h1_bs_Bin[j][i]->Draw();
      	  h1_bs_Bin[j][i]->Fit("gaus","","",LL_fit,UL_fit);
      	  TF1 *gaus_BS = h1_bs_Bin[j][i]->GetFunction("gaus");
      	  t_bs_fit_mean[j][i][1] = gaus_BS->GetParameter(1);
      	  t_bs_fit_mean_err[j][i][1] = gaus_BS->GetParError(1);
      	  t_bs_fit_sigma[j][i][2] = gaus_BS->GetParameter(2);
      	  t_bs_fit_sigma_err[j][i][2] = gaus_BS->GetParError(2);
	  
      	  SC_PTC <<"\t"<< j+1 <<"\t"<< i+0.5 <<"\t"<< 0.0 <<"\t"<< t_bs_fit_mean[j][i][1] <<"\t"<< t_bs_fit_mean_err[j][i][1]<<"\t"<< t_bs_fit_sigma[j][i][2] << "\t"<< t_bs_fit_sigma_err[j][i][2] <<endl;
      	}
      
      // //*****************************************************************
      // //***********  NS intervals                      ******************
      // //*****************************************************************
      Int_t k =0;
      for (Int_t i = 43; i < 61; i+=2)
      	{ 
	  k = (i - 42) - k;
      	  PTC_NS[j]->cd(k);
      	  bin1ns[j][i] = h2All->GetXaxis()->FindBin(i);
      	  bin2ns[j][i] = h2All->GetXaxis()->FindBin(i+2);
      	  h1_ns_Bin[j][i] = h2All->ProjectionY(Form("h1_ns_Bin_%i_%i",i,j+1),bin1ns[j][i],bin2ns[j][i],"");
	  
      	  t_max_NS[j][i]=  h1_ns_Bin[j][i]->GetMaximumBin();
      	  double LL_fit = (t_max_NS[j][i] -100)*0.1 - 0.5 ;
      	  double UL_fit = (t_max_NS[j][i] -100)*0.1 + 0.5 ;
      	  h1_ns_Bin[j][i]->Draw();
      	  h1_ns_Bin[j][i]->Fit("gaus","","",LL_fit,UL_fit);
      	  TF1 *gaus_NS = h1_ns_Bin[j][i]->GetFunction("gaus");
      	  t_ns_fit_mean[j][i][1] = gaus_NS->GetParameter(1);
      	  t_ns_fit_mean_err[j][i][1] = gaus_NS->GetParError(1);
      	  t_ns_fit_sigma[j][i][2] = gaus_NS->GetParameter(2);
      	  t_ns_fit_sigma_err[j][i][2] = gaus_NS->GetParError(2);
	  
      	  SC_PTC <<"\t"<< j+1 <<"\t"<< i+1 <<"\t"<< 0.0 <<"\t"<< t_ns_fit_mean[j][i][1] <<"\t"<< t_ns_fit_mean_err[j][i][1]<<"\t"<< t_ns_fit_sigma[j][i][2] << "\t"<< t_ns_fit_sigma_err[j][i][2] <<endl;
      	}
    }
  
  SC_PTC.close();
}

