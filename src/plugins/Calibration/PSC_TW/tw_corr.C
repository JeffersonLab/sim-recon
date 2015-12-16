// ROOT fitting macro for the PSC time-walk corrections
// Created 8/18/2015
// Alex Barnes
// University of Connecticut
//
// Notes:
// 	1) the input file must be changed if not using the current directory's
//	   hd_root.root file
//	2) Change #define GLOBAL_FIT to 0 or 1 if you want individual
//	   or global fits, respectively
//
// Example usage:
// 	> root -l
// 	root> .L tw_fit-rf.C
// 	root> tw_corr()
//
// Functions:
//    tw_plot(): This function takes the original root file branches and fills
//               the initial histograms
//    tw_proj(): This function takes the initial histograms and makes y projections
//               and fits these projections with a Gaussian. The means and sigmas
//               are used to make points on a TGraph
//    tw_fit():  This function takes the tw_proj() TGraph and fits it following Elton's
//               method
//    tw_corr(): This function uses the results of tw_fit() to perform a time-walk
//               correction on the data and fill corrected histograms
//    fitf():	 This function provides the timewalk fit function for tw_fit()

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TProfile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TString.h>

#include <fstream>

#define DEBUG false
#define OUTPUT false

// Declare files
TFile *infile;				// input ROOT file, 
TFile *outfile;				// output ROOT file

// Declare constants
const Int_t NMODULES = 8;		// Modules per arm
const Double_t THRESHOLD = 16; 		// discriminator threshold in ADC counts (40 mV)
const Int_t BINS_PP = 1000;		// number of bins for the pulse peak axis
const Int_t PP_MAX = 1000;		// maximum pulse peak for histograms
const Int_t PP_MIN = 0;			// minimum pulse peak for histograms
const Int_t BINS_DT = 100;		// number of bins for the time difference axis
const Int_t DT_MAX = 5;			// maximum time difference for histograms
const Int_t DT_MIN = -5;		// minimum time difference for histograms
const Int_t PROF_WIDTH = 82;		// number of bins for the TProfiles
const Int_t PROF_MIN = 180;		// minimum TProfile axis value 
const Int_t PROF_MAX = 1000;		// maximum TProfile axis value

// Declare variables
Double_t dtMeanL[NMODULES] = {};	// used to center the time dist. for the left modules
Double_t dtMeanR[NMODULES] = {};	// used to center the time dist. for the right modules
Double_t mpv[2][NMODULES] = {},{};	// used to find the most probable value for the pulse peaks

// Individual time walk fit parameters
Double_t c0[2][NMODULES] = {},{};	// constant paramter for each module
Double_t c1[2][NMODULES] = {},{};	// linear parameter for each module
Double_t c2[2][NMODULES] = {},{};	// exponent parameter for each module

// Declare 1D histograms
TH1D *h_dt_l[NMODULES];			// timing dist for each left module against the RF
TH1D *h_dt_r[NMODULES];			// timing dist for each right module against the RF
TH1D *h_dt_l_corr[NMODULES];		// timing dist for each left module against the RF
TH1D *h_dt_r_corr[NMODULES];		// timing dist for each right module against the RF
TH1D *h_pp_l[NMODULES];			// pulse peak distribution for the left modules
TH1D *h_pp_r[NMODULES];			// pulse peak distribution for the right modules

// Declare 2D histograms
TH2F *h_dt_vs_pp_l[NMODULES];		// dt vs pp for each left module
TH2F *h_dt_vs_pp_r[NMODULES];		// dt vs pp for each right module
TH2F *h_dt_vs_pp_l_corr[NMODULES];	// dt vs pp for each left module after corrections
TH2F *h_dt_vs_pp_r_corr[NMODULES];	// dt vs pp for each right module after corrections

//TGraphErrors *g_dt_vs_pp;		
TProfile *g_dt_vs_pp_l[NMODULES];	// dt vs pp plot converted into TProfile
TProfile *g_dt_vs_pp_r[NMODULES];	// dt vs pp plot converted into TProfile

void tw_plot(char* inputFile) {

   // Get input file, make output file
   infile = new TFile(inputFile);
   outfile = new TFile("results.root","RECREATE");

   // Create directories
#if DEBUG
   std::cout << "Creating output file directories" << std::endl;
#endif
   outfile->mkdir("dt_vs_pp");
   outfile->mkdir("dt_vs_pp/dt_vs_pp_L");
   outfile->mkdir("dt_vs_pp/dt_vs_pp_R");
   outfile->mkdir("dt_vs_pp_corr");
   outfile->mkdir("dt_vs_pp_corr/dt_vs_pp_corr_L");
   outfile->mkdir("dt_vs_pp_corr/dt_vs_pp_corr_R");
   outfile->mkdir("dt_vs_pp_profile");
   outfile->mkdir("dt_vs_pp_profile/dt_vs_pp_profile_L");
   outfile->mkdir("dt_vs_pp_profile/dt_vs_pp_profile_R");

   // Make canvas
   TCanvas *c = (TCanvas*)gROOT->FindObject("c");
   if (c == 0) {
      c = new TCanvas("c","c",0,20,600,500);
   }

   // Make histograms
#if DEBUG
   std::cout << "Naming histograms" << std::endl;
#endif
   // Per module histograms
   for (Int_t i = 0; i < NMODULES; ++i) {
      // left module
      h_dt_l[i] = new TH1D(Form("dt_l_%i",i+1),
                           Form("Time difference;PSC_L_%i - RF_L (ns)",i+1),
                           BINS_DT,DT_MIN,DT_MAX);
      h_pp_l[i] = new TH1D(Form("pp_l_%i",i+1),
                           Form("Pulse peak;PSC_L_%i pulse peak (adc counts)",i+1),
                           BINS_PP,PP_MIN,PP_MAX);
      h_dt_vs_pp_l[i] = (TH2F*)infile->Get(Form("h_dt_vs_pp_l_%i",i+1));
      h_dt_vs_pp_l_corr[i] = new TH2F(Form("h_dt_vs_pp_L_corr_%i",i+1),
                                      Form("Time difference vs pulse peak for PSC_L_%i",i+1),
                                      BINS_PP,PP_MIN,PP_MAX,BINS_DT,DT_MIN,DT_MAX);
      g_dt_vs_pp_l[i] = new TProfile(Form("g_dt_vs_pp_L_%i",i+1),
                                     Form("Time difference vs pulse peak for PSC_L_%i",i+1),
                                     PROF_WIDTH,PROF_MIN,PROF_MAX,"");
      g_dt_vs_pp_l[i]->SetMinimum(DT_MIN);
      g_dt_vs_pp_l[i]->SetMaximum(DT_MAX);

      // right module
      h_dt_r[i] = new TH1D(Form("dt_r_%i",i+1),
                           Form("Time difference;PSC_R_%i - RF_R (ns)",i+1),
                           BINS_DT,DT_MIN,DT_MAX);
      h_pp_r[i] = new TH1D(Form("pp_r_%i",i+1),
                           Form("Pulse peak;PSC_R_%i pulse peak (adc counts)",i+1),
                           BINS_PP,PP_MIN,PP_MAX);
      h_dt_vs_pp_r[i] = (TH2F*)infile->Get(Form("h_dt_vs_pp_r_%i",i+1));
      h_dt_vs_pp_r_corr[i] = new TH2F(Form("h_dt_vs_pp_R_corr_%i",i+1),
                                      Form("Time difference vs pulse peak for PSC_R_%i",i+1),
                                      BINS_PP,PP_MIN,PP_MAX,BINS_DT,DT_MIN,DT_MAX);
      g_dt_vs_pp_r[i] = new TProfile(Form("g_dt_vs_pp_R_%i",i+1),
                                     Form("Time difference vs pulse peak for PSC_R_%i",i+1),
                                     PROF_WIDTH,PROF_MIN,PROF_MAX,"");
      g_dt_vs_pp_r[i]->SetMinimum(DT_MIN);
      g_dt_vs_pp_r[i]->SetMaximum(DT_MAX);
   }

#if DEBUG
   std::cout << "Fitting initial timing distributions" << std::endl;
#endif
   for (Int_t i = 0; i < NMODULES; ++i) {
      // left module
      h_dt_l[i] = h_dt_vs_pp_l[i]->ProjectionY();
      TFitResultPtr resultpL0 = h_dt_l[i]->Fit("gaus","sq");
      double meanL = resultpL0->Parameter(1);
      double sigmaL = resultpL0->Parameter(2);

      // fit a second time using the mean and sigma of the initial fit
      TFitResultPtr resultpL = h_dt_l[i]->Fit("gaus","sq","",meanL-2*sigmaL,meanL+2*sigmaL);
      dtMeanL[i] = resultpL->Parameter(1);
#if OUTPUT 
      std::cout << "Sigma for left module " << i+1 << " is " << 1000*resultpL->Parameter(2) << std::endl;
#endif

      h_dt_l[i]->Reset();

      // right module
      h_dt_r[i] = h_dt_vs_pp_r[i]->ProjectionY();
      TFitResultPtr resultpR0 = h_dt_r[i]->Fit("gaus","sq");
      double meanR = resultpR0->Parameter(1);
      double sigmaR = resultpR0->Parameter(2);

      // fit a second time using the mean and sigma of the initial fit
      TFitResultPtr resultpR = h_dt_r[i]->Fit("gaus","sq","",meanR-2*sigmaR,meanR+2*sigmaR);
      dtMeanR[i] = resultpR->Parameter(1);
#if OUTPUT
      std::cout << "Sigma for right module " << i+1 << " is " << 1000*resultpR->Parameter(2) << std::endl;
#endif

      h_dt_r[i]->Reset();
   }

#if DEBUG
   std::cout << "Making TProfile" << std::endl;
#endif
   for (Int_t m = 0; m < NMODULES; ++m) {
      for (Int_t i = 1; i <= 1000; ++i) {
         for (Int_t j = 1; j <= 100; ++j) {
            // left module
            Double_t content_l = h_dt_vs_pp_l[m]->GetBinContent(i,j);
            double t_shift_l = -5 + 0.1*j - dtMeanL[m];
            if (content_l > 0 && t_shift_l < -2)
               g_dt_vs_pp_l[m]->Fill(i,t_shift_l+4,content_l);
            else if (content_l > 0 && t_shift_l > 2)
               g_dt_vs_pp_l[m]->Fill(i,t_shift_l-4,content_l);
            else
               g_dt_vs_pp_l[m]->Fill(i,t_shift_l,content_l);

            // right module
            Double_t content_r = h_dt_vs_pp_r[m]->GetBinContent(i,j);
            double t_shift_r = -5 + 0.1*j - dtMeanR[m];
            if (content_r > 0 && t_shift_r < -2)
               g_dt_vs_pp_r[m]->Fill(i,t_shift_r+4,content_r);
            else if (content_r > 0 && t_shift_r > 2)
               g_dt_vs_pp_r[m]->Fill(i,t_shift_r-4,content_r);
            else
               g_dt_vs_pp_r[m]->Fill(i,t_shift_r,content_r);
         }
      }
   }

#if DEBUG
   std::cout << "Getting most probable value" << std::endl;
#endif
   // Get the most probable value of the pulse peak distribution for each module
   for (Int_t i = 0; i < NMODULES; ++i) {
      // left module
      h_pp_l[i] = h_dt_vs_pp_l[i]->ProjectionX();
      TFitResultPtr ptrL1 = h_pp_l[i]->Fit("landau","sq");
      double mpv_L = ptrL1->Parameter(1);
      double sigma_pp_L = ptrL1->Parameter(2);
      TFitResultPtr ptrL2 = h_pp_l[i]->Fit("landau","sq","",mpv_L-1.5*sigma_pp_L,mpv_L+1.5*sigma_pp_L);
      mpv[0][i] = ptrL2->Parameter(1);

      // right module
      h_pp_r[i] = h_dt_vs_pp_r[i]->ProjectionX();
      TFitResultPtr ptrR1 = h_pp_r[i]->Fit("landau","sq");
      double mpv_R = ptrR1->Parameter(1);
      double sigma_pp_R = ptrR1->Parameter(2);
      TFitResultPtr ptrR2 = h_pp_r[i]->Fit("landau","sq","",mpv_R-1.5*sigma_pp_R,mpv_R+1.5*sigma_pp_R);
      mpv[1][i] = ptrR2->Parameter(1);
   }

#if DEBUG
   std::cout << "Writing to file ..." << std::endl;
#endif
   for (Int_t i = 0; i < NMODULES; ++i) {
      // left module
      outfile->cd("dt_vs_pp/dt_vs_pp_L");
      h_dt_vs_pp_l[i]->Write();

      outfile->cd("dt_vs_pp_profile/dt_vs_pp_profile_L");
      g_dt_vs_pp_l[i]->Write();

      // right module
      outfile->cd("dt_vs_pp/dt_vs_pp_R");
      h_dt_vs_pp_r[i]->Write();

      outfile->cd("dt_vs_pp_profile/dt_vs_pp_profile_L");
      g_dt_vs_pp_l[i]->Write();
   }

   for (int j = 0; j < NMODULES; ++j) {
#if DEBUG
      std::cout << "Calling tw_fit() for modules " << j << std::endl;
#endif
      // left module
      tw_fit(g_dt_vs_pp_l[j],0,j);
      
      // right module
      tw_fit(g_dt_vs_pp_r[j],1,j);
   }
}

// Individual tw_fit()
void tw_fit(TProfile *h_tw, int arm, int module) {

   // fitf defined at the bottom of the file
   TF1 *ftw = new TF1("ftw",fitf,100.,1000.,3);
#if DEBUG
   std::cout << "Defined time-walk fit function" << std::endl;
#endif

   ftw->SetParameter(0,-1.9);
   ftw->SetParameter(1,4.8);
   ftw->SetParameter(2,-0.4);
   ftw->SetParName(0,"c0");
   ftw->SetParName(1,"c1");
   ftw->SetParName(2,"c2");

   h_tw->Fit("ftw","RWq");
#if DEBUG
   std::cout << "Fit histogram" << std::endl;
#endif

   c0[arm][module] = ftw->GetParameter("c0");
   c1[arm][module] = ftw->GetParameter("c1");
   c2[arm][module] = ftw->GetParameter("c2");
#if OUTPUT
   std::cout << "c0: " << c0[arm][module] << std::endl;
   std::cout << "c1: " << c1[arm][module] << std::endl;
   std::cout << "c2: " << c2[arm][module] << std::endl;
#endif
}

void tw_corr(char* inputFile) {

   std::cout << "Input file: " << inputFile << std::endl;

   tw_plot(inputFile);

#if DEBUG
   std::cout << "Writing psc_tw_parms.out" << std::endl;
#endif
   // Writing parameters to psc_tw_parms.out for use in CCDB
   std::ofstream flog("psc_tw_parms.out");
   for (Int_t i = 0; i < 2; ++i) {
      for (Int_t j = 0; j < NMODULES; ++j) {
         flog << i         << std::setw(15)
              << j+1       << std::setw(15)
              << c0[i][j]  << std::setw(15)
              << c1[i][j]  << std::setw(15)
              << c2[i][j]  << std::setw(15)
              << THRESHOLD << std::setw(15)
              << mpv[i][j] << std::endl;
      }
   }

#if DEBUG
   std::cout << "About to fill corrected histograms" << std::endl;
#endif
   for (Int_t i = 0; i < NMODULES; ++i) {
      for (Int_t j = 0; j < 1000; ++j) {
         for (Int_t k = 0; k < 100; ++k) {
            // left module
            double content_l = h_dt_vs_pp_l[i]->GetBinContent(j,k);
            double t_shift_l = -5 + 0.1*k - dtMeanL[i];
            if (content_l > 0 && t_shift_l > 2)
               t_shift_l -= 4;
            else if (content_l > 0 && t_shift_l < -2)
               t_shift_r += 4;
            t_shift_l -= c1[0][i]*pow((j/THRESHOLD),c2[0][i]) - c1[0][i]*pow(mpv[0][i]/THRESHOLD,c2[0][i]);
            h_dt_vs_pp_l_corr[i]->Fill(j,t_shift_l,content_l);

            // right module
            double content_r = h_dt_vs_pp_r[i]->GetBinContent(j,k);
            double t_shift_r = -5 + 0.1*k - dtMeanR[i];
            if (content_r > 0 && t_shift_r > 2)
               t_shift_r -= 4;
            else if (content_r > 0 && t_shift_r < -2)
               t_shift_r += 4;
            t_shift_r -= c1[1][i]*pow((j/THRESHOLD),c2[1][i]) - c1[1][i]*pow(mpv[1][i]/THRESHOLD,c2[1][i]);
            h_dt_vs_pp_r_corr[i]->Fill(j,t_shift_r,content_r);
         }
      }
   }

   for (Int_t i = 0; i < NMODULES; ++i) {
      // left module
      outfile->cd("dt_vs_pp_corr/dt_vs_pp_corr_L");
      h_dt_vs_pp_l_corr[i]->Write();

      // right module
      outfile->cd("dt_vs_pp_corr/dt_vs_pp_corr_R");
      h_dt_vs_pp_r_corr[i]->Write();
   }

   // calculate new sigmas
   double avg = 0;
   double avgL = 0;
   double avgR = 0;

   std::ofstream flog("sigmas.out");
   flog << "Arm\t" << "Module\t" << "Sigma (ps)" << std::endl;

   for (Int_t i = 0; i < NMODULES; ++i) {
      // left module
      h_dt_l_corr[i] = h_dt_vs_pp_l_corr[i]->ProjectionY();
      TFitResultPtr ptrL = h_dt_l_corr[i]->Fit("gaus","sq");
#if OUTPUT
      std::cout << "Sigma for left module " << i+1 << ": " << (ptrL->Parameter(2))*1000 << std::endl;
#endif
      avgL += (ptrL->Parameter(2))*1000;
      flog << "0\t" << i+1 << "\t" << (ptrL->Parameter(2))*1000 << std::endl;

      // right module
      h_dt_r_corr[i] = h_dt_vs_pp_r_corr[i]->ProjectionY();
      TFitResultPtr ptrR = h_dt_r_corr[i]->Fit("gaus","sq");
#if OUTPUT
      std::cout << "Sigma for right module " << i+1 << ": " << (ptrR->Parameter(2))*1000 << std::endl;
#endif
      avgR += (ptrR->Parameter(2))*1000;
      flog << "1\t" << i+1 << "\t" << (ptrR->Parameter(2))*1000 << std::endl;

      avg += (ptrL->Parameter(2) + ptrR->Parameter(2))*1000;
   }
   avgL /= 8;
   avgR /= 8;
   avg  /= 16;

#if OUTPUT
   std::cout << "avgL: " << avgL << std::endl;
   std::cout << "avgR: " << avgR << std::endl;
   std::cout << "avg: "  << avg  << std::endl;
#endif

   delete c;
}

// Define the fit function, fitf
Double_t fitf(Double_t *x, Double_t *par) {
   Float_t xx = x[0];
   Double_t f = par[0] + par[1]*pow((xx/THRESHOLD),par[2]);
   return f;
}
