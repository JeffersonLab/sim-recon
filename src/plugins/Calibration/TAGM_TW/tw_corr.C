// ROOT fitting macro for the TAGM time-walk corrections
// Created 8/18/2015
// Alex Barnes
// University of Connecticut
//
// Example usage:
// 	> root -l -b -q 'tw_corr.C("hd_root.root")'
//
// Functions:
//    tw_plot(): This function takes the original root file branches and fills
//               the initial histograms
//    tw_fit():  This function takes the tw_proj() TGraph and fits it following Elton's
//               method
//    tw_corr(): This function uses the results of tw_fit() to pepscorm a time-walk
//               correction on the data and fill corrected histograms

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <sstream>
#include <iostream>

#define Ecut false
#define OUTPUT false
#define DEBUG false

// Declare files
TFile *infile;
TFile *outfile;
std::ofstream flog("tagm_tw_parms.out");

// Declare constants
const Int_t NCOLUMNS = 100;	// number of groups for the calibration
const Int_t NROWS = 5;		// number of rows in the TAGM
const Int_t BINS_PP = 1000;	// number of bins for the pulse peak axis
const Int_t PP_MAX = 1000;	// maximum pulse peak for histograms
const Int_t PP_MIN = 0;		// minimum pulse peak for histograms
const Int_t BINS_DT = 100;	// number of bins for the time difference axis
const Int_t DT_MAX = 5;		// maximum time difference for the adjusted histograms
const Int_t DT_MIN = -5;	// minimum time difference for the adjusted histograms
const Int_t BINS_PROF = 46;	// number of bins for the TProfile
const Int_t PROF_MIN = 80;	// minimum TProfile axis value
const Int_t PROF_MAX = 1000;	// maximum TProfile axis value
const Int_t PROF_MAX = PROF_MIN+BINS_PROF*20;	// maximum TProfile axis value
const Double_t THRESHOLD = 8;	// discriminator threshold in ADC counts (20 mV)

double dtMean[NCOLUMNS] = {};	// used to center the time dist.
double sigmas[2][100] = {}, {};	// 
double mpv[NCOLUMNS] = {};
double mpv_sig[NCOLUMNS] = {};

// Time-walk fit parameters
double c0[NCOLUMNS] = {};
double c1[NCOLUMNS] = {};
double c2[NCOLUMNS] = {};

// Declare 1D histograms
TH1D *h_dt[NCOLUMNS];		// tm - rf
TH1D *h_dt_corr[NCOLUMNS];	// tm - rf centered about 0 with t.w. corrections
TH1D *h_pp[NCOLUMNS];		// pulse height spectrum

// Declare 2D histograms
TH2F *h_dt_vs_pp[NCOLUMNS];	// tm - rf vs pulse peak
TH2F *h_dt_vs_pp_0[NCOLUMNS];	// tm - rf vs pulse peak centered about delta t = 0
TH2F *h_dt_vs_pp_corr[NCOLUMNS];// tm - rf vs pulse peak, t.w. corrected
TProfile *g_dt_vs_pp[NCOLUMNS];	// tm - rf vs pulse peak projection plot

void tw_plot(char* inputFile) {

   // Get input file, make output file
   infile = new TFile(inputFile);
   outfile = new TFile("results.root","RECREATE");

   // Create directories
#if DEBUG
   std::cout << "Creating output file directories" << std::endl;
#endif
   outfile->mkdir("dt_vs_pp");
   outfile->mkdir("dt_vs_pp_0");
   outfile->mkdir("dt_vs_pp_corr");
   outfile->mkdir("dt_vs_pp_profile");

   // Make canvas
   TCanvas *c = (TCanvas*)gROOT->FindObject("c");
   if (c == 0) {
      c = new TCanvas("c","c",0,20,600,500);
   }

#if DEBUG
   std::cout << "Naming histograms" << std::endl;
#endif
   for (Int_t i = 0; i < NCOLUMNS; ++i) {
      // Make histograms
      h_dt[i] = new TH1D(Form("dt_%i",i+1),
                         Form("Time difference;Col_%i - RF (ns)",i+1),
                         2*BINS_DT,-37,-17);

      h_dt_corr[i] = new TH1D(Form("dt_corr_%i",i+1),
                              Form("Time difference;Col_%i - RF (ns)",i+1),
                              2*BINS_DT,DT_MIN,DT_MAX);

      h_pp[i] = new TH1D(Form("pp_%i",i+1),
                         Form("Pulse Peak;Col_%i Pulse Peak (adc counts)",i+1),
                         BINS_PP,PP_MIN,PP_MAX);

      h_dt_vs_pp[i] = (TH2F*)infile->Get(Form("h_dt_vs_pp_%i",i+1));

      h_dt_vs_pp_0[i] = new TH2F(Form("h_dt_vs_pp_0_%i",i+1),
                                 Form("Time difference vs pulse peak for column %i;\
                                 Pulse peak (adc counts);TAGM - RF (ns)",i+1),
                                 BINS_PP,PP_MIN,PP_MAX,BINS_DT,DT_MIN,DT_MAX);

      h_dt_vs_pp_corr[i] = new TH2F(Form("h_dt_vs_pp_corr_%i",i+1),
                                    Form("Time difference vs pulse peak for column %i;\
                                    Pulse peak (adc counts);TAGM - RF (ns)",i+1),
                                    BINS_PP,PP_MIN,PP_MAX,BINS_DT,DT_MIN,DT_MAX);

      g_dt_vs_pp[i] = new TProfile(Form("g_dt_vs_pp_%i",i+1),
                                   Form("Time difference vs pulse peak for col_%i;\
                                   Pulse peak (adc counts);TAGM - RF (ns)",i+1),
                                   BINS_PROF,PROF_MIN,PROF_MAX, "s");
      g_dt_vs_pp[i]->SetMinimum(DT_MIN);
      g_dt_vs_pp[i]->SetMaximum(DT_MAX);
   }

#if DEBUG
   std::cout << "Fitting initial timing distributions" << std::endl;
#endif
   for (Int_t i = 0; i < NCOLUMNS; ++i) {
      h_dt[i] = h_dt_vs_pp[i]->ProjectionY();
      Int_t binMax = h_dt[i]->GetMaximumBin();
      TFitResultPtr resultp0 = h_dt[i]->Fit("gaus","sq","",-37+binMax*0.1-2,-37+binMax*0.1+2);
      if (resultp0->IsEmpty()) {
         dtMean[i] = 0;
#if DEBUG
         std::cout << "Missing channel: " << i+1 << std::endl;
#endif
         continue;
      } 
      double mean = resultp0->Parameter(1);
      double sigma = resultp0->Parameter(2);
      TFitResultPtr resultp = h_dt[i]->Fit("gaus","sq","",mean-2*sigma,mean+2*sigma);
      if (resultp->IsEmpty()) {
         dtMean[i] = 0;
#if DEBUG
         std::cout << "Something went wrong with the second fit for column: " << i+1 << std::endl;
#endif
         continue;
      }
      dtMean[i] = resultp->Parameter(1);

#if OUTPUT
      std::cout << "Time shift for column " << i+1 << ": "
                << dtMean[i] << std::endl;

      std::cout << "Initial sigma for column " << i+1 << ": "
                << (resultp->Parameter(2))*1000 << std::endl;
#endif

      sigmas[0][i] = (resultp->Parameter(2))*1000;

      h_dt[i]->Reset();
   }

   // Fill histograms with adjusted timing which centers around 0
#if DEBUG
   std::cout << "Making TProfile" << std::endl;
#endif
   for (Int_t col = 0; col < NCOLUMNS; ++col) {
      for (Int_t i = 1; i <= 1000; ++i) {
         for (Int_t j = 1; j <= 200; ++j) {
            Double_t content = h_dt_vs_pp[col]->GetBinContent(i,j);
            Double_t t_shift = -37 + 0.1*(j-1) - dtMean[col];
            if (content > 0 && (t_shift < -2.5 || t_shift > 2.5))
               content = 0;
            else if (content > 0 && t_shift < -2)
               t_shift+4;
            else if (content > 0 && t_shift > 2)
               t_shift-4;
            g_dt_vs_pp[col]->Fill(i,t_shift,content);
            h_dt_vs_pp_0[col]->Fill(i,t_shift,content);
         }
      }
      for (Int_t i = 1; i <= BINS_PROF; ++i) {
         Double_t binEntries = g_dt_vs_pp[col]->GetBinEntries(i);
         if (binEntries < 4) {
            g_dt_vs_pp[col]->SetBinEntries(i,0);
         }
      }
   }

   // Get the most probable value of the pulse peak distribution
#if DEBUG
   std::cout << "Getting most probable value" << std::endl;
#endif
   for (Int_t i = 0; i < NCOLUMNS; ++i) {
      h_pp[i] = h_dt_vs_pp[i]->ProjectionX();
      double integral = h_pp[i]->Integral(0,1000);
      if (integral < 100) {
#if OUTPUT
         std::cout << "Channel missing or low statistics: " << i+1 << std::endl;
#endif
         continue;
      }
      h_pp[i]->GetXaxis()->SetRange(150,1000);
      double mean_i = h_pp[i]->GetMean();
      double rms_i = h_pp[i]->GetRMS();
      integral = h_pp[i]->Integral(150,1000);
      if (integral > 100) {
         TFitResultPtr ptr = h_pp[i]->Fit("gaus","sq","",mean_i-1.5*rms_i,mean_i+1.5*rms_i);
         if (ptr->Parameter(1) < 200) {
            h_pp[i]->GetXaxis()->SetRange(0,400);
            mean_i = h_pp[i]->GetMean();
            rms_i = h_pp[i]->GetRMS();
            ptr = h_pp[i]->Fit("gaus","sq");
         }
         mpv[i] = ptr->Parameter(1);
         mpv_sig[i] = ptr->Parameter(2);
      }
      else {
         h_pp[i]->GetXaxis()->SetRange(0,400);
         mean_i = h_pp[i]->GetMean();
         rms_i = h_pp[i]->GetRMS();
         if (mean_i > 0) {
            TFitResultPtr ptr = h_pp[i]->Fit("gaus","sq");
            mpv[i] = ptr->Parameter(1);
            mpv_sig[i] = ptr->Parameter(2);
         }
         else {
            mpv[i] = THRESHOLD;
            mpv_sig[i] = 0;
         }
      }
      h_pp[i]->GetXaxis()->SetRange(0,1000);
   }

#if DEBUG
   std::cout << "Writing to file ..." << std::endl;
#endif
   for (Int_t i = 0; i < NCOLUMNS; ++i) {
      outfile->cd("dt_vs_pp");
      h_dt_vs_pp[i]->Write();

      outfile->cd("dt_vs_pp_0");
      h_dt_vs_pp_0[i]->Write();
   }

#if DEBUG
   std::cout << "Calling tw_fit()" << std::endl;
#endif

   for (Int_t i = 0; i < NCOLUMNS; ++i) {
      tw_fit(g_dt_vs_pp[i], i);
   }
}

void tw_fit(TProfile *h_tw, int col) {
#if DEBUG
   std::cout << "Making fit function for column " << col+1 << std::endl;
#endif
   // 50 is the discriminator threshold in units of adc counts. 50 is 15 mV
   // The RF had a threshold of 40 counts for the 3185 run.
   // The TM had a threshold of 20 counts for the 3185 run.
   TF1 *ftw = new TF1("ftw",fitf,50.,1000.,3); // fitf defined at the bottom
#if DEBUG
   std::cout << "Defined time-walk fit function for column " << col+1 << std::endl;
#endif

   ftw->SetParameter(0,-0.13);
   ftw->SetParameter(1,28.7);
   ftw->SetParameter(2,-1.93);
   ftw->SetParName(0,"c0");
   ftw->SetParName(1,"c1");
   ftw->SetParName(2,"c2");

   TFitResultPtr ptr = h_tw->Fit("ftw","sRWq");
#if DEBUG
   std::cout << "Fit histogram" << std::endl;
#endif

   c0[col] = ftw->GetParameter("c0");
   c1[col] = ftw->GetParameter("c1");
   c2[col] = ftw->GetParameter("c2");

   if (ptr->IsEmpty()) {
      c0[col] = 1;
      c1[col] = -1;
      c2[col] = 0;
   }
#if OUTPUT
   std::cout << "c0: " << c0[col] << std::endl;
   std::cout << "c1: " << c1[col] << std::endl;
   std::cout << "c2: " << c2[col] << std::endl;
#endif

   flog << "0"       << std::setw(15) 
        << col+1     << std::setw(15) 
        << c0[col]   << std::setw(15)
        << c1[col]   << std::setw(15)
        << c2[col]   << std::setw(15)
        << THRESHOLD << std::setw(15)
        << mpv[col]  << std::endl;
   if (col+1 == 9 || col+1 == 27 || col+1 == 81 || col+1 == 99) {
      for (int row = 1; row <= NROWS; ++row) {
         flog << row		<< std::setw(15)
              << col+1		<< std::setw(15)
              << c0[col]	<< std::setw(15)
              << c1[col]	<< std::setw(15)
              << c2[col]	<< std::setw(15)
              << THRESHOLD	<< std::setw(15)
              << mpv[col] 	<< std::endl;
      }
   }
   outfile->cd("dt_vs_pp_profile");
   h_tw->Write();
}

void tw_corr(char* inputFile) {
   // Use the results of tw_fit() to make time-walk corrected histograms
   std::cout << "Input file: " << inputFile << std::endl;

   tw_plot(inputFile);

   // Add default parameters for unused columns
   for (int i = 1; i < 3; ++i) {
      flog << "0"	<< std::setw(15)
           << NCOLUMNS+i<< std::setw(15)
           << "1"	<< std::setw(15)
           << "-1"	<< std::setw(15)
           << "0"	<< std::setw(15)
           << THRESHOLD << std::setw(15)
           << THRESHOLD	<< std::endl;
   }

#if DEBUG
   std::cout << "About to fill corrected histograms" << std::endl;
#endif
   for (Int_t col = 0; col < NCOLUMNS; ++col) {
      for (Int_t i = 1; i <= 1000; ++i) {
         for (Int_t j = 1; j <= 200; ++j) {
            Double_t content = h_dt_vs_pp[col]->GetBinContent(i,j);
            Double_t t_shift = -37 + 0.1*(j-1) - dtMean[col];
            if (content > 0 && (t_shift < -2.5 || t_shift > 2.5))
               content = 0;
            else if (content > 0 && t_shift < -2)
               t_shift += 4;
            else if (content > 0 && t_shift > 2)
               t_shift -= 4;
            // Use time-walk parameters
            if (mpv[col] > 150)
               t_shift -= c1[col]*pow((i/THRESHOLD),c2[col]) - c1[col]*pow((mpv[col]/THRESHOLD),c2[col]);
            else
               t_shift -= c1[col]*pow((i/THRESHOLD),c2[col]) - c1[col]*pow((110/THRESHOLD),c2[col]);
            h_dt_vs_pp_corr[col]->Fill(i,t_shift,content);
         }
      }
   }

   for (Int_t i = 0; i < NCOLUMNS; ++i) {
      outfile->cd("dt_vs_pp_corr");
      h_dt_vs_pp_corr[i]->Write();
   }

   double avg = 0;
   int tot = 0;
   std::ofstream flog("sigmas.out");
   flog << "Column" << "\t"
        << "Initial (ps)" << "\t"
        << "Final (ps)" << std::endl;
   for (Int_t i = 0; i < NCOLUMNS; ++i) {
#if DEBUG
      std::cout << "Calculating sigma for column " << i+1 << std::endl;
#endif
      h_dt_corr[i] = h_dt_vs_pp_corr[i]->ProjectionY();
      Int_t binMax = h_dt_corr[i]->GetMaximumBin();
      TFitResultPtr ptr = h_dt_corr[i]->Fit("gaus","sq","",-5+binMax*0.1-2,-5+binMax*0.1+2);
      if ((int)ptr < 0) continue;
      double mean = ptr->Parameter(1);
      double sigma = ptr->Parameter(2);
      if (mean < -5) {
         avg+= 0;
         sigmas[1][i] = 0;
      }
      else {
         TFitResultPtr ptr2 = h_dt_corr[i]->Fit("gaus","sq","",mean-2*sigma,mean+2*sigma);
         avg += (ptr2->Parameter(2))*1000;
         sigmas[1][i] = (ptr2->Parameter(2))*1000;
         tot++;
      }

      flog << i+1 << "\t" << sigmas[0][i]
           << "\t\t" << sigmas[1][i]
           << std::endl;
   }
   avg  /= tot;
   flog << "avg\t\t\t" << avg << std::endl;

   delete c;
}

Double_t fitf(Double_t *x, Double_t *par) {
   Float_t xx = x[0];
   Double_t f = par[0] + par[1]*pow((xx/THRESHOLD),par[2]);
   return f;
}
