// ROOT macro to grab a histogram and make projections per bin
// Created by Alex Barnes, 6/30/2015

#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TProfile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TMath.h>

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

#define CORRECTIONS false	// when true the correction parameters will be applied
#define TAGM true		// include TAGM when true 
#define TAGH false		// include TAGH when true

// Declare constants
#if TAGM
const int MAX_COLUMNS			= 100;	// Number of TAGM channels
int bad_channels_tm			= 0;	// Number of channels excluded in fit
double fitResults_tm[MAX_COLUMNS][3]	= {},{};// fit parameters, 3 pol2 params
double max_E_tm[MAX_COLUMNS]		= {};	// maximum of fit
double fit_tm[3]			= {};	// average fit parameters
#endif

#if TAGH
const int MAX_COUNTERS			= 274;	// Number of TAGH channels
int bad_channels_th			= 0;	// Number of channels excluded in fit
double fitResults_th[MAX_COUNTERS][3]	= {},{};// fit parameters, 3 pol2 params
double max_E_th[MAX_COUNTERS]		= {};	// maximum of fit
double fit_th[3]			= {};	// average fit parameters
#endif

double Ebw_PS = 0.013;				// PS energy bin width
double Ebl_PS = 2.3;				// PS low energy
double Ebh_PS = 4.9;				// PS high energy
double NEb_PS = (Ebh_PS - Ebl_PS)/Ebw_PS;	// Number of bins for PS energy

void PSEcorr(char* inputFile) {
   TCanvas *c1 = (TCanvas*)gROOT->FindObject("c1");
   if (c1 == 0)
      c1 = new TCanvas("c1","c1",0,20,600,500);
   c1->cd();

   std::cout << "Getting ROOT file" << std::endl;
   TFile *infile = new TFile(inputFile);
   TFile *outfile = new TFile("profiles.root", "recreate");

   // Create directories
   std::cout << "Creating output file directories" << std::endl;
   outfile->mkdir("TAGM");
   outfile->mkdir("TAGH");

#if TAGM

#if !CORRECTIONS
   std::ofstream fout_tm("Eparms-TAGM.out");
#endif // !CORRECTIONS

#if CORRECTIONS
   std::fstream parmsFile("Eparms-TAGM.out", ios_base::in);
   string parm0, parm1, parm2;
   double p0, p1, p2;
   int index = 0;
   if (parmsFile) {
      std::string line;
      while (parmsFile >> std::ws && std::getline(parmsFile, line))
         ;
      std::istringstream input(line);
      input >> parm0 >> parm1 >> parm2;
      p0 = atof(parm0.c_str());
      p1 = atof(parm1.c_str());
      p2 = atof(parm2.c_str());
   }
#endif // CORRECTIONS

   std::cout << "Getting original TAGM histograms" << std::endl;

   TH2F *h_psE_vs_psEl_tm[MAX_COLUMNS];
   TProfile *psE_vs_psEl_tm[MAX_COLUMNS];
   TProfile *psE_vs_psEl_tm_corr[MAX_COLUMNS];
   TH1D *h_tmp;
   for (int col = 0; col < MAX_COLUMNS; ++col) {
      h_psE_vs_psEl_tm[col] = (TH2F*)infile->Get(Form("TAGM/h_psE_vs_psEl_tm_%i",col+1));
      psE_vs_psEl_tm_corr[col] = new TProfile(Form("psE_vs_psEl_tm_corr_%i",col+1),
                                              Form("PS E vs PS left, TAGM col %i;\
                                              Energy asymmetry;PS energy (GeV)",col+1),
                                              50,0,1,2.3,4.9);
      psE_vs_psEl_tm_corr[col]->SetMaximum(4.9);
      psE_vs_psEl_tm_corr[col]->SetMinimum(2.3);
      h_tmp = h_psE_vs_psEl_tm[col]->ProjectionY();
      int bMax = h_tmp->GetMaximumBin();
      int yMax = bMax;
      int yMin = bMax;
      for (int i = 0; i < (NEb_PS-bMax); ++i) {
         if (h_tmp->GetBinContent(bMax+i) < 2) {
            yMax = bMax + i;
            break;
         }
      }
      for (int i = 0; i < bMax; ++i) {
         if (h_tmp->GetBinContent(bMax-i) < 2) {
            yMin = bMax - i;
            break;
         }
      }
      psE_vs_psEl_tm[col] = h_psE_vs_psEl_tm[col]->ProfileX(Form("col_%i",col+1),yMin,yMax);
      int bad_bins = 0;
      for (int i = 1; i <= 50; ++i) {
         double bEntries = psE_vs_psEl_tm[col]->GetBinEntries(i);
         if (bEntries < 5) {
            psE_vs_psEl_tm[col]->SetBinEntries(i,0);
            bad_bins++;
         }
      }
      int nentries = psE_vs_psEl_tm[col]->GetEntries();
      if (nentries == 0 || bad_bins > 45) {
         for (int i = 0; i < 3; ++i) {
            fitResults_tm[col][i] = 0;
            bad_channels_tm++;
         }
      }
      else {
#if CORRECTIONS
         for (int i = 1; i <= 50; ++i) {
            double content = psE_vs_psEl_tm[col]->GetBinContent(i);
            double asym = (i+0.5)*0.02;
            //content /= (p0 + p1*asym + p2*asym*asym);
            psE_vs_psEl_tm_corr[col]->Fill(asym,content);
         }
#endif // CORRECTIONS

#if !CORRECTIONS
         TFitResultPtr ptr = psE_vs_psEl_tm[col]->Fit("pol2","sq");
         fitResults_tm[col][0] = ptr->Parameter(0);
         fitResults_tm[col][1] = ptr->Parameter(1);
         fitResults_tm[col][2] = ptr->Parameter(2);

         // maximum of parabola at c0 - c1^2/(4c2)
         max_E_tm[col] = -1*fitResults_tm[col][1]*fitResults_tm[col][1];
         max_E_tm[col] /= 4*fitResults_tm[col][2];
         max_E_tm[col] += fitResults_tm[col][0];

         for (int i = 0; i < 3; ++i) {
            fitResults_tm[col][i] /= max_E_tm[col];
            fit_tm[i] += fitResults_tm[col][i];
         }
#endif // !CORRECTIONS
      }

#if CORRECTIONS
      outfile->cd("TAGM");
      psE_vs_psEl_tm_corr[col]->Write();
#endif // CORRECTIONS

#if !CORRECTIONS
      outfile->cd("TAGM");
      psE_vs_psEl_tm[col]->Write();
#endif // !CORRECTIONS
   }
#if !CORRECTIONS
   for (int i = 0; i < 3; ++i) {
      fit_tm[i] /= (MAX_COLUMNS - bad_channels_tm);
   }

   fout_tm << fit_tm[0] << std::setw(15)
           << fit_tm[1] << std::setw(15)
           << fit_tm[2] << std::endl;

#endif // !CORRECTIONS
#endif // TAGM

#if TAGH

#if !CORRECTIONS
   std::ofstream fout_th("Eparms-TAGH.out");
#endif // !CORRECTIONS

#if CORRECTIONS
   std::fstream parmsFile("Eparms-TAGH.out", ios_base::in);
   string parm0, parm1, parm2;
   double p0, p1, p2;
   int index = 0;
   if (parmsFile) {
      std::string line;
      while (parmsFile >> std::ws && std::getline(parmsFile, line))
         ;
      std::istringstream input(line);
      input >> parm0 >> parm1 >> parm2;
      p0 = atof(parm0.c_str());
      p1 = atof(parm1.c_str());
      p2 = atof(parm2.c_str());
   }
   TProfile *psE_vs_psEl_th_corr[MAX_COUNTERS];
#endif // CORRECTIONS

   std::cout << "Getting original TAGH histograms" << std::endl;

   TH2F *h_psE_vs_psEl_th[MAX_COUNTERS];
   TH2F *h_psE_vs_psEl_th_side[MAX_COUNTERS];
   TProfile *psE_vs_psEl_th[MAX_COUNTERS];
   TH1D *h_tmp;
   for (int hodo = 0; hodo < MAX_COUNTERS; ++hodo) {
      h_psE_vs_psEl_th[hodo] = (TH2F*)infile->Get(Form("TAGH/psE_vs_psEl_th_%i",hodo+1));
      h_tmp = h_psE_vs_psEl_th[hodo]->ProjectionY();
      int bMax = h_tmp->GetMaximumBin();
      int yMax = bMax;
      int yMin = bMax;
      for (int i = 0; i < (NEb_PS-bMax); ++i) {
         if (h_tmp->GetBinContent(bMax+i) < 2) {
            yMax = bMax + i;
            break;
         }
      }
      for (int i = 0; i < bMax; ++i) {
         if (h_tmp->GetBinContent(bMax-i) < 2) {
            yMin = bMax - i;
            break;
         }
      }
      psE_vs_psEl_th[hodo] = h_psE_vs_psEl_th[hodo]->ProfileX(Form("hodo_%i",hodo+1),yMin,yMax);
      int bad_bins = 0;
      for (int i = 1; i <= 50; ++i) {
         double bEntries = psE_vs_psEl_th[hodo]->GetBinEntries(i);
         if (bEntries < 5) {
            psE_vs_psEl_th[hodo]->SetBinEntries(i,0);
            bad_bins++;
         }
      }
      int nentries = psE_vs_psEl_th[hodo]->GetEntries();
      if (nentries == 0 || bad_bins > 45) {
         for (int i = 0; i < 3; ++i) {
            fitResults_th[hodo][i] = 0;
            bad_channels_th++;
         }
      }
      nentries = 0;
      else {
#if CORRECTIONS
         psE_vs_psEl_th_corr[hodo] = new TProfile(Form("psE_vs_psEl_th_corr_%i",hodo+1),
                                                  Form("PS E vs PS left, TAGH counter %i;\
                                                  Energy asymmetry;PS energy (GeV)",hodo+1),
                                                  50,0,1,2.3,4.9);
         psE_vs_psEl_th_corr[hodo]->SetMaximum(4.9);
         psE_vs_psEl_th_corr[hodo]->SetMinimum(2.3);
         for (int i = 1; i <= 50; ++i) {
            double content = psE_vs_psEl_th[hodo]->GetBinContent(i);
            double asym = (i+0.5)*0.02;
            content /= (p0 + p1*asym + p2*asym*asym);
            psE_vs_psEl_th_corr[hodo]->Fill(asym,content);
         }
#endif // CORRECTIONS

#if !CORRECTIONS
         TFitResultPtr ptr = psE_vs_psEl_th[hodo]->Fit("pol2","sq");
         fitResults_th[hodo][0] = ptr->Parameter(0);
         fitResults_th[hodo][1] = ptr->Parameter(1);
         fitResults_th[hodo][2] = ptr->Parameter(2);

         // maximum of parabola at c0 - c1^2/(4c2)
         max_E_th[hodo] = -1*fitResults_th[hodo][1]*fitResults_th[hodo][1];
         max_E_th[hodo] /= 4*fitResults_th[hodo][2];
         max_E_th[hodo] += fitResults_th[hodo][0];

         for (int i = 0; i < 3; ++i) {
            fitResults_th[hodo][i] /= max_E_th[hodo];
            fit_th[i] += fitResults_th[hodo][i];
         }
#endif // !CORRECTIONS
      }

#if CORRECTIONS
      outfile->cd("TAGH");
      psE_vs_psEl_th_corr[hodo]->Write();
#endif // CORRECTIONS

#if !CORRECTIONS
      outfile->cd("TAGH");
      h_psE_vs_psEl_th[hodo]->Write();
      psE_vs_psEl_th[hodo]->Write();
#endif // !CORRECTIONS
   }

#if !CORRECTIONS
   std::cout << "bad channels " << bad_channels_th << std::endl;

   for (int i = 0; i < 3; ++i) {
      fit_th[i] /= (MAX_COUNTERS - bad_channels_th);
   }

   fout_th << fit_th[0] << std::setw(15)
           << fit_th[1] << std::setw(15) 
           << fit_th[2] << std::endl;

#endif // !CORRECTIONS
#endif // TAGH
}

