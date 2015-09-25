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

#include <TTree.h>
#include <TBranch.h>
#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>

#include <sstream>
#include <iostream>

#define GLOBAL_FIT 0	// Set to 1 for global fits, set to 0 for individual fits

// Declare files
TFile *infile;		// input ROOT file, 
TFile *outfile;		// output ROOT file

// Declare constants
const Int_t NMODULES = 8;		// Modules per arm
const Double_t THRESHOLD = 16; 		// discriminator threshold in ADC counts (40 mV)

// Declare trees
TTree *PSC_tree;

// Declare variables
Int_t PSC_tree_nentries;
Int_t psc_mod_l;	// left PSC module number
Int_t psc_mod_r;	// right PSC module number
Double_t tdc_l;		// left PSC tdc time
Double_t tdc_r;		// right PSC tdc time
Double_t rf_l;		// RF time associated with tdc_l
Double_t rf_r;		// RF time associated with tdc_r
Double_t pp_l;		// pulse peak associated with tdc_l
Double_t pp_r;		// pulse peak associated with tdc_r

double means[80] = {}; 		// used in tw_proj() to create TGraphError data points
double errors[100] = {};	// used in tw_proj() to create TGraphError error bars
double dtMeanL[NMODULES] = {};	// used to center the time dist. for the left modules
double dtMeanR[NMODULES] = {};	// used to center the time dist. for the right modules
double mpv[2][NMODULES] = {},{};// used to find the most probable value for the pulse peaks

#if GLOBAL_FIT == 1
// Global time walk fit parameters
double c0 = 0;
double c1 = 0;
double c2 = 0;
#endif

#if GLOBAL_FIT == 0
// Individual time walk fit parameters
double c0[2][NMODULES] = {},{};
double c1[2][NMODULES] = {},{};
double c2[2][NMODULES] = {},{};
#endif

// Declare 1D histograms
TH1F *h_dt_l[NMODULES];		// timing dist for each left module against the RF
TH1F *h_dt_r[NMODULES];		// timing dist for each right module against the RF
TH1F *h_dt_all;			// all timing dist. after being centered at 0
TH1F *h_dt_l_offset[NMODULES];	// timing dist for each left module after centered at 0
TH1F *h_dt_r_offset[NMODULES];	// timing dist for each right module after centered at 0
TH1F *h_dt_l_corr[NMODULES];	// timing dist for each left module after corrections
TH1F *h_dt_r_corr[NMODULES];	// timing dist for each right module after corrections
TH1F *h_dt_all_corr;		// timing dist for all modules after corrections
TH1F *h_pp_l[NMODULES];		// pulse peak for the left modules
TH1F *h_pp_r[NMODULES];		// pulse peak for the right modules

// Declare 2D histograms
TH2F *h_dt_vs_pp_l[NMODULES];		// dt vs pp for each left module
TH2F *h_dt_vs_pp_r[NMODULES];		// dt vs pp for each right module
TH2F *h_dt_vs_pp_l_corr[NMODULES];	// dt vs pp with corrections for each left module
TH2F *h_dt_vs_pp_r_corr[NMODULES];	// dt vs pp with corrections for each right module
TH2F *h_dt_vs_pp_all;			// dt vs pp for all modules
TH2F *h_dt_vs_pp_all_corr;		// dt vs pp with corrections for all modules

TGraphErrors *g_dt_vs_pp;		// dt vs pp plot converted into data points

void tw_plot() {
   // Get input file, make output file
   std::cout << "Getting input file" << std::endl;
   infile = new TFile("hd_root.root");
   std::cout << "Making output file" << std::endl;
   outfile = new TFile("results.root","RECREATE");

   // Create directories
   std::cout << "Creating output file directories" << std::endl;
   outfile->mkdir("dt_init");
   outfile->mkdir("dt_init/dt_init_L");
   outfile->mkdir("dt_init/dt_init_R");
   outfile->mkdir("dt_offset");
   outfile->mkdir("dt_offset/dt_offset_L");
   outfile->mkdir("dt_offset/dt_offset_R");
   outfile->mkdir("dt_corr");
   outfile->mkdir("dt_corr/dt_corr_L");
   outfile->mkdir("dt_corr/dt_corr_R");
   outfile->mkdir("pp");
   outfile->mkdir("pp/pp_L");
   outfile->mkdir("pp/pp_R");
   outfile->mkdir("dt_vs_pp");
   outfile->mkdir("dt_vs_pp/dt_vs_pp_L");
   outfile->mkdir("dt_vs_pp/dt_vs_pp_R");
   outfile->mkdir("dt_vs_pp_corr");
   outfile->mkdir("dt_vs_pp_corr/dt_vs_pp_corr_L");
   outfile->mkdir("dt_vs_pp_corr/dt_vs_pp_corr_R");
   outfile->mkdir("dt_vs_pp_graph");
   outfile->mkdir("dt_vs_pp_graph/dt_vs_pp_graph_L");
   outfile->mkdir("dt_vs_pp_graph/dt_vs_pp_graph_R");

   // Make canvas
   TCanvas *c = (TCanvas*)gROOT->FindObject("c");
   if (c == 0) {
      c = new TCanvas("c","c",0,20,600,500);
   }

   // Make histograms
   std::cout << "Naming histograms" << std::endl;
   for (Int_t i = 0; i < NMODULES; ++i) {
      h_dt_l[i] = new TH1F(Form("dt_l_%i",i+1),Form("Time difference;PSC_L_%i - RF_L (ns)",i+1),100,-5,5);
      h_dt_r[i] = new TH1F(Form("dt_r_%i",i+1),Form("Time difference;PSC_R_%i - RF_R (ns)",i+1),100,-5,5);
      h_dt_l_offset[i] = new TH1F(Form("dt_l_offset_%i",i+1),Form("Time difference after offset;PSC_L_%i - RF_L (ns)",i+1),100,-5,5);
      h_dt_r_offset[i] = new TH1F(Form("dt_r_offset_%i",i+1),Form("Time difference after offset;PSC_R_%i - RF_R (ns)",i+1),100,-5,5);
      h_dt_l_corr[i] = new TH1F(Form("dt_l_corr_%i",i+1),Form("Time difference after corrections;PSC_L_%i - RF_L (ns)",i+1),100,-5,5);
      h_dt_r_corr[i] = new TH1F(Form("dt_r_corr_%i",i+1),Form("Time difference after corrections;PSC_R_%i - RF_R (ns)",i+1),100,-5,5);
      h_pp_l[i] = new TH1F(Form("pp_l_%i",i+1),Form("Pulse Peak;PSC_L_%i Pulse Peak (adc counts)",i+1),1000,0,1000);
      h_pp_r[i] = new TH1F(Form("pp_r_%i",i+1),Form("Pulse Peak;PSC_R_%i Pulse Peak (adc counts)",i+1),1000,0,1000);
      h_dt_vs_pp_l[i] = new TH2F(Form("dt_vs_pp_l_%i",i+1),Form("dt vs pulse peak for PSC_L %i;Pulse peak (adc counts); PSC_L %i - RF_L (ns)",i+1,i+1),1000,0,1000,100,-5,5);
      h_dt_vs_pp_r[i] = new TH2F(Form("dt_vs_pp_r_%i",i+1),Form("dt vs pulse peak for PSC_R %i;Pulse peak (adc counts); PSC_R %i - RF_R (ns)",i+1,i+1),1000,0,1000,100,-5,5);
      h_dt_vs_pp_l_corr[i] = new TH2F(Form("dt_vs_pp_l_corr_%i",i+1),Form("dt vs pulse peak corrected for PSC_L %i;Pulse peak (adc counts); PSC_L %i - RF_L (ns)",i+1,i+1),1000,0,1000,100,-5,5);
      h_dt_vs_pp_r_corr[i] = new TH2F(Form("dt_vs_pp_r_corr_%i",i+1),Form("dt vs pulse peak corrected for PSC_R %i;Pulse peak (adc counts); PSC_R %i - RF_R (ns)",i+1,i+1),1000,0,1000,100,-5,5);
   }

   h_dt_all = new TH1F("dt_all","Time difference PSC - RF;PSC - RF (ns)",100,-5,5);
   h_dt_all_corr = new TH1F("dt_all_corr","Time difference corrected PSC - RF;PSC - RF (ns)",100,-5,5);
   h_dt_vs_pp_all = new TH2F("dt_vs_pp_all","dt vs pulse peak for all PSC;Pulse peak (adc counts);PSC - RF (ns)",1000,0,1000,100,-5,5);
   h_dt_vs_pp_all_corr = new TH2F("dt_vs_pp_all_corr","dt vs pulse peak corrected for PSC;Pulse peak (adc counts); PSC - RF (ns)",1000,0,1000,100,-5,5);

   
   // Make trees from input file
   PSC_tree = new TTree();
   PSC_tree = (TTree*)infile->Get("PSC_tree");

   // Grab branches
   PSC_tree->SetBranchAddress("tdc_l", &tdc_l);
   PSC_tree->SetBranchAddress("tdc_r", &tdc_r);
   PSC_tree->SetBranchAddress("rf_l", &rf_l);
   PSC_tree->SetBranchAddress("rf_r", &rf_r);
   PSC_tree->SetBranchAddress("pp_l", &pp_l);
   PSC_tree->SetBranchAddress("pp_r", &pp_r);
   PSC_tree->SetBranchAddress("psc_mod_l", &psc_mod_l);
   PSC_tree->SetBranchAddress("psc_mod_r", &psc_mod_r);

   // Get the number of entries
   PSC_tree_nentries = PSC_tree->GetEntries();

   // Get initial timing offsets
   std::cout << "Getting initial timing offsets" << std::endl;
   for (Int_t i = 0; i < PSC_tree_nentries; ++i) {
      PSC_tree->GetEntry(i);
      if (psc_mod_l < 1 || psc_mod_r < 1) continue;
      h_dt_l[psc_mod_l - 1]->Fill(tdc_l - rf_l);
      h_dt_r[psc_mod_r - 1]->Fill(tdc_r - rf_r);
   }
   std::cout << "Writing initial tming distributions" << std::endl;
   for (Int_t i = 0; i < NMODULES; ++i) {
      outfile->cd("dt_init/dt_init_L");
      h_dt_l[i]->Write();

      outfile->cd("dt_init/dt_init_R");
      h_dt_r[i]->Write();
   }

   std::cout << "Fitting initial timing distributions" << std::endl;
   for (Int_t i = 0; i < NMODULES; ++i) {
      TFitResultPtr resultpL0 = h_dt_l[i]->Fit("gaus","sq");
      double meanL = resultpL0->Parameter(1);
      double sigmaL = resultpL0->Parameter(2);
      // fit a second time using the mean and sigma of the initial fit
      TFitResultPtr resultpL = h_dt_l[i]->Fit("gaus","sq","",meanL-2*sigmaL,meanL+2*sigmaL);
      dtMeanL[i] = resultpL->Parameter(1);

      TFitResultPtr resultpR0 = h_dt_r[i]->Fit("gaus","sq");
      double meanR = resultpR0->Parameter(1);
      double sigmaR = resultpR0->Parameter(2);
      // fit a second time using the mean and sigma of the initial fit
      TFitResultPtr resultpR = h_dt_r[i]->Fit("gaus","sq","",meanR-2*sigmaR,meanR+2*sigmaR);
      dtMeanR[i] = resultpR->Parameter(1);

      h_dt_l[i]->Reset();
      h_dt_r[i]->Reset();
   }

   // Fill histograms with adjusted timing which centers around 0
   std::cout << "Adjusting time distributions to be centered about 0" << std::endl;
   for (Int_t i = 0; i < PSC_tree_nentries; ++i) {
      PSC_tree->GetEntry(i);
      if (psc_mod_l < 1 || psc_mod_r < 1) continue;

      // check for 4 ns wrap around due to RF, add 4 ns if necessary
      double tdiff_l = tdc_l - rf_l - dtMeanL[psc_mod_l - 1];
      double tdiff_r = tdc_r - rf_r - dtMeanR[psc_mod_r - 1];
      if (tdiff_l < -2)  tdiff_l += 4;
      if (tdiff_r < -2)  tdiff_r += 4;

      // left modules
      h_pp_l[psc_mod_l - 1]->Fill(pp_l);
      h_dt_l_offset[psc_mod_l - 1]->Fill(tdiff_l);
      h_dt_all->Fill(tdiff_l);
      h_dt_vs_pp_l[psc_mod_l - 1]->Fill(pp_l,tdiff_l);
      h_dt_vs_pp_all->Fill(pp_l,tdiff_l);
      
      // right modules
      h_pp_r[psc_mod_r - 1]->Fill(pp_r);
      h_dt_r_offset[psc_mod_r - 1]->Fill(tdiff_r);
      h_dt_all->Fill(tdiff_r);
      h_dt_vs_pp_r[psc_mod_r - 1]->Fill(pp_r,tdiff_r);
      h_dt_vs_pp_all->Fill(pp_r, tdiff_r);
   }

   // Get the most probable value of the pulse peak distribution for each module
   for (Int_t i = 0; i < NMODULES; ++i) {
      // left module
      TFitResultPtr ptrL1 = h_pp_l[i]->Fit("landau","sq");
      double mpv_L = ptrL1->Parameter(1);
      double sigma_pp_L = ptrL1->Parameter(2);
      TFitResultPtr ptrL2 = h_pp_l[i]->Fit("landau","sq","",mpv_L-1.5*sigma_pp_L,mpv_L+1.5*sigma_pp_L);
      mpv[0][i] = ptrL2->Parameter(1);

      // right module
      TFitResultPtr ptrR1 = h_pp_r[i]->Fit("landau","sq");
      double mpv_R = ptrR1->Parameter(1);
      double sigma_pp_R = ptrR1->Parameter(2);
      TFitResultPtr ptrR2 = h_pp_r[i]->Fit("landau","sq","",mpv_R-1.5*sigma_pp_R,mpv_R+1.5*sigma_pp_R);
      mpv[1][i] = ptrR2->Parameter(1);
   }

   std::cout << "Writing to file ..." << std::endl;
   for (Int_t i = 0; i < NMODULES; ++i) {
      outfile->cd("dt_offset/dt_offset_L");
      h_dt_l_offset[i]->Write();

      outfile->cd("dt_offset/dt_offset_R");
      h_dt_r_offset[i]->Write();

      outfile->cd("dt_offset");
      h_dt_all->Write();

      outfile->cd("pp/pp_L");
      h_pp_l[i]->Write();

      outfile->cd("pp/pp_R");
      h_pp_r[i]->Write();

      outfile->cd("dt_vs_pp/dt_vs_pp_L");
      h_dt_vs_pp_l[i]->Write();

      outfile->cd("dt_vs_pp/dt_vs_pp_R");
      h_dt_vs_pp_r[i]->Write();
   }
   outfile->cd("dt_vs_pp");
   h_dt_vs_pp_all->Write();

#if GLOBAL_FIT == 1
   std::cout << "Calling tw_proj()" << std::endl;
   tw_proj(h_dt_vs_pp_all);
#endif

#if GLOBAL_FIT == 0
   for (int j = 0; j < NMODULES; ++j) {
      std::cout << "Calling tw_fit() for arm" << i 
                << " and module " << j << std::endl;
      tw_proj(h_dt_vs_pp_l[j],0,j);
      tw_proj(h_dt_vs_pp_r[j],1,j);
   }
#endif

}

// Global fit tw_proj()
void tw_proj(TH2F *h_tw) {
   // Take a dt_vs_pp histogram and make y-axis projections (dt) and fit with a gaussian.
   // Use the means and errors to produce a TGraph with error bars

   std::cout << "Making projections..." << std::endl;

   TH1D *h_dt_psc[80];
   double bins[80] = {};
   double binErr[80] = {};

   for (int i = 0; i < 80; ++i) {
      bins[i] = (i+0.5)*10 + 200; // 200 based on the range of the PSC pulse heights
      binErr[i] = 0.25;
      errors[i] = 0.25;
      h_dt_psc[i] = new TH1D(Form("h_dt_psc_%i",i+1),Form("Y projection %i;dt (ns)",i+1),100,-5,5);
      h_dt_psc[i] = h_tw->ProjectionY(Form("h_dt_psc_%i",i+1),i*10+200,(i+1)*10+200);
      if (h_dt_psc[i]->GetEntries() < 20) {
         if (i < 30) {
            means[i] = 1;
            errors[i] = 1;
            continue;
         }
         else if (i > 50) {
            means[i] = means[i-13] + means[i-12] + means[i-11]; // predicts large pulse peaks
            means[i] /= 3;
            errors[i] = 0.5;
            continue;
         }
      }
      TFitResultPtr resultp = h_dt_psc[i]->Fit("gaus","sq","",-1.5,1.5);
      if (resultp->Parameter(1) > 1 || resultp->Parameter(1) < -1) {
         if (i < 30) {
            means[i] = 1;
            errors[i] = 1;
         }
         else if (i > 50) {
            means[i] = means[i-13] + means[i-12] + means[i-11]; // predicts large pulse peaks
            means[i] /= 3;
            errors[i] = 0.5;
         }
      }
      else {
         means[i] = resultp->Parameter(1);
         errors[i] = resultp->Parameter(2);
      }
   }

   std::cout << "Making TGraphErrors" << std::endl;
   g_dt_vs_pp = new TGraphErrors(80,bins,means,binErr,errors);
   g_dt_vs_pp->SetTitle("Time difference vs pulse peak");
   g_dt_vs_pp->GetXaxis()->SetTitle("Pulse peak (adc counts)");
   g_dt_vs_pp->GetYaxis()->SetTitle("dt (ns)");
   g_dt_vs_pp->GetXaxis()->CenterTitle();
   g_dt_vs_pp->GetYaxis()->CenterTitle();
   TAxis *axis = g_dt_vs_pp->GetXaxis();
   axis->SetLimits(-5,1005);			// x axis
   g_dt_vs_pp->GetHistogram()->SetMaximum(5.);	// y axis
   g_dt_vs_pp->GetHistogram()->SetMinimum(-5.);	// y axis
   g_dt_vs_pp->Draw("ap");

   std::cout << "Calling tw_fit()" << std::endl;
   tw_fit(g_dt_vs_pp);

   for (int i = 0; i < 80; ++i) {
      delete h_dt_psc[i];
   }
   delete g_dt_vs_pp;

}

// Individual fit tw_proj()
void tw_proj(TH2F *h_tw, int arm, int module) {
   // Take a dt_vs_pp histogram and make y-axis projections (dt) and fit with a gaussian.
   // Use the means and errors to produce a TGraph with error bars

   std::cout << "Making projections..." << std::endl;

   TH1D *h_dt_psc[80];
   double bins[80] = {};
   double binErr[80] = {};

   for (int i = 0; i < 80; ++i) {
      bins[i] = (i+0.5)*10 + 200;
      binErr[i] = 0.25;
      errors[i] = 0.25;
      h_dt_psc[i] = new TH1D(Form("h_dt_psc_%i",i+1),Form("Y projection %i;dt (ns)",i+1),100,-5,5);
      h_dt_psc[i] = h_tw->ProjectionY(Form("h_dt_psc_%i",i+1),i*10+200,(i+1)*10+200);
      if (h_dt_psc[i]->GetEntries() < 20) {
         if (i < 30) {
            means[i] = 1;
            errors[i] = 1;
            continue;
         }
         else if (i > 50) {
            means[i] = means[i-13] + means[i-12] + means[i-11]; // predicts large pulse peaks
            means[i] /= 3;
            errors[i] = 0.5;
            continue;
         }
      }
      TFitResultPtr resultp = h_dt_psc[i]->Fit("gaus","sq","",-1.5,1.5);
      if (resultp->Parameter(1) > 1 || resultp->Parameter(1) < -1) {
         if (i < 30) {
            means[i] = 1;
            errors[i] = 1;
         }
         else if (i > 50) {
            means[i] = means[i-13] + means[i-12] + means[i-11]; // predicts large pulse peaks
            means[i] /= 3;
            errors[i] = 0.5;
         }
      }
      else {
         means[i] = resultp->Parameter(1);
         errors[i] = resultp->Parameter(2);
      }
   }

   std::cout << "Making TGraphErrors" << std::endl;

   char a = "";
   if (arm == 0)
      a = 'L';
   else
      a = 'R';

   g_dt_vs_pp = new TGraphErrors(80,bins,means,binErr,errors);
   g_dt_vs_pp->SetTitle(Form("Time difference vs pulse peak for PSC_%c_%i",a,module+1));
   g_dt_vs_pp->GetXaxis()->SetTitle("Pulse peak (adc counts)");
   g_dt_vs_pp->GetYaxis()->SetTitle("dt (ns)");
   g_dt_vs_pp->GetXaxis()->CenterTitle();
   g_dt_vs_pp->GetYaxis()->CenterTitle();
   TAxis *axis = g_dt_vs_pp->GetXaxis();
   axis->SetLimits(-5,1005);			// x axis
   g_dt_vs_pp->GetHistogram()->SetMaximum(5.);	// y axis
   g_dt_vs_pp->GetHistogram()->SetMinimum(-5.); // y axis
   g_dt_vs_pp->Draw("ap");

   std::cout << "Calling tw_fit()" << std::endl;
   tw_fit(g_dt_vs_pp, arm, module);

   for (int i = 0; i < 80; ++i) {
      delete h_dt_psc[i];
   }
   delete g_dt_vs_pp;

}

// Global fit tw_fit()
void tw_fit(TGraphErrors *h_tw) {
   std::cout << "Making fit function" << std::endl;
   TF1 *ftw = new TF1("ftw",fitf,100.,1000.,3);		// fitf defined at the bottom
   std::cout << "Defined time-walk fit function" << std::endl;

   // Set initial parameters
   ftw->SetParameter(0,-1.9);
   ftw->SetParameter(1,4.8);
   ftw->SetParameter(2,-0.4);
   ftw->SetParName(0,"c0");
   ftw->SetParName(1,"c1");
   ftw->SetParName(2,"c2");

   std::cout << "Going to fit histograms" << std::endl;
   h_tw->Fit("ftw","RWq");
   std::cout << "Fit histogram" << std::endl;

   c0 = ftw->GetParameter("c0");
   c1 = ftw->GetParameter("c1");
   c2 = ftw->GetParameter("c2");
   std::cout << "c0: " << c0 << std::endl;
   std::cout << "c1: " << c1 << std::endl;
   std::cout << "c2: " << c2 << std::endl;

   // Write fit parameters to tw_parameters.out for use in CCDB
   std::ofstream flog("tw_parameters.out");
   for (Int_t i = 0; i < 2; ++i) {
      for (Int_t j = 0; j < NMODULES; ++j) {
         flog << i         << std::setw(15)
              << j+1       << std::setw(15)
              << c0        << std::setw(15)
              << c1        << std::setw(15)
              << c2        << std::setw(15)
              << THRESHOLD << std::setw(15)
              << mpv[i][j] << std::endl;
      }
   }

   outfile->cd("dt_vs_pp_graph");
   h_tw->Write();

}

// Individual tw_fit()
void tw_fit(TGraphErrors *h_tw, int arm, int module) {
   std::cout << "Making fit function" << std::endl;
   TF1 *ftw = new TF1("ftw",fitf,100.,1000.,3);		// fitf is defined at the bottom
   std::cout << "Defined time-walk fit function" << std::endl;

   ftw->SetParameter(0,-1.9);
   ftw->SetParameter(1,4.8);
   ftw->SetParameter(2,-0.4);
   ftw->SetParName(0,"c0");
   ftw->SetParName(1,"c1");
   ftw->SetParName(2,"c2");

   std::cout << "Going to fit histograms" << std::endl;
   h_tw->Fit("ftw","RWq");
   std::cout << "Fit histogram" << std::endl;

   c0[arm][module] = ftw->GetParameter("c0");
   c1[arm][module] = ftw->GetParameter("c1");
   c2[arm][module] = ftw->GetParameter("c2");
   std::cout << "c0: " << c0[arm][module] << std::endl;
   std::cout << "c1: " << c1[arm][module] << std::endl;
   std::cout << "c2: " << c2[arm][module] << std::endl;

   // Writing parameters to tw_parameters.out for use in CCDB
   std::ofstream flog("tw_parameters.out");
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

   if (arm == 0) {
      outfile->cd("dt_vs_pp_graph/dt_vs_pp_graph_L");
      stringstream title;
      title << "TGraph_L_" << module + 1;
      std::cout << "Saving " << title.str().c_str() << std::endl;
      h_tw->Write(title.str().c_str());
   }
   else {
      outfile->cd("dt_vs_pp_graph/dt_vs_pp_graph_R");
      stringstream title;
      title << "TGraph_R_" << module + 1;
      std::cout << "Saving " << title.str().c_str() << std::endl;
      h_tw->Write(title.str().c_str());
   }
}

void tw_corr() {
   std::cout << "Calling tw_plot()" << std::endl;
   tw_plot();

   std::cout << "Making corrected histograms" << std::endl;
   
   std::cout << "About to fill histograms" << std::endl;
   for (Int_t i = 0; i < PSC_tree_nentries; ++i) {
      PSC_tree->GetEntry(i);
      if (psc_mod_l < 1 || psc_mod_r < 1) continue;

      // check for 4 ns wrap and correct for it
      double tdiff_l = tdc_l - rf_l - dtMeanL[psc_mod_l - 1];
      double tdiff_r = tdc_r - rf_r - dtMeanR[psc_mod_r - 1];
      if (tdiff_l < -2)  tdiff_l += 4;
      if (tdiff_r < -2) tdiff_r += 4;

      // left modules
#if GLOBAL_FIT == 1
      tdiff_l = tdiff_l - c1*pow((pp_l/THRESHOLD),c2) + c1*pow(mpv[0][psc_mod_l - 1]/THRESHOLD,c2);
#else
      tdiff_l = tdiff_l - c1[0][psc_mod_l - 1]*pow((pp_l/THRESHOLD),c2[0][psc_mod_l - 1]) + c1[0][psc_mod_l - 1]*pow(mpv[0][psc_mod_l - 1]/THRESHOLD,c2[0][psc_mod_l - 1]);
#endif
      h_dt_l_corr[psc_mod_l - 1]->Fill(tdiff_l);
      h_dt_all_corr->Fill(tdiff_l);
      h_dt_vs_pp_l_corr[psc_mod_l - 1]->Fill(pp_l,tdiff_l);
      h_dt_vs_pp_all_corr->Fill(pp_l,tdiff_l);

      // right modules
#if GLOBAL_FIT == 1
      tdiff_r = tdiff_r - c1*pow((pp_r/THRESHOLD),c2) + c1*pow(mpv[1][psc_mod_r - 1]/THRESHOLD,c2);
#else
      tdiff_r = tdiff_r - c1[1][psc_mod_r - 1]*pow((pp_r/THRESHOLD),c2[1][psc_mod_r - 1]) + c1[1][psc_mod_r - 1]*pow(mpv[1][psc_mod_r - 1]/THRESHOLD,c2[1][psc_mod_r - 1]);
#endif
      h_dt_r_corr[psc_mod_r - 1]->Fill(tdiff_r);
      h_dt_all_corr->Fill(tdiff_r);
      h_dt_vs_pp_r_corr[psc_mod_r - 1]->Fill(pp_r,tdiff_r);
      h_dt_vs_pp_all_corr->Fill(pp_r,tdiff_r);
   }
   for (Int_t i = 0; i < NMODULES; ++i) {
      outfile->cd("dt_corr/dt_corr_L");
      h_dt_l_corr[i]->Write();

      outfile->cd("dt_corr/dt_corr_R");
      h_dt_r_corr[i]->Write();

      outfile->cd("dt_corr");
      h_dt_all_corr->Write();

      outfile->cd("dt_vs_pp_corr/dt_vs_pp_corr_L");
      h_dt_vs_pp_l_corr[i]->Write();

      outfile->cd("dt_vs_pp_corr/dt_vs_pp_corr_R");
      h_dt_vs_pp_r_corr[i]->Write();
   }
   outfile->cd("dt_vs_pp_corr");
   h_dt_vs_pp_all_corr->Write();

   double avg = 0;
   double avgL = 0;
   double avgR = 0;

   for (Int_t i = 0; i < NMODULES; ++i) {
      TFitResultPtr ptrL = h_dt_l_corr[i]->Fit("gaus","sq");
      std::cout << "Sigma for left module " << i+1 << ": " << (ptrL->Parameter(2))*1000 << std::endl;
      TFitResultPtr ptrR = h_dt_r_corr[i]->Fit("gaus","sq");
      std::cout << "Sigma for right module " << i+1 << ": " << (ptrR->Parameter(2))*1000 << std::endl;
      avgL += (ptrL->Parameter(2))*1000;
      avgR += (ptrR->Parameter(2))*1000;
      avg += (ptrL->Parameter(2) + ptrR->Parameter(2))*1000;
   }
   avgL /= 8;
   avgR /= 8;
   avg  /= 16;

   std::cout << "avgL: " << avgL << std::endl;
   std::cout << "avgR: " << avgR << std::endl;
   std::cout << "avg: "  << avg << std::endl;

   delete c;
}

Double_t fitf(Double_t *x, Double_t *par) {
   Float_t xx = x[0];
   Double_t f = par[0] + par[1]*pow((xx/THRESHOLD),par[2]);
   return f;
}
