#include <TProfile.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TPaveLabel.h>

#include <sstream>
#include <iostream>

TFile *f;
TH2I *h_tw[102];
TH2I *h_tw_ind[5][4];

void display(char const *inputFile) {
   
   TCanvas *c1[25];
   for (Int_t i = 0; i < 25; ++i)
   {
      c1[i] = new TCanvas(Form("c%i",i),Form("c%i",i),1200,800);
      gStyle->SetOptStat(0);
      c1[i]->Divide(2,3,0.0001,0.0001);
   }

   f = new TFile(inputFile);
   std::cout << "file: " << inputFile << std::endl;

   int canvas = 0;
   for (Int_t i = 0; i < 102; ++i)
   {
      h_tw[i] = (TH2I*)f->Get(Form("TAGM_TW/t-rf/h_dt_vs_pp_%i",i+1));
      h_tw[i]->SetTitle(Form("Timewalk Col %i",i+1));
      canvas = i/5;
      int cell = i%5 + 1;
      c1[canvas]->cd(cell);
      h_tw[i]->Draw("colz");
   }
   for (Int_t i = 0; i < 5; ++i)
   {
      for (Int_t j = 0; j < 4; ++j)
      {
         h_tw_ind[i][j] = (TH2I*)f->Get(Form("TAGM_TW/t-rf/h_dt_vs_pp_ind_%i_%i",i+1,j+1));
         h_tw_ind[i][j]->SetTitle(Form("Timewalk ind. Row %i Col %i",i+1,j+1));
         canvas = 21 + j;
         int cell = i+1;
         c1[canvas]->cd(cell);
         h_tw_ind[i][j]->Draw("colz");
      }
   }
}      
