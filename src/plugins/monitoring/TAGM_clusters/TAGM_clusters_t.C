// This macro displays the pre-merge time difference
// between adjacent columns

#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TPaveLabel.h>

#include <sstream>
#include <iostream>

TFile *f;
TH1I *h_t[10];

void TAGM_clusters_t(char* inputFile) {
   
   TCanvas *c1;
   c1 = new TCanvas("c1","c1",1200,800);
   gStyle->SetOptStat(0);
   c1->Divide(5,2);

   f = new TFile(inputFile);

   for (Int_t i = 0; i < 10; ++i)
   {
      h_t[i] = (TH1I*)f->Get(Form("Before/deltaT_%i",i+1));
      h_t[i]->SetTitle(Form("delta t, %i-%i",10*i+1,10*i+10));
      h_t[i]->GetXaxis()->SetTitle("#deltat (ns)");
      c1->cd(i+1);
      h_t[i]->Draw();
   }
}      
