// This macro will plot the difference between before and after merging
// for the occupancy and energy plots

#include <TH1.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TPaveText.h>
#include <TPaveLabel.h>

#include <sstream>
#include <iostream>

TFile *f;
TH1I *hEb;
TH1I *hEa;

void TAGM_clusters_E(char const *inputFile) {
   
   TCanvas *c1;
   c1 = new TCanvas("c1","c1",600,500);
   gStyle->SetOptStat(0);

   f = new TFile(inputFile);

   hEb = (TH1I*)f->Get("Before/E_b")->Clone();
   hEa = (TH1I*)f->Get("After/E_a")->Clone();

   hEb->SetTitle("Energy occupancy before - after");
   hEb->GetXaxis()->SetTitle("Energy (GeV)");

   hEb->Add(hEa,-1);
   hEb->Draw();

}      
