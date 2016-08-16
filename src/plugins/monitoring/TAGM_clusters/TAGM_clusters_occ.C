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
TH1I *h_occ_b;
TH1I *h_occ_a;

void TAGM_clusters_occ(char const *inputFile) {
   
   TCanvas *c1;
   c1 = new TCanvas("c1","c1",600,500);
   gStyle->SetOptStat(0);

   f = new TFile(inputFile);

   h_occ_b = (TH1I*)f->Get("Before/occupancy_b")->Clone();
   h_occ_a = (TH1I*)f->Get("After/occupancy_a")->Clone();

   h_occ_b->SetTitle("Occupancy before - after");
   h_occ_b->GetXaxis()->SetTitle("Column number");

   h_occ_b->Add(h_occ_a,-1);
   h_occ_b->Draw();

}      
