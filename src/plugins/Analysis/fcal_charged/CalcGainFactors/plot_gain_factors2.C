void plot_gain_factors2(void)
{
// test program to read a file with a list of file names
//

#include <TRandom.h>

gROOT->Reset();
//TTree *Bfield = (TTree *) gROOT->FindObject("Bfield");
gStyle->SetPalette(1,0);
gStyle->SetOptStat(kTRUE);
gStyle->SetOptStat(11111111);
gStyle->SetOptFit(kTRUE);
gStyle->SetPadRightMargin(0.15);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
//

   char string[256];
   Int_t j,jj;
   Int_t ndx;
   Int_t nbins=100;

   const int nchan=2800;

   Double_t g2931[nchan];
   Double_t g3079[nchan];
   Double_t g3179[nchan];
   Double_t g3180[nchan];
   Double_t g3182[nchan];
   Double_t g2931_r20_60[nchan];
   Double_t g3079_r20_60[nchan];
   Double_t g3179_r20_60[nchan];
   Double_t g3180_r20_60[nchan];
   Double_t g3182_r20_60[nchan];

   Double_t g3079_3182_r20[nchan];
   Double_t g3079_3182_r20_60[nchan];
   Double_t g2931_3182_r20[nchan];
   Double_t g2931_3182_r20_60[nchan];

   // Look at runs 2931-3182;

   sprintf(string,"GainFactorIter1_2931-3182_r20-60.root");
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered2931_3182_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx2931_3182_r20_60 = (TH2D*)in->Get("h_yvsx");

   ifstream gainfile;
   ofstream outgain;
   outgain.open ("plot_gain_factors2.list");

   // Repeat for Run2931-3182
      
   gainfile.close();

   TString infile = "GainFactorIter1_2931-3182_r20-60.txt";
   cout << "Opening file: " << infile.Data() << endl;
   gainfile.open (infile.Data());
   if (!gainfile) {
       cout << "ERROR: Failed to open data file= " << infile.Data() << endl;
       return;
       }

   TString line;
   while (line.ReadLine(gainfile)){

     TObjArray *tokens = line.Tokenize(" ");
     Float_t gain = (((TObjString*)tokens->At(0))->GetString()).Atof();
     TString s = ((TObjString*)tokens->At(1))->GetString().Remove(0,8);
     Int_t channel = s.Atoi();
     g2931_3182_r20_60[channel] = gain;
     outgain << " R3182_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }      
   gainfile.close();



   TCanvas *c10 = new TCanvas("c10","c10 plot_gain_factors2",200,10,1000,500);
   c10->Divide(2,1);
   // c10->SetGridx();
   // c10->SetGridy();
   // c10->SetLogz();

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c10->cd(1);
   runnum = 3182;

   sprintf(string,"Electron Gain Factors Run 2931-3182, Jul27");
   h_yvsx2931_3182_r20_60->SetTitle(string);
   h_yvsx2931_3182_r20_60->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx2931_3182_r20_60->GetYaxis()->SetTitle("y");
   h_yvsx2931_3182_r20_60->GetXaxis()->SetTitle("x");
   h_yvsx2931_3182_r20_60->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c10->cd(2);

   h_unordered2931_3182_r20_60->SetTitle(string);;
   h_unordered2931_3182_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered2931_3182_r20_60->GetYaxis()->SetTitle("Counts");
   h_unordered2931_3182_r20_60->GetXaxis()->SetTitle("Gain Factor");
   h_unordered2931_3182_r20_60->Fit("gaus");
   h_unordered2931_3182_r20_60->Draw("");

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c10->cd(3);

   Float_t ymin=0;
   Float_t ymax=150;
   c10->cd(4);

   // print histograms

   sprintf(string,"plot_gain_factors2.pdf");
   c10->SaveAs(string);

   // close files

   outgain.close();

   // in->Close();


}

