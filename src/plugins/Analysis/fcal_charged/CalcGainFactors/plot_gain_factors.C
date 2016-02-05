void plot_gain_factors(void)
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

   // read histograms from file

   Int_t runnum=2931;

   sprintf(string,"GainFactorIter1_%4d_r20.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered2931_r20 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx2931_r20 = (TH2D*)in->Get("h_yvsx");

   sprintf(string,"GainFactorIter1_%4d_r20-60.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered2931_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx2931_r20_60 = (TH2D*)in->Get("h_yvsx");

   Int_t runnum=3079;

   sprintf(string,"GainFactorIter1_%4d_r20.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3079_r20 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3079_r20 = (TH2D*)in->Get("h_yvsx");

   sprintf(string,"GainFactorIter1_%4d_r20-60.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3079_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3079_r20_60 = (TH2D*)in->Get("h_yvsx");

   Int_t runnum=3179;

   sprintf(string,"GainFactorIter1_%4d_r20.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3179_r20 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3179_r20 = (TH2D*)in->Get("h_yvsx");

   sprintf(string,"GainFactorIter1_%4d_r20-60.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3179_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3179_r20_60 = (TH2D*)in->Get("h_yvsx");

   Int_t runnum=3180;

   sprintf(string,"GainFactorIter1_%4d_r20.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3180_r20 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3180_r20 = (TH2D*)in->Get("h_yvsx");

   sprintf(string,"GainFactorIter1_%4d_r20-60.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3180_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3180_r20_60 = (TH2D*)in->Get("h_yvsx");

   Int_t runnum=3182;

   sprintf(string,"GainFactorIter1_%4d_r20.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3182_r20 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3182_r20 = (TH2D*)in->Get("h_yvsx");

   sprintf(string,"GainFactorIter1_%4d_r20-60.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3182_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3182_r20_60 = (TH2D*)in->Get("h_yvsx");

   // Look at runs 3079-3182;

   sprintf(string,"GainFactorIter1_3079-3182_r20.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3079_3182_r20 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3079_3182_r20 = (TH2D*)in->Get("h_yvsx");

   sprintf(string,"GainFactorIter1_3079-3182_r20-60.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered3079_3182_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx3079_3182_r20_60 = (TH2D*)in->Get("h_yvsx");

   // Look at runs 2931-3182;

   sprintf(string,"GainFactorIter1_2931-3182_r20.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered2931_3182_r20 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx2931_3182_r20 = (TH2D*)in->Get("h_yvsx");

   sprintf(string,"GainFactorIter1_2931-3182_r20-60.root",runnum);
   printf ("Histogram input filename=%s\n",string);
   TFile *in = new TFile(string,"read");

   TH1D *h_unordered2931_3182_r20_60 = (TH1D*)in->Get("h_unordered");
   TH2D *h_yvsx2931_3182_r20_60 = (TH2D*)in->Get("h_yvsx");

   // define correlation histograms

   TH2D *h3079_2931 = new TH2D ("h3079_2931","Gains R3079 vs R2931",nbins,0,2,nbins,0,2);
   TH2D *h3179_2931 = new TH2D ("h3179_2931","Gains R3179 vs R2931",nbins,0,2,nbins,0,2);
   TH2D *h3180_2931 = new TH2D ("h3180_2931","Gains R3180 vs R2931",nbins,0,2,nbins,0,2);
   TH2D *h3182_2931 = new TH2D ("h3182_2931","Gains R3182 vs R2931",nbins,0,2,nbins,0,2);

   TH2D *h3079_2931_r20_60 = new TH2D ("h3079_2931_r20_60","Gains R3079 vs R2931 _r20_60",nbins,0,2,nbins,0,2);
   TH2D *h3179_2931_r20_60 = new TH2D ("h3179_2931_r20_60","Gains R3179 vs R2931 _r20_60",nbins,0,2,nbins,0,2);
   TH2D *h3180_2931_r20_60 = new TH2D ("h3180_2931_r20_60","Gains R3180 vs R2931 _r20_60",nbins,0,2,nbins,0,2);
   TH2D *h3182_2931_r20_60 = new TH2D ("h3182_2931_r20_60","Gains R3182 vs R2931 _r20_60",nbins,0,2,nbins,0,2);

   TH2D *h3079_3182_2931_r20 = new TH2D ("3079_3182_2931_r20","Gains R3079-3182 vs R2931 _r20",nbins,0,2,nbins,0,2);
   TH2D *h3079_3182_2931_r20_60 = new TH2D ("3079_3182_2931_r20_60","Gains R3079-3182 vs R2931 _r20_60",nbins,0,2,nbins,0,2);

   TH1D *gRatio_r20 = new TH1D ("gRatio","Ratio of R3079-3182/R2931 _r20",nbins,0,2);
   TH1D *gRatio_r20_60 = new TH1D ("gRatio","Ratio of R3079-3182/R2931 _r20_60",nbins,0,2);

   ifstream gainfile;
   ofstream outgain;
   outgain.open ("plot_gain_factors.list");

   TString infile = "GainFactorIter1_2931_r20.txt";
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
     g2931[channel] = gain;
     outgain << " R2931, channel=" << channel << " gain=" << gain  << endl;  
   }  
   gainfile.close();

   TString infile = "GainFactorIter1_3079_r20.txt";
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
     g3079[channel] = gain;
     outgain << " R3079, channel=" << channel << " gain=" << gain  << endl;
   
   }  
   gainfile.close();

   TString infile = "GainFactorIter1_3179_r20.txt";
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
     g3179[channel] = gain;
     outgain << " R3179, channel=" << channel << " gain=" << gain  << endl;
   
   }     
   gainfile.close();

   TString infile = "GainFactorIter1_3180_r20.txt";
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
     g3180[channel] = gain;
     outgain << " R3180, channel=" << channel << " gain=" << gain  << endl;
   
   }       
   gainfile.close();

   TString infile = "GainFactorIter1_3182_r20.txt";
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
     g3182[channel] = gain;
     outgain << " R3182, channel=" << channel << " gain=" << gain  << endl;
   
   }

   // Repeat for 20_60
      
   gainfile.close();

   TString infile = "GainFactorIter1_2931_r20-60.txt";
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
     g2931_r20_60[channel] = gain;
     outgain << " R2931_r20_60, channel=" << channel << " gain=" << gain  << endl;  
   }  
   gainfile.close();

   TString infile = "GainFactorIter1_3079_r20-60.txt";
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
     g3079_r20_60[channel] = gain;
     outgain << " R3079_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }  
   gainfile.close();

   TString infile = "GainFactorIter1_3179_r20-60.txt";
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
     g3179_r20_60[channel] = gain;
     outgain << " R3179_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }     
   gainfile.close();

   TString infile = "GainFactorIter1_3180_r20-60.txt";
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
     g3180_r20_60[channel] = gain;
     outgain << " R3180_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }       
   gainfile.close();

   TString infile = "GainFactorIter1_3182_r20-60.txt";
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
     g3182_r20_60[channel] = gain;
     outgain << " R3182_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }  


   // Repeat for Run3079-3182
      
   gainfile.close();

   TString infile = "GainFactorIter1_3079-3182_r20-60.txt";
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
     g3079_3182_r20_60[channel] = gain;
     outgain << " R3182_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }      
   gainfile.close();

   TString infile = "GainFactorIter1_3079-3182_r20.txt";
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
     g3079_3182_r20[channel] = gain;
     outgain << " R3182_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }      
     

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

   TString infile = "GainFactorIter1_2931-3182_r20.txt";
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
     g2931_3182_r20[channel] = gain;
     outgain << " R3182_r20_60, channel=" << channel << " gain=" << gain  << endl;
   
   }      
     



   // histogram gains

   for (Int_t jj=0;jj<nchan;jj++) {
     if (g2931[jj]!=1 || g3079[jj]!=1) h3079_2931->Fill(g2931[jj],g3079[jj]);
     if (g2931[jj]!=1 || g3179[jj]!=1) h3179_2931->Fill(g2931[jj],g3179[jj]);
     if (g2931[jj]!=1 || g3180[jj]!=1) h3180_2931->Fill(g2931[jj],g3180[jj]);
     if (g2931[jj]!=1 || g3182[jj]!=1) h3182_2931->Fill(g2931[jj],g3182[jj]);

     // cout << " jj=" << jj << " g3182_r20_60=" << g3182_r20_60[jj] << endl;

     if (g2931_r20_60[jj]!=1 || g3079_r20_60[jj]!=1) h3079_2931_r20_60->Fill(g2931_r20_60[jj],g3079_r20_60[jj]);
     if (g2931_r20_60[jj]!=1 || g3179_r20_60[jj]!=1) h3179_2931_r20_60->Fill(g2931_r20_60[jj],g3179_r20_60[jj]);
     if (g2931_r20_60[jj]!=1 || g3180_r20_60[jj]!=1) h3180_2931_r20_60->Fill(g2931_r20_60[jj],g3180_r20_60[jj]);
     if (g2931_r20_60[jj]!=1 || g3182_r20_60[jj]!=1) h3182_2931_r20_60->Fill(g2931_r20_60[jj],g3182_r20_60[jj]);

     if (g2931[jj]!=1 || g3079_3182_r20[jj]!=1) h3079_3182_2931_r20->Fill(g2931[jj],g3079_3182_r20[jj]);
     if (g2931_r20_60[jj]!=1 || g3079_3182_r20_60[jj]!=1) h3079_3182_2931_r20_60->Fill(g2931_r20_60[jj],g3079_3182_r20_60[jj]);

     if (g2931[jj]!=1 && g3079_3182_r20[jj]!=1) gRatio_r20->Fill(g3079_3182_r20[jj]/g2931[jj]);
     if (g2931_r20_60[jj]!=1 && g3079_3182_r20_60[jj]!=1) gRatio_r20_60->Fill(g3079_3182_r20_60[jj]/g2931_r20_60[jj]);
   }
    

   // plot histograms

   TCanvas *c1 = new TCanvas("c1","c1 plot_gain_factors",200,10,700,700);
   c1->Divide(2,2);
   // c1->SetGridx();
   // c1->SetGridy();
   // c1->SetLogz();
   Int_t runnum=2931;

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c1->cd(1);

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx2931_r20_60->SetTitle(string);
   h_yvsx2931_r20_60->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx2931_r20_60->GetYaxis()->SetTitle("y");
   h_yvsx2931_r20_60->GetXaxis()->SetTitle("x");
   h_yvsx2931_r20_60->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c1->cd(2);

   h_unordered2931_r20_60->SetTitle(string);;
   h_unordered2931_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered2931_r20_60->GetYaxis()->SetTitle("Counts");
   h_unordered2931_r20_60->GetXaxis()->SetTitle("Gain Factor");
   h_unordered2931_r20_60->Fit("gaus");
   h_unordered2931_r20_60->Draw("");

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c1->cd(3);

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx2931_r20->SetTitle(string);
   h_yvsx2931_r20->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx2931_r20->GetYaxis()->SetTitle("y");
   h_yvsx2931_r20->GetXaxis()->SetTitle("x");
   h_yvsx2931_r20->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c1->cd(4);

   h_unordered2931_r20->SetTitle(string);;
   h_unordered2931_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered2931_r20->GetYaxis()->SetTitle("Counts");
   h_unordered2931_r20->GetXaxis()->SetTitle("Gain Factor");
   h_unordered2931_r20->Draw("");

   TCanvas *c2 = new TCanvas("c2","c2 plot_gain_factors",200,10,700,700);
   c2->Divide(2,2);
   // c2->SetGridx();
   // c2->SetGridy();
   // c2->SetLogz();

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c2->cd(1);
   runnum = 3079;

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3079_r20_60->SetTitle(string);
   h_yvsx3079_r20_60->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3079_r20_60->GetYaxis()->SetTitle("y");
   h_yvsx3079_r20_60->GetXaxis()->SetTitle("x");
   h_yvsx3079_r20_60->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c2->cd(2);

   h_unordered3079_r20_60->SetTitle(string);;
   h_unordered3079_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3079_r20_60->GetYaxis()->SetTitle("Counts");
   h_unordered3079_r20_60->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3079_r20_60->Draw("");

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c2->cd(3);

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3079_r20->SetTitle(string);
   h_yvsx3079_r20->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3079_r20->GetYaxis()->SetTitle("y");
   h_yvsx3079_r20->GetXaxis()->SetTitle("x");
   h_yvsx3079_r20->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c2->cd(4);

   h_unordered3079_r20->SetTitle(string);;
   h_unordered3079_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3079_r20->GetYaxis()->SetTitle("Counts");
   h_unordered3079_r20->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3079_r20->Draw("");   

   TCanvas *c3 = new TCanvas("c3","c3 plot_gain_factors",200,10,700,700);
   c3->Divide(2,2);
   // c3->SetGridx();
   // c3->SetGridy();
   // c3->SetLogz();

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c3->cd(1);
   runnum = 3179;

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3179_r20_60->SetTitle(string);
   h_yvsx3179_r20_60->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3179_r20_60->GetYaxis()->SetTitle("y");
   h_yvsx3179_r20_60->GetXaxis()->SetTitle("x");
   h_yvsx3179_r20_60->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c3->cd(2);

   h_unordered3179_r20_60->SetTitle(string);;
   h_unordered3179_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3179_r20_60->GetYaxis()->SetTitle("Counts");
   h_unordered3179_r20_60->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3179_r20_60->Draw("");

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c3->cd(3);

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3179_r20->SetTitle(string);
   h_yvsx3179_r20->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3179_r20->GetYaxis()->SetTitle("y");
   h_yvsx3179_r20->GetXaxis()->SetTitle("x");
   h_yvsx3179_r20->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c3->cd(4);

   h_unordered3179_r20->SetTitle(string);;
   h_unordered3179_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3179_r20->GetYaxis()->SetTitle("Counts");
   h_unordered3179_r20->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3179_r20->Draw("");


   TCanvas *c4 = new TCanvas("c4","c4 plot_gain_factors",200,10,700,700);
   c4->Divide(2,2);
   // c4->SetGridx();
   // c4->SetGridy();
   // c4->SetLogz();

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c4->cd(1);
   runnum = 3180;

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3180_r20_60->SetTitle(string);
   h_yvsx3180_r20_60->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3180_r20_60->GetYaxis()->SetTitle("y");
   h_yvsx3180_r20_60->GetXaxis()->SetTitle("x");
   h_yvsx3180_r20_60->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c4->cd(2);

   h_unordered3180_r20_60->SetTitle(string);;
   h_unordered3180_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3180_r20_60->GetYaxis()->SetTitle("Counts");
   h_unordered3180_r20_60->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3180_r20_60->Draw("");

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c4->cd(3);

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3180_r20->SetTitle(string);
   h_yvsx3180_r20->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3180_r20->GetYaxis()->SetTitle("y");
   h_yvsx3180_r20->GetXaxis()->SetTitle("x");
   h_yvsx3180_r20->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c4->cd(4);

   h_unordered3180_r20->SetTitle(string);;
   h_unordered3180_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3180_r20->GetYaxis()->SetTitle("Counts");
   h_unordered3180_r20->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3180_r20->Draw("");


   TCanvas *c5 = new TCanvas("c5","c5 plot_gain_factors",200,10,700,700);
   c5->Divide(2,2);
   // c5->SetGridx();
   // c5->SetGridy();
   // c5->SetLogz();

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c5->cd(1);
   runnum = 3182;

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3182_r20_60->SetTitle(string);
   h_yvsx3182_r20_60->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3182_r20_60->GetYaxis()->SetTitle("y");
   h_yvsx3182_r20_60->GetXaxis()->SetTitle("x");
   h_yvsx3182_r20_60->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c5->cd(2);

   h_unordered3182_r20_60->SetTitle(string);;
   h_unordered3182_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3182_r20_60->GetYaxis()->SetTitle("Counts");
   h_unordered3182_r20_60->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3182_r20_60->Draw("");

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c5->cd(3);

   sprintf(string,"Electron Gain Factors Run%4d\n",runnum);
   h_yvsx3182_r20->SetTitle(string);
   h_yvsx3182_r20->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3182_r20->GetYaxis()->SetTitle("y");
   h_yvsx3182_r20->GetXaxis()->SetTitle("x");
   h_yvsx3182_r20->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c5->cd(4);

   h_unordered3182_r20->SetTitle(string);;
   h_unordered3182_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3182_r20->GetYaxis()->SetTitle("Counts");
   h_unordered3182_r20->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3182_r20->Draw("");

   TCanvas *c6 = new TCanvas("c6","c6 plot_gain_factors",200,10,700,700);
   c6->Divide(2,2);
   // c6->SetGridx();
   // c6->SetGridy();
   // c6->SetLogz();

   c6->cd(1);

   h3079_2931->SetTitle("");
   // h3079_2931->GetYaxis()->SetRangeUser(ymin,ymax);
   h3079_2931->GetYaxis()->SetTitle("Gain R3079");
   h3079_2931->GetXaxis()->SetTitle("Gain R2931");
   h3079_2931->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   c6->cd(2);

   h3179_2931->SetTitle("");
   // h3179_2931->GetYaxis()->SetRangeUser(ymin,ymax);
   h3179_2931->GetYaxis()->SetTitle("Gain R3179");
   h3179_2931->GetXaxis()->SetTitle("Gain R2931");
   h3179_2931->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   c6->cd(3);

   h3180_2931->SetTitle("");
   // h3180_2931->GetYaxis()->SetRangeUser(ymin,ymax);
   h3180_2931->GetYaxis()->SetTitle("Gain R3180");
   h3180_2931->GetXaxis()->SetTitle("Gain R2931");
   h3180_2931->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   c6->cd(4);

   h3182_2931->SetTitle("");
   // h3182_2931->GetYaxis()->SetRangeUser(ymin,ymax);
   h3182_2931->GetYaxis()->SetTitle("Gain R3182");
   h3182_2931->GetXaxis()->SetTitle("Gain R2931");
   h3182_2931->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();


   TCanvas *c7 = new TCanvas("c7","c7 plot_gain_factors",200,10,700,700);
   c7->Divide(2,2);
   // c7->SetGridx();
   // c7->SetGridy();
   // c7->SetLogz();

   c7->cd(1);

   h3079_2931_r20_60->SetTitle("");
   // h3079_2931_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h3079_2931_r20_60->GetYaxis()->SetTitle("Gain R3079");
   h3079_2931_r20_60->GetXaxis()->SetTitle("Gain R2931");
   h3079_2931_r20_60->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   c7->cd(2);

   h3179_2931_r20_60->SetTitle("");
   // h3179_2931_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h3179_2931_r20_60->GetYaxis()->SetTitle("Gain R3179");
   h3179_2931_r20_60->GetXaxis()->SetTitle("Gain R2931");
   h3179_2931_r20_60->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   c7->cd(3);

   h3180_2931_r20_60->SetTitle("");
   // h3180_2931_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h3180_2931_r20_60->GetYaxis()->SetTitle("Gain R3180");
   h3180_2931_r20_60->GetXaxis()->SetTitle("Gain R2931");
   h3180_2931_r20_60->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   c7->cd(4);

   h3182_2931_r20_60->SetTitle("");
   // h3182_2931_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h3182_2931_r20_60->GetYaxis()->SetTitle("Gain R3182");
   h3182_2931_r20_60->GetXaxis()->SetTitle("Gain R2931");
   h3182_2931_r20_60->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();


   TCanvas *c8 = new TCanvas("c8","c8 plot_gain_factors",200,10,700,700);
   c8->Divide(2,2);
   // c8->SetGridx();
   // c8->SetGridy();
   // c8->SetLogz();

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c8->cd(1);
   runnum = 3182;

   sprintf(string,"Electron Gain Factors Run 3079-3182");
   h_yvsx3079_3182_r20_60->SetTitle(string);
   h_yvsx3079_3182_r20_60->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3079_3182_r20_60->GetYaxis()->SetTitle("y");
   h_yvsx3079_3182_r20_60->GetXaxis()->SetTitle("x");
   h_yvsx3079_3182_r20_60->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c8->cd(2);

   h_unordered3079_3182_r20_60->SetTitle(string);;
   h_unordered3079_3182_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3079_3182_r20_60->GetYaxis()->SetTitle("Counts");
   h_unordered3079_3182_r20_60->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3079_3182_r20_60->Fit("gaus");
   h_unordered3079_3182_r20_60->Draw("");

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c8->cd(3);

   sprintf(string,"Electron Gain Factors Run 3079-3182");
   h_yvsx3079_3182_r20->SetTitle(string);
   h_yvsx3079_3182_r20->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx3079_3182_r20->GetYaxis()->SetTitle("y");
   h_yvsx3079_3182_r20->GetXaxis()->SetTitle("x");
   h_yvsx3079_3182_r20->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c8->cd(4);

   h_unordered3079_3182_r20->SetTitle(string);;
   h_unordered3079_3182_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered3079_3182_r20->GetYaxis()->SetTitle("Counts");
   h_unordered3079_3182_r20->GetXaxis()->SetTitle("Gain Factor");
   h_unordered3079_3182_r20->Draw("");


   TCanvas *c9 = new TCanvas("c9","c9 plot_gain_factors",200,10,700,700);
   c9->Divide(2,2);
   // c9->SetGridx();
   // c9->SetGridy();
   // c9->SetLogz();

   c9->cd(1);

   h3079_3182_2931_r20_60->SetTitle("");
   // h3079_3182_2931_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   h3079_3182_2931_r20_60->GetYaxis()->SetTitle("Gain R3079-3182");
   h3079_3182_2931_r20_60->GetXaxis()->SetTitle("Gain R2931");
   h3079_3182_2931_r20_60->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   c9->cd(2);

   h3079_3182_2931_r20->SetTitle("");
   // h3079_3182_2931_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h3079_3182_2931_r20->GetYaxis()->SetTitle("Gain R3079-3182");
   h3079_3182_2931_r20->GetXaxis()->SetTitle("Gain R2931");
   h3079_3182_2931_r20->Draw("colz");
   TLine *diag = new TLine (0,0,2,2);
   diag->SetLineColor(2);
   diag->SetLineWidth(2);
   diag->Draw();

   Float_t ymin=0;
   Float_t ymax=150;
   c9->cd(3);

   gRatio_r20_60->SetTitle(string);;
   // gRatio_r20_60->GetYaxis()->SetRangeUser(ymin,ymax);
   gRatio_r20_60->GetYaxis()->SetTitle("Counts");
   gRatio_r20_60->GetXaxis()->SetTitle("Ratio R3079-3182/R2931 _r20_60");
   gRatio_r20_60->Fit("gaus");
   gRatio_r20_60->Draw("");

   Float_t ymin=0;
   Float_t ymax=150;
   c9->cd(4);

   gRatio_r20->SetTitle(string);;
   // gRatio_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   gRatio_r20->GetYaxis()->SetTitle("Counts");
   gRatio_r20->GetXaxis()->SetTitle("Ratio R3079-3182/R2931 _r20");
   gRatio_r20->Fit("gaus");
   gRatio_r20->Draw("");


   TCanvas *c10 = new TCanvas("c10","c10 plot_gain_factors",200,10,700,700);
   c10->Divide(2,2);
   // c10->SetGridx();
   // c10->SetGridy();
   // c10->SetLogz();

   Float_t zmin=0.8;
   Float_t zmax=1.6;
   c10->cd(1);
   runnum = 3182;

   sprintf(string,"Electron Gain Factors Run 2931-3182");
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

   sprintf(string,"Electron Gain Factors Run 2931-3182");
   h_yvsx2931_3182_r20->SetTitle(string);
   h_yvsx2931_3182_r20->GetZaxis()->SetRangeUser(zmin,zmax);
   h_yvsx2931_3182_r20->GetYaxis()->SetTitle("y");
   h_yvsx2931_3182_r20->GetXaxis()->SetTitle("x");
   h_yvsx2931_3182_r20->Draw("colz");

   Float_t ymin=0;
   Float_t ymax=150;
   c10->cd(4);

   h_unordered2931_3182_r20->SetTitle(string);;
   h_unordered2931_3182_r20->GetYaxis()->SetRangeUser(ymin,ymax);
   h_unordered2931_3182_r20->GetYaxis()->SetTitle("Counts");
   h_unordered2931_3182_r20->GetXaxis()->SetTitle("Gain Factor");
   h_unordered2931_3182_r20->Draw("");

     /*sprintf (string,"#sigma=%.2f\n",sigma);
     printf("string=%s",string);
     t1 = new TLatex(0.4,0.5,string);
     t1->SetNDC();
     t1->SetTextSize(0.07);
     t1->Draw();*/

   // print histograms

   sprintf(string,"plot_gain_factors.pdf(");
   c1->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c2->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c3->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c4->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c5->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c6->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c7->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c8->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf");
   c10->SaveAs(string);
   sprintf(string,"plot_gain_factors.pdf)");
   c9->SaveAs(string);


   // close files

   outgain.close();

   // in->Close();


}

