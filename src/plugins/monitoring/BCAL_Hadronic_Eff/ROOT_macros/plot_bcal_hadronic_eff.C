void plot_bcal_hadronic_eff(void)
{
// read output of Read_bcal_hadronic_eff2 dat files with efficiency information and plot vs run number.
//

#include <TRandom.h>

gROOT->Reset();
//TTree *Bfield = (TTree *) gROOT->FindObject("Bfield");
gStyle->SetPalette(1,0);
gStyle->SetOptStat(kFALSE);
// gStyle->SetOptStat(11111111);
gStyle->SetPadRightMargin(0.15);
gStyle->SetPadLeftMargin(0.15);
gStyle->SetPadBottomMargin(0.15);
//

   char string[256];
    const Int_t nruns=6;
    
    /*vector <float> eff_up1(nruns,0.);
    vector <float> eff_up2(nruns,0.);
    vector <float> eff_up3(nruns,0.);
    vector <float> eff_up4(nruns,0.);
    
    vector <float> eff_down1(nruns,0.);
    vector <float> eff_down2(nruns,0.);
    vector <float> eff_down3(nruns,0.);
    vector <float> eff_down4(nruns,0.);*/
    
    Float_t eff_up1[nruns]={0,0,0,0,0,0};
    Float_t eff_up2[nruns]={0,0,0,0,0,0};
    Float_t eff_up3[nruns]={0,0,0,0,0,0};
    Float_t eff_up4[nruns]={0,0,0,0,0,0};
    
    Float_t eff_down1[nruns]={0,0,0,0,0,0};
    Float_t eff_down2[nruns]={0,0,0,0,0,0};
    Float_t eff_down3[nruns]={0,0,0,0,0,0};
    Float_t eff_down4[nruns]={0,0,0,0,0,0};
    
    Float_t runnum[nruns]={10492,11529,11658,30734,30801};
    
    map<Int_t, TString> runs;
    runs[0]="010492";
    runs[1]="011529";
    runs[2]="011659";
    runs[3]="030734";
    runs[4]="030801";
    runs[5]="030890";
    
    Int_t layer=1;
    Int_t coinc_cut=3;

    TString datfile;
    
    for (Int_t jrun=0; jrun<nruns; jrun++) {
    
    
    for (layer=1; layer<5; layer++) {
    
    datfile = "dat/R"+runs[jrun]+"_layer"+TString::Itoa(layer,10)+"_cut"+TString::Itoa(coinc_cut,10)+".dat";
    
    cout << "Opening file: " << datfile.Data() << endl;

   ifstream file;
   file.open (datfile.Data());
   if (!file) {
       cout << "ERROR: Failed to open data file= " << datfile.Data() << endl;
       return;
       }
        
    Int_t ndx=0;
    TString line;
   while (line.ReadLine(file)){
       
       cout << "line=" << line << endl;
       
       TObjArray *tokens = line.Tokenize(" ");
       Int_t ntokens = tokens->GetEntries();
       // Float_t  = (((TObjString*)tokens->At(1))->GetString()).Atof();
       // cout << " token1=" << ((TObjString*)tokens->At(0))->GetString()  << " token2=" << ((TObjString*)tokens->At(1))->GetString() << endl;
       if (ndx == 0) {
           Float_t run = ((TObjString*)tokens->At(1))->GetString().Atof();
           Float_t xlayer = ((TObjString*)tokens->At(3))->GetString().Atof();
           Float_t xcoinc_cut = ((TObjString*)tokens->At(5))->GetString().Atof();
           cout << " run=" << run << " xlayer=" << xlayer << " xcoinc_cut=" << xcoinc_cut << endl;
           runnum[jrun] = run;
       }
       if (ndx == 1) {
           Float_t eff_up = ((TObjString*)tokens->At(7))->GetString().Atof();
           Float_t eff_down = ((TObjString*)tokens->At(9))->GetString().Atof();
           cout << " eff_up=" << eff_up << " eff_down=" << eff_down << endl;
           
           switch (layer) {
               case 1:
                   eff_up1[jrun] = eff_up;
                   eff_down1[jrun] = eff_down;
                   break;
               case 2:
                   eff_up2[jrun] = eff_up;
               		eff_down2[jrun] = eff_down;
                   break;
               case 3:
                   eff_up3[jrun] = eff_up;
                   eff_down3[jrun] = eff_down;
                   break;
               case 4:
                   eff_up4[jrun] = eff_up;
                   eff_down4[jrun] = eff_down;
                   break;
               default:
                   cout << "*** Illegal layer=" << layer << endl;
                   return;
           }  // end switch
       } // end if
       
       ndx++;
       
   }   // end reading over dat file

   file.close();
        
    }
    }
    
    cout << "up1=" << eff_up1[0] << endl;
    cout << "up2=" << eff_up1[1] << endl;
    cout << "up3=" << eff_up1[2] << endl;
    
    TGraph *gr_eff_up1 = new TGraph(nruns,runnum,eff_up1);
    TGraph *gr_eff_up2 = new TGraph(nruns,runnum,eff_up2);
    TGraph *gr_eff_up3 = new TGraph(nruns,runnum,eff_up3);
    TGraph *gr_eff_up4 = new TGraph(nruns,runnum,eff_up4);
    
    TGraph *gr_eff_down1 = new TGraph(nruns,runnum,eff_down1);
    TGraph *gr_eff_down2 = new TGraph(nruns,runnum,eff_down2);
    TGraph *gr_eff_down3 = new TGraph(nruns,runnum,eff_down3);
    TGraph *gr_eff_down4 = new TGraph(nruns,runnum,eff_down4);
    
    
    TCanvas *c0 = new TCanvas("c0", "c0",200,10,1000,700);
       
       
    gr_eff_up1->SetTitle("");
    // gr_eff_up1->GetXaxis()->SetRangeUser(xmin,xmax);
    gr_eff_up1->GetYaxis()->SetRangeUser(0.5,1);
    gr_eff_up1->GetXaxis()->SetTitleSize(0.05);
    gr_eff_up1->GetYaxis()->SetTitleSize(0.05);
    gr_eff_up1->GetYaxis()->SetTitle("Efficiency");
    gr_eff_up1->GetXaxis()->SetTitle("Run Number");
    gr_eff_up1->SetMarkerColor(4);
    gr_eff_up1->SetMarkerStyle(20);
    gr_eff_up1->Draw("Ap");
    
    gr_eff_up2->SetMarkerColor(2);
    gr_eff_up2->SetMarkerStyle(20);
    gr_eff_up2->Draw("psame");
    
    gr_eff_up3->SetMarkerColor(1);
    gr_eff_up3->SetMarkerStyle(20);
    gr_eff_up3->Draw("psame");
    
    gr_eff_up4->SetMarkerColor(3);
    gr_eff_up4->SetMarkerStyle(20);
    gr_eff_up4->Draw("psame");
    
    TLegend *leg = new TLegend(0.6,0.4,0.8,0.6);
    leg->AddEntry(gr_eff_up1,"Layer 1","p");
    leg->AddEntry(gr_eff_up2,"Layer 2","p");
    leg->AddEntry(gr_eff_up3,"Layer 3","p");
    leg->AddEntry(gr_eff_up4,"Layer 4","p");
    leg->Draw();
    
    
    c0->SaveAs("plot_bcal_hadronic_eff.pdf");

}

