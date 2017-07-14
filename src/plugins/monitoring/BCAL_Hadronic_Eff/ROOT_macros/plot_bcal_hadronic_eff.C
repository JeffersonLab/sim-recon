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
    const Int_t nruns=58;
    
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
    
    Float_t runnum[nruns]={10492,10498,10590,10598,10659,10778,10873,10907,11082,11128,11157,11264,11300,11405,11436,11529,11659,
        11366,11384,11407,11432,11436,11450,11473,11483,11510,11520,11553,
        30734,30801,30890,30830,30833,30834,30835,30836,30838,30900,30902,30903,30926,30927,30928,30930,
        30274,30280.30300,30350,30405,30450,30497,30592,30610,30650,30961,31000,31029,31056};
    
    map<Int_t, TString> runs;
    runs[0]="010492";
    runs[1]="010498";
    runs[2]="010590";
    runs[3]="010598";
    runs[4]="010659";
    runs[5]="010778";
    runs[6]="010873";
    runs[7]="010907";
    runs[8]="011082";
    runs[9]="011128";
    runs[10]="011157";
    runs[11]="011264";
    runs[12]="011300";
    runs[13]="011405";
    runs[14]="011436";
    runs[15]="011529";
    runs[16]="011659";
    runs[17]="011366";
    runs[18]="011384";
    runs[19]="011407";
    runs[20]="011432";
    runs[21]="011436";
    runs[22]="011450";
    runs[23]="011473";
    runs[24]="011483";
    runs[25]="011510";
    runs[26]="011520";
    runs[27]="011553";
    runs[28]="030734";
    runs[29]="030801";
    runs[30]="030890";
    runs[31]="030830";
    runs[32]="030833";
    runs[33]="030834";
    runs[34]="030835";
    runs[35]="030836";
    runs[36]="030838";
    runs[37]="030900";
    runs[38]="030902";
    runs[39]="030903";
    runs[40]="030926";
    runs[41]="030927";
    runs[42]="030928";
    runs[43]="030930";
    runs[44]="030274";
    runs[45]="030280";
    runs[46]="030300";
    runs[47]="030350";
    runs[48]="030405";
    runs[49]="030450";
    runs[50]="030497";
    runs[51]="030592";
    runs[52]="030610";
    runs[53]="030650";
    runs[54]="030961";
    runs[55]="031000";
    runs[56]="031029";
    runs[57]="031056";
    
    
    
    
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
    
    TLegend *leg = new TLegend(0.6,0.3,0.8,0.5);
    leg->AddEntry(gr_eff_up1,"Layer 1","p");
    leg->AddEntry(gr_eff_up2,"Layer 2","p");
    leg->AddEntry(gr_eff_up3,"Layer 3","p");
    leg->AddEntry(gr_eff_up4,"Layer 4","p");
    leg->Draw();
    
    
    TCanvas *c1 = new TCanvas("c1", "c1",200,10,1000,700);
    
    TGraph *gr_eff_up1_copy = (TGraph*)gr_eff_up1->Clone("gr_eff_up1_copy");
    
    gr_eff_up1_copy->GetYaxis()->SetRangeUser(0.9,1.0);
    gr_eff_up1_copy->GetXaxis()->SetRangeUser(30200,31050);
    
    gr_eff_up1_copy->Draw("Ap");
    gr_eff_up2->Draw("psame");
    gr_eff_up3->Draw("psame");
    gr_eff_up4->Draw("psame");
    leg->Draw();
    
    // 2017 ranges
    
    xmin = 10000;
    xmax = 12000;
    ymin = 0.9;
    ymax = 1.0;
    
    TLine *linea = new TLine(30795,ymin,30795,ymax);
    linea->Draw();
    
    TLatex *t1 = new TLatex(30274,1.001,"2017      Low Intensity");    // t1->SetNDC();
    t1->SetTextSize(0.03);
    t1->Draw();
    t1->DrawLatex(30800,1.001,"Hi Intensity");
    
    TCanvas *c2 = new TCanvas("c2", "c2",200,10,1000,700);
    
    TGraph *gr_eff_up1_copy2 = (TGraph*)gr_eff_up1->Clone("gr_eff_up1_copy");
    gr_eff_up1_copy2->GetYaxis()->SetRangeUser(ymin,ymax);
    gr_eff_up1_copy2->GetXaxis()->SetRangeUser(xmin,xmax);
    
    gr_eff_up1_copy2->Draw("Ap");
    gr_eff_up2->Draw("psame");
    gr_eff_up3->Draw("psame");
    gr_eff_up4->Draw("psame");
    leg->Draw();
    
    // 2016 ranges
    
    TLine *line1 = new TLine(11059,ymin,11059,ymax);
    line1->Draw();
    TLine *line2 = new TLine(11366,ymin,11366,ymax);
    line2->Draw();
    TLine *line3 = new TLine(11555,ymin,11555,ymax);
    line3->Draw();

    t1->DrawLatex(11059,1.001,"CDC Hole");    // t1->SetNDC();
    t1->SetTextSize(0.03);
    t1->Draw();
    t1->DrawLatex(11366,1.001,"Golden");
    t1->DrawLatex(11655,1.001,"low rate");
    t1->DrawLatex(10495,1.001,"2016");
    
    c0->SaveAs("plot_bcal_hadronic_eff.pdf(");
    c1->SaveAs("plot_bcal_hadronic_eff.pdf");
    c2->SaveAs("plot_bcal_hadronic_eff.pdf)");
    
}

