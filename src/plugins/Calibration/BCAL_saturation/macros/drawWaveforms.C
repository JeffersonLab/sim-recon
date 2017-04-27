#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TAxis.h"

double par0 = 55800;
double par1 = -5.73e-6;
double par2 = -1.03e-10;
double par0err, par1err, par2err;

Double_t myfunction(Double_t *x, Double_t *par)
{
   Float_t xx =x[0];
   Double_t f = 1;
   if(xx > par[0]) 
     f += (par[1]*(xx-par[0]) + par[2]*pow((xx - par[0]),2));

   return f;
}

void drawWaveform(TString bcalEnd = "DS", int layer = 1, int run = 30894)
{
    double lineWidth = 2.0;
    double textSize = 0.07;
    
    bool fall2016 = false;
    
    TString figsDir = "./figs/"; figsDir+=run; figsDir+="/";
    
    TH2D* hIntegralRatio_UnsaturatedIntegral = new TH2D("IntegralRatio_UnsaturatedIntegral", "; Unsaturated Integral (fADC); Saturated/Unsaturated Integral", 750, 0, 120000, 100, 0.5, 1.05);
    TH2D* hIntegralRatio_SaturatedIntegral = new TH2D("IntegralRatio_SaturatedIntegral", "; Saturated Integral (fADC); Saturated/Unsaturated Integral", 750, 0, 120000, 100, 0.5, 1.05);
    
    TH2D* hIntegralResidual_SaturatedIntegral = new TH2D("IntegralResidual_SaturatedIntegral", "; Saturated Integral (fADC); Integral Residual", 750, 0, 120000, 100, -0.4, 0.4);
    
    TF1 *fknown = new TF1("fKnown",myfunction,0,160000,3);
    fknown->SetParameters(par0, par1, par2);

    TH2D* hIntegralRatio_UnsaturatedIntegralFull = new TH2D("IntegralRatio_UnsaturatedIntegralFull", "; Unsaturated Full Integral (fADC); Saturated/Unsaturated Full Integral", 500, 0, 200000, 100, 0.5, 1.05);
    
    TString fileName = "hd_root_";
    fileName+=run; fileName+="_saturate.root";
    TFile *f = TFile::Open(fileName);
    
    TH2I* hbcal_waveform_saturate_2D = (TH2I*)f->Get(Form("bcal_saturation/bcal%s_waveform_saturate_layer%d", bcalEnd.Data(), layer));
    TH2I* hbcal_waveform_2D = (TH2I*)f->Get(Form("bcal_saturation/bcal%s_waveform_noSaturate_layer%d", bcalEnd.Data(), layer));
    //h->GetYaxis()->SetNdivisions(4);
    //h->GetXaxis()->SetNdivisions(6);
    //h->SetLineWidth(lineWidth);
    
    TCanvas* aa = new TCanvas("aa", "aa", 600, 400);
    hbcal_waveform_saturate_2D->Draw("colz");
    hbcal_waveform_saturate_2D->SetMaximum(4096);
    hbcal_waveform_saturate_2D->GetXaxis()->SetRangeUser(0,400);
    aa->Print(Form("%sSaturate_%s_layer%d.pdf",figsDir.Data(),bcalEnd.Data(),layer));
    
    TCanvas* bb = new TCanvas("bb", "bb", 600, 400);
    hbcal_waveform_2D->Draw("colz");
    
    TCanvas* cc = new TCanvas("cc", "cc", 600, 400);
    TH1F* hBase = new TH1F("hBase", ";Samples; fADC", 100, 0, 100);
    hBase->SetMaximum(6000.);
    hBase->Draw();
    
    int nbins = hbcal_waveform_2D->GetXaxis()->GetNbins();
    for(int i=0; i<nbins; i++) {
        
        TH1I* hbcal_waveform = (TH1I*)hbcal_waveform_2D->ProjectionY(Form("Event_%d",i),i+1,i+1);
        if(hbcal_waveform->Integral() < 10) continue;
        
        hbcal_waveform->SetLineColor(kBlack);
        hbcal_waveform->SetLineWidth(2);
        int thresh = 120;
        int threshCrossBin = 0;
        for(int j=0; j<hbcal_waveform->GetXaxis()->GetNbins(); j++) {
            if(hbcal_waveform->GetBinContent(j) > thresh){
                threshCrossBin = j;
                break;
            }
        }
        
        TH1I* hbcal_waveform_temp = (TH1I*)hbcal_waveform->Clone();
        hbcal_waveform_temp->Reset();
        
        // put all waveforms on the same x-axis
        int nSaturate = 0;
        int sampleOffset = 20;
        int binMax = 0;
        double fADCmax = 0;
        int nSamples = 0;
        for(int j=0; j<hbcal_waveform->GetXaxis()->GetNbins(); j++) {
            double fADC = hbcal_waveform->GetBinContent(j+1);
            if(fADC > 4096) {
                fADC = 4096;
                nSaturate++;
            }
            hbcal_waveform_temp->SetBinContent(j+1 - threshCrossBin+sampleOffset, fADC);
            hbcal_waveform_temp->SetBinError(j+1 - threshCrossBin+sampleOffset, 0.); //sqrt(fADC));
            
            // find maximum bin for outlier waveform shapes
            if(fADC > fADCmax) {
                binMax = j+1 - threshCrossBin;
                fADCmax = fADC;
            }
            
            // find how many samples in the readout window
            if(fADC > thresh)
                nSamples++;
        }
        
        // remove waveforms with peaks too far after threshold crossing
        if(binMax > 5)
            continue;
        
        // remove waveforms with too few samples in readout window
        if(nSamples < 26)
            continue;
        
        cc->cd();
        hbcal_waveform_temp->Draw("h same");
        
        int pedestal = 100;
        int integralNSA = 26; // Spring 2016 = 26 and Fall 2016 = 13
        int integralNSB = 1;  // Spring 2016 = 1  and Fall 2016 = -2
        if(fall2016) {
            integralNSA = 13;
            integralNSB = -2;
        }
        
        int integralBinMin = sampleOffset - integralNSB;
        int integralBinMax = sampleOffset + integralNSA;
        
        // demonstration with single waveform
        if(i == 211) {
            hbcal_waveform_temp->SetTitle("");
            TH1I *hdemo_scale = (TH1I*)hbcal_waveform_temp->Clone();
            double scale = 1.5;
            hdemo_scale->Scale(scale);
            TH1I *hdemo_saturate = (TH1I*)hdemo_scale->Clone();
            for(int j=0; j<hdemo_saturate->GetXaxis()->GetNbins(); j++) {
                if(j+1 < integralBinMin || j+1 >= integralBinMax)
                    hdemo_saturate->SetBinContent(j+1, 0);
                if(hdemo_saturate->GetBinContent(j+1) > 4095)
                    hdemo_saturate->SetBinContent(j+1, 4096);
            }
            hdemo_saturate->SetFillColor(kRed);
            
            TCanvas *demo = new TCanvas("demo", "demo", 800, 400);
            demo->Divide(2,1);
            demo->cd(1);
            hbcal_waveform_temp->SetMaximum(6000);
            hbcal_waveform_temp->Draw();
            demo->cd(2);
            hdemo_scale->SetMaximum(6000);
            hdemo_scale->Draw();
            hdemo_saturate->Draw("same");
            
            demo->Print(Form("%sexampleSaturate_%0.1f.pdf", figsDir.Data(), scale));
        }
        
        double integral_unscaled = 0;
        double integral_unscaled_full = 0;
        double integral_unsaturated = 0;
        double integral_unsaturated_full = 0;
        
        // try different scale factors for integral
        for(int j=0; j<100; j++) {
            double scale = 1. + j*0.02;
            
            // calculate integral
            double integral = 0;
            double integral_full = 0;
            for(int k=integralBinMin; k<hbcal_waveform_temp->GetXaxis()->GetNbins(); k++) {
                double fADC = scale * (hbcal_waveform_temp->GetBinContent(k) - pedestal);
                double fADC_truncated = fADC > (4095 - pedestal) ? (4095 - pedestal) : fADC;
                
                // only calculate integral over used samples
                if(k < integralBinMax)
                    integral += fADC_truncated;
                
                // test using full integral (maybe relevant for mcsmear?)
                integral_full += fADC_truncated;
            }
            if(j==0) {
                integral_unscaled = integral;
                integral_unscaled_full = integral_full;
            }
            else {
                integral_unsaturated = integral_unscaled*scale;
                integral_unsaturated_full = integral_unscaled_full*scale;
            }
            
            // check for outlier fADC pulses
            if(integral_unsaturated < 60000 && integral/integral_unsaturated < 0.9)
                cout<<i<<endl;
            
            hIntegralRatio_UnsaturatedIntegral->Fill(integral_unsaturated, integral/integral_unsaturated);
            hIntegralRatio_SaturatedIntegral->Fill(integral, integral/integral_unsaturated);
            double integral_residual = (integral/fknown->Eval(integral) - integral_unsaturated)/integral_unsaturated;
            hIntegralResidual_SaturatedIntegral->Fill(integral, integral_residual);
            
            hIntegralRatio_UnsaturatedIntegralFull->Fill(integral_unsaturated_full, integral_full/integral_unsaturated_full);
        }
        
        //cc->Print(Form("figs/bcal_waveform_%03d.ps", i));
    }
    
    double maxXaxis = 1.2e5;
    if(fall2016)
        maxXaxis = 10e4;
    
    TCanvas* dd = new TCanvas("dd", "dd", 1000, 400);
    dd->Divide(3,1);
    dd->cd(1);
    hIntegralRatio_UnsaturatedIntegral->Draw("colz");
    
    TF1 *f1 = new TF1("myfunc",myfunction,0,160000,3);
    f1->SetParameters(55000,-4e-6,0);
    f1->SetLineColor(kRed);
    f1->SetParNames("offset","linear","quad");
    //hIntegralRatio_UnsaturatedIntegral->Fit(f1, "", "", 30000, 150000);
    hIntegralRatio_UnsaturatedIntegral->Draw("colz");
    hIntegralRatio_UnsaturatedIntegral->GetXaxis()->SetRangeUser(0., maxXaxis);
    
    dd->cd(2);
    hIntegralRatio_SaturatedIntegral->Draw("colz");
    hIntegralRatio_SaturatedIntegral->GetXaxis()->SetRangeUser(0., maxXaxis);
    
    TF1 *f2 = new TF1("myfunc",myfunction,0,160000,3);
    f2->SetParameters(55000,-4e-6,0);
    if(fall2016) {
        f2->SetParameters(35000,0,-4e-6);
        //f2->FixParameter(0, 34000); // Fall 2016
        //f2->FixParameter(1, 0); // Fall 2016
    }
    f2->SetLineColor(kRed);
    f2->SetParNames("offset","linear","quad");
    if(fall2016)
        hIntegralRatio_SaturatedIntegral->Fit(f2, "", "", 30000, 43000);  // Fall 2016
    else
        hIntegralRatio_SaturatedIntegral->Fit(f2, "", "", 30000, 150000); // Spring 2016
    
    // set parameters for plotting residual
    par0 = f2->GetParameter(0); par0err = f2->GetParError(0);
    par1 = f2->GetParameter(1); par1err = f2->GetParError(1);
    par2 = f2->GetParameter(2); par2err = f2->GetParError(2);
    
    hIntegralRatio_SaturatedIntegral->Draw("colz");
    
    dd->cd(3);
    hIntegralResidual_SaturatedIntegral->Draw("colz");
    hIntegralResidual_SaturatedIntegral->GetXaxis()->SetRangeUser(0., maxXaxis);
    
    dd->Print(Form("%s%s_layer%d.pdf",figsDir.Data(),bcalEnd.Data(),layer));
    
    return;
    
}

// main function to loop over BCAL layers and ends
void drawWaveforms(int run = 30894) {
    
    gSystem->Exec(Form("mkdir -p figs/%d", run));
    
    ofstream parameters;
    string parFileName = Form("./figs/%d/parameters.txt", run);
    //cout<<parFileName.data()<<endl;
    parameters.open(parFileName);
    
    TH1F *h0 = new TH1F("par0", "Integral Offset Parameter", 8, 0, 8);
    TH1F *h1 = new TH1F("par1", "Linear Parameter", 8, 0, 8);
    TH1F *h2 = new TH1F("par2", "Quadratic Parameter", 8, 0, 8);
    
    TString bcalEnd[2] = {"US", "DS"};
    for(int end=0; end<2; end++) {
        for(int layer=1; layer<=4; layer++) {

            drawWaveform(bcalEnd[end], layer, run); // first fit to determine parameters
            drawWaveform(bcalEnd[end], layer, run); // then fill residual histogram with fitted parameters
            
            parameters<<end<<"  "<<layer<<"  "<<par0<<"  "<<par1<<"  "<<par2<<endl;
            h0->SetBinContent(end*4 + layer, par0);
            h0->SetBinError(end*4 + layer, par0err);
            h1->SetBinContent(end*4 + layer, par1);
            h1->SetBinError(end*4 + layer, par1err);
            h2->SetBinContent(end*4 + layer, par2);
            h2->SetBinError(end*4 + layer, par2err);
        }
    }
    
    TCanvas *pars = new TCanvas("pars", "pars", 1200, 500);
    pars->Divide(3,1);
    pars->cd(1);
    h0->Draw();
    h0->Fit("pol0");
    pars->cd(2);
    h1->Draw();
    h1->SetMinimum(-8e-6); h1->SetMaximum(-3e-6);
    h1->Fit("pol0");
    pars->cd(3);
    h2->Draw();
    h2->SetMinimum(-1.6e-10); h2->SetMaximum(-1.0e-11);
    h2->Fit("pol0");
    TString parFigName = "./figs/";
    parFigName+=run; parFigName+="/pars.pdf";
    pars->Print(parFigName);
    
    parameters.close();
    
    return;
}

