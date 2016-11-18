// Script to extract time-walk constants for the BCAL

#include "TString.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TSystem.h"
#include "TF1.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TStyle.h"
#include "TString.h"
// #include "T.h"

#include <iostream>
#include <fstream>
#include <sstream>

namespace ExtractTimeOffsetsAndCEffNS {
    //Leave this global so the accesors don't need the pointer as an argument
    TFile *thisFile;

    // Accessor functions to grab histograms from our file
    // (Makes things easy with the HistogramTools fills)
    TH1I * Get1DHistogram(const char * plugin, const char * directoryName, const char * name, bool print = true){
        TH1I * histogram;
        TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
        thisFile->GetObject(fullName, histogram);
        if (histogram == 0){
            if (print) cout << "Unable to find histogram " << fullName.Data() << endl;
            return NULL;
        }
        return histogram;
    }

    TH1I * Get1DHistogramBasedir(const char * name, bool print = true){
        TH1I * histogram;
        TString fullName = TString(name);
        thisFile->GetObject(fullName, histogram);
        if (histogram == 0){
            if (print) cout << "Unable to find histogram " << fullName.Data() << endl;
            return NULL;
        }
        return histogram;
    }

    TH2I * Get2DHistogram(const char * plugin, const char * directoryName, const char * name){
        TH2I * histogram;
        TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
        thisFile->GetObject(fullName, histogram);
        if (histogram == 0){
            cout << "Unable to find histogram " << fullName.Data() << endl;
            return NULL;
        }
        return histogram;
    }
}

TH1D* GetGausFitMean(TH2I *hist2d, TString tag, bool makeplots);



// Do the extraction of the actual constants
void ExtractTimeOffsets(int run = 2931, TString filename = "hd_root.root", TString tag="default"){
    printf("ExtractTimeOffsets\n");

    // Open our input and output file
    ExtractTimeOffsetsAndCEffNS::thisFile = TFile::Open(filename);
    char outfilename[255];    
    sprintf(outfilename,"BCALTimeOffsetAndCEff_Results_%i_%s.root",run,tag.Data());
    TFile *outputFile = TFile::Open(outfilename, "RECREATE");
    outputFile->mkdir("Fits");
    outputFile->mkdir("ResultOverview");

    // Check to make sure it is open
    if (ExtractTimeOffsetsAndCEffNS::thisFile == 0) {
        cout << "Unable to open file " << filename.Data() << "...Exiting" << endl;
        return;
    }
    
    // This stream will be for outputting the results in a format suitable for the CCDB
    ofstream channel_global_offset_File;
    sprintf(outfilename,"channel_global_offset_BCAL_%i_%s.txt",run,tag.Data());
    channel_global_offset_File.open(outfilename);

    TH1I *OldTimeOffsets = new TH1I("OldTimeOffsets", "Total time offset (pre)", 2000, -4, 4);
    TH1F *OldTimeOffsetsVsChannel = new TH1F("OldTimeOffsetsVsChannel", "Total time offset (pre) vs channel;CCDB channel", 768, 0.5, 768.5);
    TH1I *NewTimeOffsets = new TH1I("NewTimeOffsets", "Total time offset (post)", 2000, -4, 4);
    TH1F *NewTimeOffsetsVsChannel = new TH1F("NewTimeOffsetsVsChannel", "Total time offset (post) vs channel;CCDB channel", 768, 0.5, 768.5);
    TH1I *NewTimeOffsetsFit = new TH1I("NewTimeOffsetsFit", "Total time offset fit (post)", 2000, -4, 4);
    TH1F *NewTimeOffsetsFitVsChannel = new TH1F("NewTimeOffsetsFitVsChannel", "Total time offset fit (post) vs channel;CCDB channel", 768, 0.5, 768.5);
    TH1I *RelativeTimeOffsetsFit = new TH1I("RelativeTimeOffsetsFit", "Change in time offset", 2000, -4, 4);
    TH1F *RelativeTimeOffsetsFitVsChannel = new TH1F("RelativeTimeOffsetsFitVsChannel", "Change in time offset vs channel;CCDB channel", 768, 0.5, 768.5);

    TH1D *priorOffsetHist = NULL;
    priorOffsetHist = (TH1D*)ExtractTimeOffsetsAndCEffNS::thisFile->Get("BCAL_Global_Offsets/Target Time/CCDB_raw_channel_global_offset");
    auto residualHist = ExtractTimeOffsetsAndCEffNS::Get2DHistogram("BCAL_Global_Offsets", "Target Time", "deltaTVsCell_q-");

    if (priorOffsetHist != NULL) {
        // for (int i = 1 ; i <= priorOffsetHist->GetNbinsX(); i++){
		// 	printf("%4i   %-.3f\n",  i,priorOffsetHist->GetBinContent(i));
		// }
    } else {
        printf("Failed to load %s\n","CCDB_raw_channel_global_offset");
    }

    if(residualHist != NULL && priorOffsetHist != NULL){
        int nBinsX = residualHist->GetNbinsX();
        int nBinsY = residualHist->GetNbinsY();

        TString newtag = tag + "_channel";
        bool makeplots=1; // print images of the fits
        TH1D *meanhist = GetGausFitMean(residualHist,tag.Data(),makeplots);

        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = residualHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 0.5; //ns (Full Width)
            int binWindow = int(timeWindow / nsPerBin);
            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + binWindow) <= projY->GetNbinsX() ? (j + binWindow) : projY->GetNbinsX();
                double sum = 0, nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++){
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin){
                        if (nEntries > maxEntries) {
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    }
                }
            }
            float priorOffset = priorOffsetHist->GetBinContent(i);
            float residualOffset = maxMean;
            float residualOffsetFit = meanhist->GetBinContent(i);
            float newOffset = priorOffset + residualOffset;
            float newOffsetFit = priorOffset + residualOffsetFit;

            OldTimeOffsets->Fill(priorOffset);
            OldTimeOffsetsVsChannel->SetBinContent(i,priorOffset);
            NewTimeOffsets->Fill(newOffset);
            NewTimeOffsetsVsChannel->SetBinContent(i,newOffset);
            NewTimeOffsetsFit->Fill(newOffsetFit);
            NewTimeOffsetsFitVsChannel->SetBinContent(i,newOffsetFit);
            RelativeTimeOffsetsFit->Fill(residualOffsetFit);
            RelativeTimeOffsetsFitVsChannel->SetBinContent(i,residualOffsetFit);

            printf("%4i   rolling max: %6.3f %+.3f = %6.3f    gaus fit: %6.3f %+.3f = %6.3f    diff: %6.3f\n",
                   i,priorOffset,residualOffset,newOffset,priorOffset,residualOffsetFit,newOffsetFit,residualOffset-residualOffsetFit);
            channel_global_offset_File <<  newOffsetFit << endl;;
        }
    }
    channel_global_offset_File.close();
    outputFile->Write();
    ExtractTimeOffsetsAndCEffNS::thisFile->Close();
}


TH1D* GetGausFitMean(TH2I *hist2d, TString tag, bool makeplots=0) {
    // Take a 2D histogram and return the Gaussian mean of each X bin (somewhat like a profile)
    double outermintime=-2, outermaxtime=2;

    TAxis * xaxis = hist2d->GetXaxis();
    TAxis * yaxis = hist2d->GetYaxis();
    int nBinsX = hist2d->GetNbinsX();
    int nBinsY = hist2d->GetNbinsY();
        
    printf("bins (%i,%i)\n",nBinsX,nBinsY);

    TH1D *xmean = hist2d->ProjectionX("xmean");
    xmean->Reset();

    gStyle->SetOptStat(2210);
    gStyle->SetOptFit(1);

    TCanvas *canvas_fit;
    canvas_fit = new TCanvas("canvas_fit","canvas_fit");
    gPad->SetLogy();

    for (int i = 1 ; i <= nBinsX; i++){
        TH1D *singlebin_projY = hist2d->ProjectionY("temp", i, i);
        float center = yaxis->GetBinCenter(singlebin_projY->GetMaximumBin());
        float width=1;
        float mintime = center-width;
        float maxtime = center+width;

        TF1 *f_gaus = new TF1("f_gaus","gaus",mintime,maxtime);
        //f_gaus->SetLineColor(kRed);
        //f_gaus->SetLineWidth(1);
        singlebin_projY->SetAxisRange(mintime,maxtime);
        double mean=0, meanerr=0;

        singlebin_projY->Fit(f_gaus,"RQ");
        mean = f_gaus->GetParameter(1);
        meanerr = f_gaus->GetParError(1);
        xmean->SetBinContent(i,mean);
        xmean->SetBinError(i,meanerr);
        if (makeplots || meanerr>0.5) {
            char plotname[255];
            sprintf(plotname,"plots/gausfitformean/fit_%s_%i.png",tag.Data(),i);
            canvas_fit->Print(plotname);
        }
        singlebin_projY->Reset();
    }
    return xmean;
}
