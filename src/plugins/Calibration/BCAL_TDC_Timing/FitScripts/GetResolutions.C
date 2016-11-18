// Script to extract time-walk constants for the BCAL

#include <TFile.h>
#include <TDirectory.h>

namespace GetResolutionsNS {
    //Leave this global so the accesors don't need the pointer as an argument
    TFile *thisFile=NULL;

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

    TH2I * Get2DHistogram(const char * plugin, const char * directoryName, const char * name){
        TH2I * histogram;
        TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
        thisFile->GetObject(fullName, histogram);
        if (histogram == 0){
            cout << "Unable to find histogram " << fullName.Data() << endl;
            return NULL;
        }
        //cout << "Found histogram " << fullName.Data() << endl;
        return histogram;
    }

    TH2D * Get2DWeightedHistogram(const char * plugin, const char * directoryName, const char * name){
        TH2D * histogram;
        TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
        thisFile->GetObject(fullName, histogram);
        if (histogram == 0){
            cout << "Unable to find histogram " << fullName.Data() << endl;
            return NULL;
        }
        //cout << "Found histogram " << fullName.Data() << endl;
        return histogram;
    }
}

void GetXYfromSize(int size, int &x, int &y) {
    x=1;
    y=1;
    while (x*y<size) {
        if (x==y) x++;
        else 
            if (x>y) y++;
    }
}


void GetResolutionHistograms(TH2 *hist2d, TH1D *rmshist, TH1D *fithist, TH1D *meanhist, TH1D *fitmeanhist, int minhits, TString tag, bool print);
int GetColor(int number); 
void FormatHist(TH1* hist, int color, int markerstyle);

void makeplot(std::string plot_title, TString tag,  std::string plot_label,
              std::vector<std::string> histname_plot, 
              std::vector<std::string> histdir_plot, std::string plot_type,
              std::vector<std::string> legend_plot, bool logx=1) {
    char canvasname[255];
    sprintf(canvasname,"canvas_%s_%s",plot_label.c_str(),plot_type.c_str());
    TCanvas *canvas_showers_sigma = new TCanvas(canvasname,canvasname);
    gPad->SetLogx(logx);
    gPad->SetGridx();
    gPad->SetGridy();

    char fullhistname[255], plotfilename[255];
    TLegend *leg_showers_sigma = new TLegend(0.75,0.75,0.99,0.99);
    for (int i=0; i<histname_plot.size(); i++) {
        sprintf(fullhistname,"%s_%s_%s",histdir_plot[i].c_str(),histname_plot[i].c_str(),plot_type.c_str());
        TH1D* histpointer = (TH1D*)gDirectory->Get(fullhistname);
        if (histpointer==NULL) {
            printf("    no %s\n",fullhistname);
        } else {
            if (i==0) {
                histpointer->SetTitle(plot_title.c_str());
                histpointer->Draw();
            }
            else histpointer->Draw("same");
            if (logx) histpointer->SetAxisRange(0.02,5);
            leg_showers_sigma->AddEntry(histpointer,legend_plot[i].c_str(),"ep");
        }
    }
    leg_showers_sigma->Draw();
    sprintf(plotfilename,"plots/results/%s_%s_%s.png",plot_label.c_str(),plot_type.c_str(),tag.Data());
    canvas_showers_sigma->Print(plotfilename);
}


// Do the extraction of the actual constants
void GetResolutions(int run = 2931, TString filename = "hd_root.root", TString tag = "default", bool summary=1){
    printf("GetResolutions\n");
    // Open our input and output file
     printf("Opening for input \"%s\"\n",filename.Data());
    //TFile *fil = new TFile(filename.Data());
       GetResolutionsNS::thisFile = TFile::Open(filename.Data());
    // Check to make sure it is open
    if (    GetResolutionsNS::thisFile == NULL) {
        cout << "Unable to open file " << filename.Data() << "...Exiting" << endl;
        return;
    } else {
        cout << "Opened file " << filename.Data() << "...good" << endl;
    }
    //GetResolutionsNS::thisFile = fil;

    char outfilename[255];
    sprintf(outfilename,"Resolutions_Results_%s.root",tag.Data());
      printf("Opening for output %s\n",outfilename);
    TFile *outputFile = TFile::Open(outfilename, "RECREATE");
    //    outputFile->mkdir("Fits");
    outputFile->mkdir("ResultOverview");


    gStyle->SetOptStat(2220);
    gStyle->SetOptFit(1);

    int pad=0, color=0;

    // Showers

    std::vector<std::string> histname;
    std::vector<std::string> histdir;
    std::vector<int> minhits;

    // Make list of histograms to get resolutions from
    histdir.push_back("Showers");
    histname.push_back("PionShowers_q-");
    minhits.push_back(1000);
    histdir.push_back("Showers");
    histname.push_back("AllShowers_q0");
    minhits.push_back(1000);
    if (!summary) { // only do all the fits not in summary mode
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("AllPoints_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer1_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer2_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer3_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer4_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("AllPoints_q0");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer1_q0");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer2_q0");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer3_q0");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsEnergy");
        histname.push_back("Layer4_q0");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsShowerEnergy");
        histname.push_back("AllPoints_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsShowerEnergy");
        histname.push_back("Layer1_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsShowerEnergy");
        histname.push_back("Layer2_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsShowerEnergy");
        histname.push_back("Layer3_q-");
        minhits.push_back(1000);
        histdir.push_back("Points_deltaTVsShowerEnergy");
        histname.push_back("Layer4_q-");
        minhits.push_back(1000);
        histdir.push_back("Showers");
        histname.push_back("PionShowersVsP_q-");
        minhits.push_back(1000);
        histdir.push_back("Showers");
        histname.push_back("PionShowersVsZ_q-");
        minhits.push_back(1000);
        histdir.push_back("Showers");
        histname.push_back("AllShowersVsZ_q0");
        minhits.push_back(1000);


        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("AllPoints_TDC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("AllPoints_ADC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("Layer1_TDC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("Layer2_TDC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("Layer3_TDC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("Layer1_ADC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("Layer2_ADC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("Layer3_ADC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Points_deltaTVsEnergy");
        // histname.push_back("Layer4_ADC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Hits_deltaTVsE");
        // histname.push_back("AllHits_TDC_q-");
        // minhits.push_back(1000);
        // histdir.push_back("Hits_deltaTVsE");
        // histname.push_back("AllHits_ADC_q-");
        // minhits.push_back(1000);

        histdir.push_back("Target Time");
        histname.push_back("deltaTVsCell_q-");
        minhits.push_back(10);
        histdir.push_back("Target Time");
        histname.push_back("deltaTVsCell_q-_Eweight");
        minhits.push_back(10);
        histdir.push_back("Target Time");
        histname.push_back("deltaTVsCell_q0");
        minhits.push_back(10);
        histdir.push_back("Target Time");
        histname.push_back("deltaTVsCell_q0_Eweight");
        minhits.push_back(10);
    }

    std::vector<TH2*> TimeVsEnergy;
    std::vector<TH1D*> TimeSigma;
    std::vector<TH1D*> TimeRMS;
    std::vector<TH1D*> TimeMean;
    std::vector<TH1D*> TimeFitMean;

    // Get the resolutions for the histograms
    int xpads, ypads;
    GetXYfromSize(histname.size(),xpads,ypads);
    TCanvas *canvas_resolutions = new TCanvas("canvas_resolutions","canvas_resolutions",1200,800);
    canvas_resolutions->Divide(xpads,ypads,0.001,0.001);
    TCanvas *canvas_means = new TCanvas("canvas_means","canvas_means",1200,800);
    canvas_means->Divide(xpads,ypads,0.001,0.001);
    for (int i=0; i<histname.size(); i++) {
        auto locTimeVsEnergy = GetResolutionsNS::Get2DHistogram("BCAL_Global_Offsets", histdir[i].c_str(), histname[i].c_str());
        if (locTimeVsEnergy != NULL) {
            TimeVsEnergy.push_back(locTimeVsEnergy);
        } else {
            printf("Looking for a weighted histogram\n");
            auto locTimeVsEnergy = GetResolutionsNS::Get2DWeightedHistogram("BCAL_Global_Offsets", histdir[i].c_str(), histname[i].c_str());
            if (locTimeVsEnergy != NULL) {
                //printf("was able to recover\n");
                TimeVsEnergy.push_back(locTimeVsEnergy);
            }
        }
        //printf("size: %i\n",TimeVsEnergy.size());
        if (TimeVsEnergy[i] == NULL) {
            TimeSigma.push_back(NULL);
            TimeRMS.push_back(NULL);
            TimeMean.push_back(NULL);
            TimeFitMean.push_back(NULL);
        } else {
            //printf("in loop %s %s %s\n","BCAL_Global_Offsets", histdir[i].c_str(), histname[i].c_str());
            //gPad->SetLogx();
            color = GetColor(i+1);
            TAxis * xaxis = TimeVsEnergy[i]->GetXaxis();
            int nBinsX = TimeVsEnergy[i]->GetNbinsX();
            char fullhistname[255];
            sprintf(fullhistname,"%s_%s_sigma",histdir[i].c_str(),histname[i].c_str());
            //printf("%s\n",fullhistname);
            TH1D* locTimeSigma = TimeVsEnergy[i]->ProjectionX(fullhistname);
            TH1D* locTimeRMS = TimeVsEnergy[i]->ProjectionX(Form("%s_%s_RMS",histdir[i].c_str(),histname[i].c_str()));
            TH1D* locTimeMean = TimeVsEnergy[i]->ProjectionX(Form("%s_%s_mean",histdir[i].c_str(),histname[i].c_str()));
            TH1D* locTimeFitMean = TimeVsEnergy[i]->ProjectionX(Form("%s_%s_fitmean",histdir[i].c_str(),histname[i].c_str()));
            TimeSigma.push_back(locTimeSigma);
            TimeRMS.push_back(locTimeRMS);
            TimeMean.push_back(locTimeMean);
            TimeFitMean.push_back(locTimeFitMean);

            // TimeSigma[i] = new TH1D("TimeSigma","Charged shower time resolution;E  [GeV];T_{Target}-t_{RF}  [ns]",
            //                         nBinsX,xaxis->GetXmin(),xaxis->GetXmax());
            // TimeRMS[i] = new TH1D("TimeRMS","Charged shower time resolution;E  [GeV];T_{Target}-t_{RF}  [ns]",
            //                       nBinsX,xaxis->GetXmin(),xaxis->GetXmax());
            // TimeMean[i] = new TH1D("TimeMean","Charged shower time resolution;E  [GeV];T_{Target}-t_{RF}  [ns]",
            //                         nBinsX,xaxis->GetXmin(),xaxis->GetXmax());
            // TimeFitMean[i] = new TH1D("TimeFitMean","Charged shower time resolution;E  [GeV];T_{Target}-t_{RF}  [ns]",
            //                       nBinsX,xaxis->GetXmin(),xaxis->GetXmax());

            FormatHist(TimeSigma[i],color,20);
            TimeSigma[i]->GetYaxis()->SetTitle("#sigmat [ns]");
            FormatHist(TimeRMS[i],color,22);
            TimeRMS[i]->GetYaxis()->SetTitle("#sigmat [ns]");
            FormatHist(TimeMean[i],color,20);
            TimeMean[i]->GetYaxis()->SetTitle("#deltat [ns]");
            FormatHist(TimeFitMean[i],color,22);
            TimeFitMean[i]->GetYaxis()->SetTitle("#deltat [ns]");
            printf("%2i:  ",i);
            TString newtag = tag + "_" + histname[i].c_str();
            //printf("%s\n",newtag.Data());
            GetResolutionHistograms(TimeVsEnergy[i],TimeRMS[i],TimeSigma[i],TimeMean[i],TimeFitMean[i],minhits[i],newtag.Data(),0);
            canvas_resolutions->cd(i+1);
            TimeRMS[i]->Draw();
            TimeSigma[i]->Draw("same");
            canvas_means->cd(i+1);
            TimeFitMean[i]->Draw();
            //TimeMean[i]->Draw("same");
        } 
    }
    char plotfilename[255];
    if (!summary) {
        sprintf(plotfilename,"plots/results/resolution_pad_%s.png",tag.Data());
        canvas_resolutions->Print(plotfilename);
        sprintf(plotfilename,"plots/results/means_pad_%s.png",tag.Data());
        canvas_means->Print(plotfilename);
    }

    /*
      histdir_plot.push_back("Points_deltaTVsShowerEnergy");
        histname_plot.push_back("AllPoints_q-");
        legend_plot.push_back("");
        histdir_plot.push_back("Points_deltaTVsShowerEnergy");
        histname_plot.push_back("Layer1_q-");
        legend_plot.push_back("");
        histdir_plot.push_back("Points_deltaTVsShowerEnergy");
        histname_plot.push_back("Layer2_q-");
        legend_plot.push_back("");
        histdir_plot.push_back("Points_deltaTVsShowerEnergy");
        histname_plot.push_back("Layer3_q-");
        legend_plot.push_back("");
        histdir_plot.push_back("Points_deltaTVsShowerEnergy");
        histname_plot.push_back("Layer4_q-");
        legend_plot.push_back("");
*/

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);
    char title[255], fullhistname[255];
    std::string title_plot;

    std::vector<std::string> histname_plot;
    std::vector<std::string> histdir_plot;
    std::vector<std::string> legend_plot;

    histdir_plot.push_back("Showers");
    histname_plot.push_back("PionShowers_q-");
    legend_plot.push_back("Charged showers");
    histdir_plot.push_back("Showers");
    histname_plot.push_back("AllShowers_q0");
    legend_plot.push_back("Neutral showers");
    title_plot = "BCAL Shower time resolution";
    makeplot(title_plot,tag,"shower",histname_plot,histdir_plot,"sigma",legend_plot);
    title_plot = "BCAL Shower time offset";
    makeplot(title_plot,tag,"shower",histname_plot,histdir_plot,"fitmean",legend_plot);
    histdir_plot.clear();
    histname_plot.clear();
    legend_plot.clear();

    if (!summary) {
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("AllPoints_q-");
        legend_plot.push_back("All points");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer1_q-");
        legend_plot.push_back("layer 1");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer2_q-");
        legend_plot.push_back("layer 2");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer3_q-");
        legend_plot.push_back("layer 3");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer4_q-");
        legend_plot.push_back("layer 4");
        title_plot = "Charged points time resolution";
        makeplot(title_plot,tag,"charged_points",histname_plot,histdir_plot,"sigma",legend_plot);
        title_plot = "Charged points time offset";
        makeplot(title_plot,tag,"charged_points",histname_plot,histdir_plot,"fitmean",legend_plot);
        histdir_plot.clear();
        histname_plot.clear();
        legend_plot.clear();


        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("AllPoints_q0");
        legend_plot.push_back("All points");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer1_q0");
        legend_plot.push_back("Layer 1");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer2_q0");
        legend_plot.push_back("Layer 2");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer3_q0");
        legend_plot.push_back("Layer 3");
        histdir_plot.push_back("Points_deltaTVsEnergy");
        histname_plot.push_back("Layer4_q0");
        legend_plot.push_back("Layer 4");
        title_plot = "Neutral points time resolution";
        makeplot(title_plot,tag,"neutral_points",histname_plot,histdir_plot,"sigma",legend_plot);
        title_plot = "Neutral points time offset";
        makeplot(title_plot,tag,"neutral_points",histname_plot,histdir_plot,"fitmean",legend_plot);
        histdir_plot.clear();
        histname_plot.clear();
        legend_plot.clear();

        histdir_plot.push_back("Showers");
        histname_plot.push_back("AllShowersVsZ_q0");
        legend_plot.push_back("Neutral showers");
        histdir_plot.push_back("Showers");
        histname_plot.push_back("PionShowersVsZ_q-");
        legend_plot.push_back("Charged showers");
        title_plot = "Shower time resolution";
        makeplot(title_plot,tag,"shower_vsZ",histname_plot,histdir_plot,"sigma",legend_plot,0);
        title_plot = "Shower time offset";
        makeplot(title_plot,tag,"shower_vsZ",histname_plot,histdir_plot,"fitmean",legend_plot,0);
        histdir_plot.clear();
        histname_plot.clear();
        legend_plot.clear();

        histdir_plot.push_back("Showers");
        histname_plot.push_back("PionShowersVsP_q-");
        legend_plot.push_back("Charged showers");
        title_plot = "Shower time resolution";
        makeplot(title_plot,tag,"shower_vsP",histname_plot,histdir_plot,"sigma",legend_plot);
        title_plot = "Shower time offset";
        makeplot(title_plot,tag,"shower_vsP",histname_plot,histdir_plot,"fitmean",legend_plot);
        histdir_plot.clear();
        histname_plot.clear();
        legend_plot.clear();
    }

    outputFile->Write();
    GetResolutionsNS::thisFile->Close();
    return;
}






void GetResolutionHistograms(TH2 *hist2d, TH1D *rmshist, TH1D *fitsigmahist, TH1D *meanhist, TH1D *fitmeanhist, int minhits, 
                             TString tag, bool print) { 

    TAxis * xaxis = hist2d->GetXaxis();
    TAxis * yaxis = hist2d->GetYaxis();
    int nBinsX = hist2d->GetNbinsX();
    int nBinsY = hist2d->GetNbinsY();

    printf("bins (%i,%i) %s %s\n",nBinsX,nBinsY,hist2d->GetName(),hist2d->GetTitle());

    TH1D *mulitbin_projY = hist2d->ProjectionY("mulitbin_projY", 1, 1); // structure to accumulate data bin by bin until enough to fit
    mulitbin_projY->Reset();
    TH1D *xmean = hist2d->ProjectionX("xmean"); // structure to accumulate average X value of bin
    xmean->Reset();
    // reset the histograms that will be returned
    rmshist->Reset();
    fitsigmahist->Reset();
    meanhist->Reset();
    fitmeanhist->Reset();

    TCanvas *canvas_fit = new TCanvas("canvas_fit","canvas_fit");
    gPad->SetLogy();

    for (int i = 1 ; i <= nBinsX; i++){
        TH1D *singlebin_projY = hist2d->ProjectionY("temp", i, i);
        mulitbin_projY->Add(singlebin_projY); // accumulate bin
        float center = yaxis->GetBinCenter(mulitbin_projY->GetMaximumBin());
        // To avoid fitting the tails, only fit a region 1 ns around the peak
        float width=1; // in ns
        float mintime = center-width;
        float maxtime = center+width;

        char title[255];
        sprintf(title,"E=%.3f GeV (%s)",xaxis->GetBinCenter(i),tag.Data());
        mulitbin_projY->SetTitle(title);

        int entries_multi = mulitbin_projY->GetEntries();
        float int_multi = mulitbin_projY->Integral(yaxis->FindBin(mintime),yaxis->FindBin(maxtime));
        float int_single = singlebin_projY->Integral(yaxis->FindBin(mintime),yaxis->FindBin(maxtime));
        xmean->SetBinContent(i,int_single); // used to find weighted center of new composite bin
        float weightedmean = xmean->GetMean();
        //        int bin = xmean->GetXaxis()->FindBin(xmean);
        int bin = xmean->FindBin(weightedmean);
        //if (int_multi<minhits || (i==nBinsX && int_multi<(minhits/4.))) {
        if (entries_multi<minhits || (i==nBinsX && entries_multi<(minhits/4.))) {
            //printf("%4i %4i %8.1f %8.1f %6i  %6i %7.2f %2i\n",i,nBinsX,int_single,int_multi,entries_multi,minhits,weightedmean,bin);
            continue;
        }
        //if (int_multi != int_single && i!=bin) printf("%4i  %i   %i  %.2f  %i\n", i,int_single,int_multi,weightedmean,bin);

        mulitbin_projY->SetAxisRange(-5,5);
        double mean = mulitbin_projY->GetMean();
        double meanerr = mulitbin_projY->GetMeanError();
        double RMS = mulitbin_projY->GetRMS();
        double RMSerr = mulitbin_projY->GetRMSError();    
        rmshist->SetBinContent(bin,RMS);
        rmshist->SetBinError(bin,RMSerr);
        meanhist->SetBinContent(bin,mean);
        meanhist->SetBinError(bin,meanerr);

        mulitbin_projY->SetAxisRange(mintime,maxtime);        
        TF1 *f_gaus = new TF1("f_gaus","gaus",mintime,maxtime);
        f_gaus->SetLineColor(kRed);
        f_gaus->SetLineWidth(1);

        mulitbin_projY->Fit(f_gaus,"RQ");
        double fitmean = f_gaus->GetParameter(1);
        double fitmeanerr = f_gaus->GetParError(1);
        double sigma = f_gaus->GetParameter(2);
        double sigmaerr = f_gaus->GetParError(2);
        //printf("mean %6.3f %.3f  RMS %6.3f %.3f    fit mean %6.3f %.3f  sigma %6.3f %.3f\n",
        //       mean,meanerr,sigma,sigmaerr,fitmean,fitmeanerr,sigma,sigmaerr);
        fitsigmahist->SetBinContent(bin,sigma);
        fitsigmahist->SetBinError(bin,sigmaerr);
        fitmeanhist->SetBinContent(bin,fitmean);
        fitmeanhist->SetBinError(bin,fitmeanerr);
        if (print) {
            char plotname[255];
            gStyle->SetPadRightMargin(0.20);
            gStyle->SetStatW(0.20);
            mulitbin_projY->SetAxisRange(-10,10);
            sprintf(plotname,"plots/resolutions/fit_%s_%i.png",tag.Data(),i);
            canvas_fit->Print(plotname);
        }
        mulitbin_projY->Reset();
        xmean->Reset();
    }
    delete canvas_fit;
}



int GetColor(int number)
{
    int color = number;
    if (color > 4) color++;
    if (color > 9) color = 791-9+number;
    //printf("%i  %i\n",number,color);
    return color;
}

void FormatHist(TH1* hist, int color, int markerstyle)
{
    hist->SetMarkerStyle(markerstyle);
    hist->SetMarkerSize(0.8);
    hist->SetMarkerColor(color);
    hist->SetLineColor(color);
    hist->SetLineWidth(2);
}
