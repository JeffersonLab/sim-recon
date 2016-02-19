using namespace RooFit;
TH1* GetHistogram(TFile *f,TString hname,TString htype,int counter) {
    stringstream ss; ss << counter;
    TH2I *h2 = (TH2I*)f->Get("PS_timing/"+hname);
    if (htype!="PSC"&&htype!="PS"&&htype!="TDCADC") cerr << "Unsupported histogram type: Choose 'PS', 'PSC', or 'TDCADC'" << endl;
    TH1 *h = h2->ProjectionY(TString(h2->GetName())+"_"+TString(ss.str()),counter,counter);
    TString title = (htype=="TDCADC") ? "PSC counter "+TString(ss.str()):htype+" counter "+TString(ss.str());
    h->SetTitle(title);
    h->GetXaxis()->SetTitle(h2->GetYaxis()->GetTitle());
    return h;
}
double GetMode(TH1 *h) {
    if (h->GetEntries() == 0.0) return 999.0;
    int max_bin = h->GetMaximumBin();
    double max = h->GetBinContent(max_bin);
    if (max < 7.0) return 999.0;
    return h->GetBinCenter(max_bin);
}
void WriteFitResults(ofstream &fout,TH1 *h,TString htype,int counter) {
    TString sep = ",";
    double max = h->GetBinContent(h->GetMaximumBin());
    if (h->GetEntries() < 100.0 || GetMode(h) == 999.0) {
        fout << counter << sep << max << sep << 0.0 << sep << 0.0 << sep << 0.0 << endl;
        return;
    }
    double w = (htype=="TDCADC") ? 1.0:0.75;
    double xlow = GetMode(h) - w;
    double xhigh = GetMode(h) + w;
    RooRealVar x("TDC time difference","TDC time difference [ns]",xlow,xhigh);
    RooDataHist data("data","data",RooArgList(x),h);
    // model: gaussian
    RooRealVar mean("mean","mean",GetMode(h),xlow,xhigh);
    RooRealVar sigma("sigma","sigma",0.2,0.01,0.5);
    RooGaussian gauss("gauss","gauss",x,mean,sigma);
    gauss.fitTo(data);
    // make plot
    TCanvas *canvas = new TCanvas("c","c",800,500);
    RooPlot *plot = x.frame();
    plot->SetTitle(h->GetTitle());
    plot->SetXTitle(h->GetXaxis()->GetTitle());
    TString ytitle = (htype=="PS") ? "PS hits":"PSC hits";
    plot->SetYTitle(ytitle);
    plot->SetTitleSize(0.0474,"XY");
    data.plotOn(plot);
    gauss.plotOn(plot,LineColor(kRed));
    plot->Draw();
    gauss.paramOn(plot,Layout(0.58,0.9,0.9),Format("NEU",AutoPrecision(1)));
    plot->Draw();
    stringstream ss; ss << counter;
    system(TString("mkdir -p fits_"+htype).Data());
    canvas->Print("fits_"+htype+"/counter_"+TString(ss.str())+".gif");
    delete canvas; delete plot;
    fout << counter << sep << max << sep << mean.getVal() << sep << mean.getError() << sep << sigma.getVal() << endl;
}
int fits(TString rootFile,bool doAllFits) {
    TFile *f = new TFile(rootFile,"read");
    TString hnames[] = {"PSC_tdcadcTimeDiffVsID","PSCRF_tdcTimeDiffVsID","PSRF_adcTimeDiffVsID"};
    TString htypes[] = {"TDCADC","PSC","PS"};
    int N_htypes = 2;
    if (doAllFits) N_htypes = 3;
    for (int htype=0;htype<N_htypes;htype++) {
        system("mkdir -p fits-csv");
        ofstream fout; fout.open("fits-csv/results_"+htypes[htype]+".txt");
        const int N = (htypes[htype]=="PS") ? 290:16;
        for (int i=1;i<=N;i++) { // counters
            TH1 *h = GetHistogram(f,hnames[htype],htypes[htype],i);
            WriteFitResults(fout,h,htypes[htype],i);
        }
        fout.close();
    }
    return 0;
}
