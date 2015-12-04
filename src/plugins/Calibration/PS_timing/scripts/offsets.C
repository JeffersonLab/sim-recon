void GetData(TString fname,double *y,double &mean,const int N) {
    ifstream fin(fname);
    double prod = 0.0; double sum = 0.0;
    for (int i=0;i<N;i++) {
        char sep;int n; double w,t,dt,sigma;
        fin >> n >> sep >> w >> sep >> t >> sep >> dt >> sep >> sigma;
        y[i] = t;
        prod += w*t;
        sum += w;
    }
    if (sum>0.0) mean = prod/sum;
    for (int i=0;i<N;i++) {if (y[i]==0.0) y[i] = mean;}
    fin.close();
}
void WriteBaseOffsets(double base_tdc,double base_tdcadc,double base_ps) {
    ofstream fout; fout.open("offsets/base_time_offset.txt");
    TString sep = "        ";
    fout << -base_tdc + base_tdcadc << sep << -base_tdc << sep << -base_ps << endl;
    fout.close();
}
void WritePSCTDCOffsets(double *y,double mean,const int N) {
    ofstream fout; fout.open("offsets/tdc_timing_offsets_psc.txt");
    TString sep = "        ";
    for (int i=0;i<N;i++) fout << y[i]-mean << endl;
    fout.close();
}
void WritePSCADCOffsets(double *y_tdcadc,double mean_tdcadc,double *y,double mean,const int N) {
    ofstream fout; fout.open("offsets/adc_timing_offsets_psc.txt");
    TString sep = "        ";
    for (int i=0;i<N;i++) fout << y[i]-mean-(y_tdcadc[i]-mean_tdcadc) << endl;
    fout.close();
}
void WritePSADCOffsets(double *y,double mean,const int N) {
    ofstream fout; fout.open("offsets/adc_timing_offsets_ps.txt");
    TString sep = "        ";
    for (int i=0;i<int(N/2);i++) fout << y[i]-mean << sep << y[i+145]-mean << endl;
    fout.close();
}
int offsets(TString dir) {
    system("mkdir -p offsets");
    const int N = 16; // counters for PSC offsets
    double y[N]; double mean_tdc = 0.0;
    GetData(dir+"/results_PSC.txt",y,mean_tdc,N);
    WritePSCTDCOffsets(y,mean_tdc,N);
    double y_tdcadc[N]; double mean_tdcadc = 0.0;
    GetData(dir+"/results_TDCADC.txt",y_tdcadc,mean_tdcadc,N);
    WritePSCADCOffsets(y_tdcadc,mean_tdcadc,y,mean_tdc,N);
    const int N_ps = 290; // counters for PS offsets
    double y_ps[N_ps]; double mean_ps = 0.0;
    GetData(dir+"/results_PS.txt",y_ps,mean_ps,N_ps);
    WritePSADCOffsets(y_ps,mean_ps,N_ps);
    WriteBaseOffsets(mean_tdc,mean_tdcadc,mean_ps);
    return 0;
}
