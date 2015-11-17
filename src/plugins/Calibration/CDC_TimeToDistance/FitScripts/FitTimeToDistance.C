double long_drift_func[3][3];
double short_drift_func[3][3];
double magnet_correction[2][2];

vector<double> cdc_drift_table;

unsigned int Locate(vector<double>&xx, double x){
    int n=xx.size();
    if (x==xx[0]) return 0;
    else if (x==xx[n-1]) return n-2;

    int jl=-1;
    int ju=n;
    int ascnd=(xx[n-1]>=xx[0]);
    while(ju-jl>1){
        int jm=(ju+jl)>>1;
        if ( (x>=xx[jm])==ascnd)
            jl=jm;
        else
            ju=jm;
    } 
    return jl;
}

Double_t TimeToDistance( Double_t *x, Double_t *par)
{
    double d=0.;
    double delta = x[1]; // yAxis
    double t = x[0]; // xAxis 
    if (t>0){ // For this fit it always will be, but just for completeness
        double f_0=0.;
        double f_delta=0.;

        if (delta > 0){
            double a1=par[0];
            double a2=par[1];
            double a3=par[2];
            double b1=par[3];
            double b2=par[4];
            double b3=par[5];
            double c1=par[6];
            double c2=par[7];
            double c3=par[8];

            // use "long side" functional form
            double my_t=0.001*t;
            double sqrt_t=sqrt(my_t);
            double t3=my_t*my_t*my_t;
            double delta_mag=fabs(delta);
            f_delta=(a1+a2*delta_mag)*sqrt_t+(b1+b2*delta_mag)*my_t
                +(c1+c2*delta_mag+c3*delta*delta)*t3;
            f_0=a1*sqrt_t+b1*my_t+c1*t3;
        }
        else{
            double my_t=0.001*t;
            double sqrt_t=sqrt(my_t);
            double delta_mag=fabs(delta);

            // use "short side" functional form
            double a1=par[9];
            double a2=par[10];
            double a3=par[11];
            double b1=par[12];
            double b2=par[13];
            double b3=par[14];
            double c1=par[15];
            double c2=par[16];
            double c3=par[17];

            double delta_sq=delta*delta;
            f_delta= (a1+a2*delta_mag+a3*delta_sq)*sqrt_t
                +(b1+b2*delta_mag+b3*delta_sq)*my_t;
            f_0=a1*sqrt_t+b1*my_t;
        }

        unsigned int max_index=cdc_drift_table.size()-1;
        if (t>cdc_drift_table[max_index]){
            d=f_delta;
            return d;
        }

        // Drift time is within range of table -- interpolate...
        unsigned int index=0;
        index=Locate(cdc_drift_table,t);
        double dt=cdc_drift_table[index+1]-cdc_drift_table[index];
        double frac=(t-cdc_drift_table[index])/dt;
        double d_0=0.01*(double(index)+frac); 

        double P=0.;
        double tcut=250.0; // ns
        if (t<tcut) {
            P=(tcut-t)/tcut;
        }
        d=f_delta*(d_0/f_0*P+1.-P);
    }
    return d;
}


void FitTimeToDistance(TString inputROOTFile = "hd_root.root", int run = 3650)
{
    // Script for fitting the time to distance relation from data
    // We need the values from the database to serve as our starting point

    //Pipe the current constants into this macro
    //NOTE: This dumps the "LATEST" values. If you need something else, modify this script.
    char command[100];
    sprintf(command, "ccdb dump /CDC/drift_parameters:%i:NoBField", run);
    FILE* locInputFile = gSystem->OpenPipe(command, "r");
    if(locInputFile == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), locInputFile) == NULL)
        return 0;
    //get the remaining lines
    double value;
    bool is_long = true;
    while(fgets(buff, sizeof(buff), locInputFile) != NULL){
        istringstream locConstantsStream(buff);
        if (is_long){
            locConstantsStream >> long_drift_func[0][0] >> long_drift_func[0][1] >> long_drift_func[0][2];
            locConstantsStream >> long_drift_func[1][0] >> long_drift_func[1][1] >> long_drift_func[1][2];
            locConstantsStream >> long_drift_func[2][0] >> long_drift_func[2][1] >> long_drift_func[2][2];
            locConstantsStream >> magnet_correction[0][0] >> magnet_correction[0][1]; // BField params
            cout << "a1_Long = " << long_drift_func[0][0] << endl;
            is_long = false;
        }
        else{
            locConstantsStream >> short_drift_func[0][0] >> short_drift_func[0][1] >> short_drift_func[0][2];
            locConstantsStream >> short_drift_func[1][0] >> short_drift_func[1][1] >> short_drift_func[1][2];
            locConstantsStream >> short_drift_func[2][0] >> short_drift_func[2][1] >> short_drift_func[2][2];
            locConstantsStream >> magnet_correction[1][0] >> magnet_correction[1][1];
            cout << "a1_Short = " << short_drift_func[0][0] << endl;
        }
    }
    //Close the pipe
    gSystem->ClosePipe(locInputFile);

    sprintf(command, "ccdb dump /CDC/cdc_drift_table:%i:NoBField", run);
    FILE* locInputFile = gSystem->OpenPipe(command, "r");
    if(locInputFile == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), locInputFile) == NULL)
        return 0;
    //get the remaining lines
    double value;
    while(fgets(buff, sizeof(buff), locInputFile) != NULL){
        istringstream locConstantsStream(buff);
        while (locConstantsStream >> value){
            cdc_drift_table.push_back(1000. * value);
            cout << "Drift Table value = " << value << endl;
        }
    }
    //Close the pipe
    gSystem->ClosePipe(locInputFile);
    // So now you have the input to your function

    const Int_t npar = 18;
    Double_t parameters[npar] =
    {long_drift_func[0][0], long_drift_func[0][1], long_drift_func[0][2],
        long_drift_func[1][0], long_drift_func[1][1], long_drift_func[1][2],
        long_drift_func[2][0], long_drift_func[2][1], long_drift_func[2][2],
        short_drift_func[0][0], short_drift_func[0][1], short_drift_func[0][2],
        short_drift_func[1][0], short_drift_func[1][1], short_drift_func[1][2],
        short_drift_func[2][0], short_drift_func[2][1], short_drift_func[2][2]};

    TF2 *f2 = new TF2("f2",TimeToDistance, 0, 1000, -0.3, 0.3, npar);
    f2->SetParameters(parameters);
    
    TFile *thisFile = TFile::Open(inputROOTFile);
    TProfile2D *profile = (TProfile2D *) thisFile->Get("/CDC_TimeToDistance/Predicted Drift Distance Vs Delta Vs t_drift");

    //profile->Fit("f2");
    gStyle->SetOptStat(0);

    TCanvas *c1 = new TCanvas ("c1", "c1", 800, 600);
    Double_t contours[21] = 
    { 0.00, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 
        0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0};

    profile->SetContour(21, contours);
    profile->Draw("colz");
    f2->SetContour(21, contours);
    profile->Draw("cont2 list same");
    f2->Draw("cont2 list same");
    c1->Update();
    c1->SaveAs("Before.png");

    TProfile2D *profileRebin = profile->Rebin2D(4,4, "Rebin");
    TFitResultPtr fr = profileRebin->Fit("f2", "S");
    TCanvas *c2 = new TCanvas ("c2", "c2", 800, 600);
    c2->cd();
    profile->Draw("colz");
    f2->SetContour(21, contours);
    profile->Draw("cont2 list same");
    f2->Draw("cont2 list same");
    c2->Update();
    c2->SaveAs("After.png");

    ofstream outputTextFile;
    outputTextFile.open("ccdb_Format.txt"); 
    outputTextFile << fr->Parameter(0) << " " << fr->Parameter(1) << " " << fr->Parameter(2) << " " ;
    outputTextFile << fr->Parameter(3) << " " << fr->Parameter(4) << " " << fr->Parameter(5) << " " ;
    outputTextFile << fr->Parameter(6) << " " << fr->Parameter(7) << " " << fr->Parameter(8) << " " ;
    outputTextFile << magnet_correction[0][0] << " " << magnet_correction[0][1] << endl; 
    outputTextFile << fr->Parameter(9) << " " << fr->Parameter(10) << " " << fr->Parameter(11) << " " ;
    outputTextFile << fr->Parameter(12) << " " << fr->Parameter(13) << " " << fr->Parameter(14) << " " ;
    outputTextFile << fr->Parameter(15) << " " << fr->Parameter(16) << " " << fr->Parameter(17) << " " ;
    outputTextFile << magnet_correction[1][0] << " " << magnet_correction[1][1] << endl;
    outputTextFile.close();

    return;
}
