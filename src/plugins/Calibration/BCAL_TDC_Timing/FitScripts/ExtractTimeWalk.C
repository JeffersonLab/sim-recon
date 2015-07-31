// Script to extract time-walk constants for the BCAL

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

TH2I * Get2DHistogram(const char * plugin, const char * directoryName, const char * name){
    TH2I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    histogram = thisFile->GetObject(fullName, histogram);
    if (histogram == 0){
        cout << "Unable to find histogram " << fullName.Data() << endl;
        return NULL;
    }
    return histogram;
}

// Do the extraction of the actual constants
void ExtractTimeWalk(TString filename = "hd_root.root"){

    // Open our input and output file
    thisFile = TFile::Open(filename);
    TFile *outputFile = TFile::Open("BCALTimewalk_Results.root", "RECREATE");
    outputFile->mkdir("Upstream");
    outputFile->mkdir("Downstream");

    // Check to make sure it is open
    if (thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }

    // This stream will be for outputting the results in a format suitable for the CCDB
    // Will wait to open until needed
    ofstream textFile;
    textFile.open("TimewalkBCAL.txt");

    // Declaration of the fit funtion
    // We will use the one suggested in the "Low-level calibration constants for BCAL" GlueX-doc-2618
    // Using the pulse peak, the functional form is
    // t_TDC - tADC = c0 + c1 / ( pulse height / threshold ) ^ c2
    // The threshold in ADC counts is ~14 so this is in the denominator
    TF1 *f1 = new TF1("f1", "[0]+[1]/TMath::Power(x/14,[2])", 15, 400);
    f1->SetParLimits(2, 0.25, 1.5);

    // Make some histograms to get the distributions of the fit parameters
    TH1I *h1_c0 = new TH1I("h1_c0", "Distribution of parameter c_{0}", 100, -10, 10);
    TH1I *h1_c1 = new TH1I("h1_c1", "Distribution of parameter c_{1}", 100, 0, 20);
    TH1I *h1_c2 = new TH1I("h1_c2", "Distribution of parameter c_{2}", 100, 0.0, 1.5);
    TH2I *h2_c0_c1 = new TH2I("h2_c0_c1", "c_{1} Vs. c_{0}; c_{0}; c_{1}", 100, -10, 10, 100, 0, 20);
    TH2I *h2_c0_c2 = new TH2I("h2_c0_c2", "c_{2} Vs. c_{0}; c_{0}; c_{2}", 100, -10, 10, 100, 0.0, 1.5);
    TH2I *h2_c1_c2 = new TH2I("h2_c1_c2", "c_{2} Vs. c_{1}; c_{1}; c_{2}", 100, 0.0, 20, 100, 0.0, 1.5);

    // Now we want to loop through all available module/layer/sector and try to make a fit of each one
    for (unsigned int iModule = 1; iModule <=48; iModule++){
        for (unsigned int iLayer = 1; iLayer <= 3; iLayer++){ // Only 3 layers with TDCs
            for (unsigned int iSector = 1; iSector <= 4; iSector++){

                // Format the string to lookup the histogram by name
                char name[200];
                sprintf(name, "Module %.2i Layer %.2i Sector %.2i", iModule, iLayer, iSector);

                // There are four histograms we are possibly interested in. 
                // For each M/L/S we have Upstream and downstream, lets grab them
                // These histograms are created on the fly in the plugin, so there is a chance that they do not exist, in which case the pointer will be NULL

                TH2I *h_UpstreamTW_E   = Get2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_NoCorrection_E", name);
                TH2I *h_UpstreamTW_PP  = Get2DHistogram ("BCAL_TDC_Timing", "BCAL_Upstream_Timewalk_NoCorrection_PP", name);
                TH2I *h_DownstreamTW_E = Get2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_NoCorrection_E", name);
                TH2I *h_DownstreamTW_PP = Get2DHistogram ("BCAL_TDC_Timing", "BCAL_Downstream_Timewalk_NoCorrection_PP", name); 

                // Use FitSlicesY routine to extract the mean of each x bin
                TObjArray ySlicesUpstream;
                TObjArray ySlicesDownstream;

                outputFile->cd("Upstream");
                if (h_UpstreamTW_PP != NULL) {
                    h_UpstreamTW_PP->FitSlicesY(0, 0, -1, 0, "QNR", &ySlicesUpstream);
                    TH1D *meanHist = (TH1D *) ySlicesUpstream.At(1);
                    TH1D *meanHistClone = (TH1D *) meanHist->Clone(); // The other one will delete when the array goes out of scope
                    f1->SetParameters(0, 5, 0.7); // Just out initial guess
                    TFitResultPtr fr = meanHistClone->Fit(f1, "SRQ");
                    Int_t fitStatus = fr;
                    if (fitStatus == 0){
                        double c0 = fr->Parameter(0);
                        double c1 = fr->Parameter(1);
                        double c2 = fr->Parameter(2);
                        h1_c0->Fill(c0); h1_c1->Fill(c1); h1_c2->Fill(c2);
                        h2_c0_c1->Fill(c0,c1); h2_c0_c2->Fill(c0,c2); h2_c1_c2->Fill(c1,c2);
                        textFile << iModule << " " << iLayer << " " << iSector << " 0 " << c0 << " " << c1 << " " << c2 << " 0.0" << endl;
                    }
                    else {
                        cout << "WARNING: Fit Status "<< fitStatus << " for Upstream " << name << endl;
                        textFile << iModule << " " << iLayer << " " << iSector << " 0  0.0  0.0  0.0  0.0" << endl;
                    }
                }
                else{
                    textFile << iModule << " " << iLayer << " " << iSector << " 0  0.0  0.0  0.0  0.0" << endl;
                }

                outputFile->cd("Downstream");
                if (h_DownstreamTW_PP != NULL) {
                    h_DownstreamTW_PP->FitSlicesY(0, 0, -1, 0, "QNR", &ySlicesDownstream);
                    TH1D *meanHist = (TH1D *) ySlicesDownstream.At(1);
                    TH1D *meanHistClone = (TH1D *) meanHist->Clone(); // The other one will delete when the array goes out of scope
                    f1->SetParameters(0, 5, 0.7); // Just out initial guess
                    TFitResultPtr fr = meanHistClone->Fit(f1, "SRQ");
                    Int_t fitStatus = fr;
                    if (fitStatus == 0){
                        double c0 = fr->Parameter(0);
                        double c1 = fr->Parameter(1);
                        double c2 = fr->Parameter(2);
                        h1_c0->Fill(c0); h1_c1->Fill(c1); h1_c2->Fill(c2);
                        h2_c0_c1->Fill(c0,c1); h2_c0_c2->Fill(c0,c2); h2_c1_c2->Fill(c1,c2);
                        textFile << iModule << " " << iLayer << " " << iSector << " 1 " << c0 << " " << c1 << " " << c2 << " 0.0" << endl;
                    }
                    else {
                        cout << "WARNING: Fit Status "<< fitStatus << " for Downstream " << name << endl;
                        textFile << iModule << " " << iLayer << " " << iSector << " 1  0.0  0.0  0.0  0.0" << endl;
                    }
                }
                else{
                    textFile << iModule << " " << iLayer << " " << iSector << " 1  0.0  0.0  0.0  0.0" << endl;
                }
            }
        }
    }
    textFile.close();
    outputFile->Write();
    thisFile->Close();
    return;
}



