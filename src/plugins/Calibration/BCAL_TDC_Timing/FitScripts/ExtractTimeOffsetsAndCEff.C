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
void ExtractTimeOffsetsAndCEff(int run = 2931, TString filename = "hd_root.root"){

    // Open our input and output file
    thisFile = TFile::Open(filename);
    TFile *outputFile = TFile::Open("BCALTimeOffsetAndCEff_Results.root", "RECREATE");
    outputFile->mkdir("Fits");
    outputFile->mkdir("ResultOverview");

    // Check to make sure it is open
    if (thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }

    // We need to have the value of C_effective used for the position determination to get the time offsets
    double C_eff = 16.75; // cm/ns

    // Since we are updating existing constants we will need their current values. We can pipe them in from the CCDB...
    // We need a place to store them...
    double tdc_offsets[1152];
    double adc_offsets[1536];

    //Pipe the current constants into this macro
    //NOTE: This dumps the "LATEST" values. If you need something else, modify this script.
    char command[100];
    sprintf(command, "ccdb dump /BCAL/TDC_offsets:%i:default", run);
    FILE* locInputFile = gSystem->OpenPipe(command, "r");
    if(locInputFile == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), locInputFile) == NULL)
        return 0;
    //get the remaining lines
    double time;
    int counter = 0;
    while(fgets(buff, sizeof(buff), locInputFile) != NULL){
        istringstream locConstantsStream(buff);
        locConstantsStream >> time;
        cout << "TDC Time Offset = " << time << endl;
        tdc_offsets[counter] = time;
        counter++;
    }
    if (counter != 1152) cout << "Wrong number of TDC entries (" << counter << " != 1152)?" << endl;
    //Close the pipe
    gSystem->ClosePipe(locInputFile);

    sprintf(command, "ccdb dump /BCAL/ADC_timing_offsets:%i:default", run);
    locInputFile = gSystem->OpenPipe(command, "r");
    if(locInputFile == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), locInputFile) == NULL)
        return 0;
    //get the remaining lines
    counter = 0;
    while(fgets(buff, sizeof(buff), locInputFile) != NULL){
        istringstream locConstantsStream(buff);
        locConstantsStream >> time;
        cout << "ADC Time Offset = " << time << endl;
        adc_offsets[counter] = time;
        counter++;
    }
    if (counter != 1536) cout << "Wrong number of ADC entries (" << counter << " != 1536)?" << endl;
    //Close the pipe
    gSystem->ClosePipe(locInputFile);

    // This stream will be for outputting the results in a format suitable for the CCDB
    // Will wait to open until needed
    ofstream adcOffsetFile, tdcOffsetFile;
    adcOffsetFile.open("ADCOffsetsBCAL.txt");
    tdcOffsetFile.open("TDCOffsetsBCAL.txt");

    // Declaration of the fit funtion
    TF1 *f1 = new TF1("f1", "[0]+[1]*x", -200, 200);
    f1->SetParLimits(0, -20.0, 20.0);
    f1->SetParLimits(1, 0.85, 1.15);
    outputFile->cd("ResultOverview");
    // Make some histograms to get the distributions of the fit parameters
    TH1I *h1_c0 = new TH1I("h1_c0", "Distribution of parameter c_{0}", 100, -15, 15);
    TH1I *h1_c1 = new TH1I("h1_c1", "Distribution of parameter c_{1}", 100, 0.85, 1.15);
    TH1F *h1_c0_all = new TH1F ("h1_c0_all", "Value of c0; CCDB Index; c0 [cm]", 768, 0.5, 768.5);
    TH1F *h1_c1_all = new TH1F ("h1_c1_all", "Value of c1; CCDB Index; c1 [cm]", 768, 0.5, 768.5);
    TH2I *h2_c0_c1 = new TH2I("h2_c0_c1", "c_{1} Vs. c_{0}; c_{0}; c_{1}", 100, -15, 15, 100, 0.85, 1.15);
    outputFile->cd("Fits");
    // Now we want to loop through all available module/layer/sector and try to make a fit of each one
    for (unsigned int iModule = 1; iModule <=48; iModule++){
        for (unsigned int iLayer = 1; iLayer <= 4; iLayer++){ // Only 3 layers with TDCs
            for (unsigned int iSector = 1; iSector <= 4; iSector++){
                int the_cell = (iModule - 1) * 16 + (iLayer - 1) * 4 + iSector;
                int the_tdc_cell = (iModule - 1) * 12 + (iLayer - 1) * 4 + iSector; // One less layer of TDCs
                // Format the string to lookup the histogram by name
                char name[200];
                sprintf(name, "Module %.2i Layer %.2i Sector %.2i", iModule, iLayer, iSector);

                // These histograms are created on the fly in the plugin, so there is a chance that they do not exist, in which case the pointer will be NULL

                TH2I *h_offsets   = Get2DHistogram ("BCAL_TDC_Offsets", "Z Position", name);

                // Use FitSlicesY routine to extract the mean of each x bin
                TObjArray ySlices;

                if (h_offsets != NULL) {
                    h_offsets->RebinX(5);
                    TProfile *profile = h_offsets->ProfileX();
                    f1->SetParameters(0, 1); // Just out initial guess
                    TFitResultPtr fr = profile->Fit(f1, "SQR");
                    Int_t fitStatus = fr;
                    if (fitStatus == 0){
                        double c0 = fr->Parameter(0);
                        double c0_err = fr->ParError(0);
                        double c1 = fr->Parameter(1);
                        double c1_err = fr->ParError(1);
                        if (c0 == 10.0 || c0 == -10.0 || c1 == 0.9 || c1 == 1.1){
                            cout << "WARNING: Parameter hit limit " << name << endl;
                        }
                        h1_c0->Fill(c0); h1_c1->Fill(c1); h2_c0_c1->Fill(c0,c1); 
                        h1_c0_all->SetBinContent(the_cell, c0); h1_c0_all->SetBinError(the_cell, c0_err);
                        h1_c1_all->SetBinContent(the_cell, c1); h1_c1_all->SetBinError(the_cell, c1_err);
                        adcOffsetFile << adc_offsets[(the_cell - 1) * 2] + 0.5 * c0 / C_eff << endl;
                        adcOffsetFile << adc_offsets[ the_cell*2 - 1] - 0.5 * c0 / C_eff << endl;
                        if (iLayer != 4){
                            tdcOffsetFile << tdc_offsets[(the_tdc_cell - 1) * 2] + 0.5 * c0 / C_eff << endl;
                            tdcOffsetFile << tdc_offsets[the_tdc_cell*2 - 1] - 0.5 * c0 / C_eff << endl;
                        }
                    }
                    else {
                        cout << "WARNING: Fit Status "<< fitStatus << " for Upstream " << name << endl;
                        adcOffsetFile << adc_offsets[(the_cell - 1) * 2] << endl;
                        adcOffsetFile << adc_offsets[the_cell*2 - 1] << endl;
                        if (iLayer != 4){
                            tdcOffsetFile << tdc_offsets[(the_tdc_cell - 1) * 2] << endl;
                            tdcOffsetFile << tdc_offsets[the_tdc_cell*2 - 1] << endl;
                        }
                    }
                }
                else{
                    adcOffsetFile << adc_offsets[ (the_cell-1) * 2] << endl;
                    adcOffsetFile << adc_offsets[  the_cell*2  - 1] << endl;
                    if (iLayer != 4){
                        tdcOffsetFile << tdc_offsets[(the_tdc_cell -1) * 2] << endl;
                        tdcOffsetFile << tdc_offsets[ the_tdc_cell*2 - 1] << endl;
                    }
                }
            }
        }
    }
    adcOffsetFile.close();
    tdcOffsetFile.close();
    outputFile->Write();
    thisFile->Close();
    return;
}



