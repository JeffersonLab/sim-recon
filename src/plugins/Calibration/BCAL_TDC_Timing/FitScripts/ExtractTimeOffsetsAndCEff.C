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
    double tdiff_u_d[768];
    double channel_global_offset[768];

    //Pipe the current constants into this macro
    //NOTE: This dumps the "LATEST" values. If you need something else, modify this script.
    char command[100];
    sprintf(command, "ccdb dump /BCAL/tdiff_u_d:%i:default", run);
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
        cout << "t_up - t_down = " << time << endl;
        tdiff_u_d[counter] = time;
        counter++;
    }
    if (counter != 768) cout << "Wrong number of tdiff_u_d entries (" << counter << " != 768)?" << endl;
    //Close the pipe
    gSystem->ClosePipe(locInputFile);

    sprintf(command, "ccdb dump /BCAL/channel_global_offset:%i:default", run);
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
        cout << "Channel Global Offset = " << time << endl;
        channel_global_offset[counter] = time;
        counter++;
    }
    if (counter != 768) cout << "Wrong number of channel global offsets (" << counter << " != 768)?" << endl;
    //Close the pipe
    gSystem->ClosePipe(locInputFile);

    // This stream will be for outputting the results in a format suitable for the CCDB
    // Will wait to open until needed
    ofstream channel_global_offset_File, tdiff_u_d_File;
    channel_global_offset_File.open("channel_global_offset_BCAL.txt");
    tdiff_u_d_File.open("tdiff_u_d_BCAL.txt");

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

    // Fit the global offset histogram to get the per channel global offset
    TH1D * selectedBCALOffset = new TH1D("selectedBCALOffset", "Selected Global BCAL Offset; CCDB Index; Offset [ns]", 768, 0.5, 768 + 0.5);
    TH1I * BCALOffsetDistribution = new TH1I("BCALOffsetDistribution", "Global BCAL Offset; Global Offset [ns]; Entries", 100, -10, 10);

    thisHist = Get2DHistogram("BCAL_Global_Offsets", "Target Time", "Target Time Minus RF Time Vs. Cell Number");
    if(thisHist != NULL){
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
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
            channel_global_offset_File <<  maxMean + channel_global_offset[i-1] << endl;;
            selectedBCALOffset->SetBinContent(i, maxMean);
            BCALOffsetDistribution->Fill(maxMean);
        }
    }

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
                        tdiff_u_d_File << tdiff_u_d[the_cell - 1] + c0 / C_eff << endl;
                    }
                    else {
                        cout << "WARNING: Fit Status "<< fitStatus << " for Upstream " << name << endl;
                        tdiff_u_d_File << tdiff_u_d[the_cell - 1] << endl;
                    }
                }
                else{
                    tdiff_u_d_File << tdiff_u_d[the_cell - 1] << endl;
                }
            }
        }
    }
    channel_global_offset_File.close();
    tdiff_u_d_File.close();
    outputFile->Write();
    thisFile->Close();
    return;
}



