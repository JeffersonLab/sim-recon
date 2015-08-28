TFile * thisFile;

TH1I * Get1DHistogram(const char * plugin, const char * directoryName, const char * name){
    TH1I * histogram;
    TString fullName = TString(plugin) + "/" + TString(directoryName) + "/" + TString(name);
    thisFile->GetObject(fullName, histogram);
    if (histogram == 0){
        cout << "Unable to find histogram " << fullName.Data() << endl;
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

int GetCCDBIndexTAGM(unsigned int column, unsigned int row){
    int CCDBIndex = column + row;
    if (column > 9) CCDBIndex += 5;
    if (column > 27) CCDBIndex += 5;
    if (column > 81) CCDBIndex += 5;
    if (column > 99) CCDBIndex += 5;

    return CCDBIndex;
}

void ExtractTrackBasedTiming(int runNumber){

    TString fileName = Form ("Run%i/TrackBasedTiming.root", runNumber);
    TString prefix = Form ("Run%i/constants/TrackBasedTiming/",runNumber);
    TString inputPrefix = Form ("Run%i/constants/TDCADCTiming/",runNumber);

    thisFile = TFile::Open( fileName , "UPDATE");
    if (thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }

    //We need the existing constants, The best we can do here is just read them from the file.
    vector<double> sc_tdc_time_offsets;
    vector<double> sc_fadc_time_offsets;
    vector<double> tof_tdc_time_offsets;
    vector<double> tof_fadc_time_offsets;
    vector<double> tagm_tdc_time_offsets;
    vector<double> tagm_fadc_time_offsets;
    vector<double> tagh_tdc_time_offsets;
    vector<double> tagh_fadc_time_offsets;

    double sc_t_base_fadc;
    double sc_t_base_tdc;
    double tof_t_base_fadc;
    double tof_t_base_tdc;
    double bcal_t_base_fadc;
    double bcal_t_base_tdc;
    double tagm_t_base_fadc;
    double tagm_t_base_tdc;
    double tagh_t_base_fadc;
    double tagh_t_base_tdc;
    double fcal_t_base;
    double cdc_t_base;

    ifstream inFile;
    inFile.open(inputPrefix + "sc_tdc_timing_offsets.txt");
    string line;
    if (inFile.is_open()){
        while (getline (inFile, line)){
            sc_tdc_time_offsets.push_back(atof(line.data()));
        }
    }
    inFile.close();

    ifstream inFile;
    inFile.open(inputPrefix + "sc_adc_timing_offsets.txt");
    string line;
    if (inFile.is_open()){
        while (getline (inFile, line)){
            sc_fadc_time_offsets.push_back(atof(line.data()));
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tof_tdc_timing_offsets.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            tof_tdc_time_offsets.push_back(atof(line.data()));
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tof_adc_timing_offsets.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            tof_fadc_time_offsets.push_back(atof(line.data()));
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tagm_tdc_timing_offsets.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            double r, c, offset;
            while (iss>>r>>c>>offset){
                //if (row != 0) continue;
                tagm_tdc_time_offsets.push_back(offset);
            }
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tagm_adc_timing_offsets.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            double r, c, offset;
            while (iss>>r>>c>>offset){
                //if (row != 0) continue;
                tagm_fadc_time_offsets.push_back(offset);
            }
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tagh_tdc_timing_offsets.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            double counter, offset;
            while (iss>>counter>>offset){
                tagh_tdc_time_offsets.push_back(offset);
            }
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tagh_adc_timing_offsets.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            double counter, offset;
            while (iss>>counter>>offset){
                tagh_fadc_time_offsets.push_back(offset);
            }
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tof_base_time.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            iss>>tof_t_base_fadc>>tof_t_base_tdc;
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "sc_base_time.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            iss>>sc_t_base_fadc>>sc_t_base_tdc;
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "bcal_base_time.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            double adc_offset, tdc_offset;
            iss>>adc_offset>>tdc_offset; // TDC not used currently
            bcal_t_base_fadc = adc_offset;
            bcal_t_base_tdc = tdc_offset;
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "tagm_base_time.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            double adc_offset, tdc_offset;
            iss>>adc_offset>>tdc_offset; // TDC not used currently
            tagm_t_base_fadc = adc_offset;
            tagm_t_base_tdc = tdc_offset;
        }
    }

    inFile.close();
    inFile.open(inputPrefix + "tagh_base_time.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            double adc_offset, tdc_offset;
            iss>>adc_offset>>tdc_offset; // TDC not used currently
            tagh_t_base_fadc = adc_offset;
            tagh_t_base_tdc = tdc_offset;
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "fcal_base_time.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            iss>>fcal_t_base; 
        }
    }
    inFile.close();

    inFile.open(inputPrefix + "cdc_base_time.txt");
    if (inFile.is_open()){
        while (getline (inFile, line)){
            istringstream iss(line);
            iss>>cdc_t_base; 
        }
    }
    inFile.close();


    // Do our final step in the timing alignment with tracking

    //When the RF is present we can try to simply pick out the correct beam bucket for each of the runs
    //First just a simple check to see if we have the appropriate data
    bool useRF = false;
    double RF_Period = 4.0080161;
    TH1I *testHist = Get1DHistogram("HLDetectorTiming", "TAGH_TDC_RF_Compare","Counter ID 001");
    if (testHist != NULL){ // Not great since we rely on channel 1 working, but can be craftier later.
        useRF = true;
    }
    ofstream outFile;
    TH2I *thisHist; 
    thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - SC Target Time");
    if (useRF) thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGM - RFBunch Time");
    if (thisHist != NULL){
        //Statistics on these histograms are really quite low we will have to rebin and do some interpolation
        outFile.open(prefix + "tagm_tdc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        outFile.open(prefix + "tagm_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        TH1D * selectedTAGMOffset = new TH1D("selectedTAGMOffset", "Selected TAGM Offset; Column; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        TH1I * TAGMOffsetDistribution = new TH1I("TAGMOffsetDistribution", "TAGM Offset; TAGM Offset [ns]; Entries", 500, -250, 250);
        for (int i = 1 ; i <= nBinsX; i++){ 
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            //chose the correct number of bins based on the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 3; //ns (Full Width)
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
            //In the case there is RF, our job is to pick just the number of the correct beam bunch, so that's really all we need.
            if(useRF) {
                int beamBucket = int((maxMean / RF_Period) + 0.5); // +0.5 to handle rounding correctly
                selectedTAGMOffset->SetBinContent(i, beamBucket);
                TAGMOffsetDistribution->Fill(beamBucket);
            }
            else{
                selectedTAGMOffset->SetBinContent(i, maxMean);
                TAGMOffsetDistribution->Fill(maxMean);
            }
        }
        /*
        if (!useRF){
            //TFitResultPtr fr1 = selectedTAGMOffset->Fit("pol1", "SQ", "", 0.5, nBinsX + 0.5);
            TFitResultPtr fr1 = selectedTAGMOffset->Fit("pol1", "SQ", "", 5, 50);

            for (int i = 1 ; i <= nBinsX; i++){
                double x0 = fr1->Parameter(0);
                double x1 = fr1->Parameter(1);
                //double x2 = fr1->Parameter(2);
                //double fitResult = x0 + i*x1 + i*i*x2;
                double fitResult = x0 + i*x1;

                double outlierCut = 20;
                double valueToUse = selectedTAGMOffset->GetBinContent(i);
                if (fabs(selectedTAGMOffset->GetBinContent(i) - fitResult) > outlierCut && valueToUse != 0.0){
                    valueToUse = fitResult;
                }

                selectedTAGMOffset->SetBinContent(i, valueToUse);
                if (valueToUse != 0 ) TAGMOffsetDistribution->Fill(valueToUse);
            }
        }
*/
        double meanOffset = TAGMOffsetDistribution->GetMean();
        // This might be in units of beam bunches, so we need to convert
        if (useRF) meanOffset *= RF_Period;
        /*
           for (int i = 1 ; i <= nBinsX; i++){
           double valueToUse = selectedTAGMOffset->GetBinContent(i);
           if (useRF) valueToUse *= RF_Period;
           if (valueToUse == 0) valueToUse = meanOffset;
           outFile.open(prefix + "tagm_tdc_timing_offsets.txt", ios::out | ios::app);
           outFile << "0 " << i << " " << valueToUse + tagm_tdc_time_offsets[i-1] - meanOffset<< endl;
           if (i == 7 || i == 25 || i == 79 || i == 97){
           for(int j = 1; j <= 5; j++){
           outFile << j << " " << i << " " << valueToUse + tagm_tdc_time_offsets[i-1] - meanOffset<< endl;
           }
           }
           outFile.close();
        // Apply the same shift to the adc offsets
        outFile.open(prefix + "tagm_adc_timing_offsets.txt", ios::out | ios::app);
        outFile << "0 " << i << " " << valueToUse + tagm_fadc_time_offsets[i-1] - meanOffset<< endl;
        if (i == 7 || i == 25 || i == 79 || i == 97){
        for(int j = 1; j <= 5; j++){
        outFile << j << " " << i << " " << valueToUse + tagm_fadc_time_offsets[i-1] - meanOffset<< endl;
        }
        }
        outFile.close();
        }
        */

        outFile.open(prefix + "tagm_adc_timing_offsets.txt", ios::out);
        //for (int i = 1 ; i <= nBinsX; i++){
        // Loop over rows
        for (unsigned int column = 1; column <= 102; column++){
            int index = GetCCDBIndexTAGM(column, 0);
            double valueToUse = selectedTAGMOffset->GetBinContent(index);
            if (useRF) valueToUse *= RF_Period;
            if (valueToUse == 0) valueToUse = meanOffset;
            outFile << "0 " << column << " " << valueToUse + tagm_fadc_time_offsets[index-1] - meanOffset<< endl;
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    valueToUse = selectedTAGMOffset->GetBinContent(index);
                    if (useRF) valueToUse *= RF_Period;
                    if (valueToUse == 0) valueToUse = meanOffset;
                    outFile << row << " " << column << " " << valueToUse + tagm_fadc_time_offsets[index-1] - meanOffset<< endl;
                }
            }
        }
        outFile.close();

        outFile.open(prefix + "tagm_tdc_timing_offsets.txt", ios::out);
        //for (int i = 1 ; i <= nBinsX; i++){
        // Loop over rows
        for (unsigned int column = 1; column <= 102; column++){
            int index = GetCCDBIndexTAGM(column, 0);
            double valueToUse = selectedTAGMOffset->GetBinContent(index);
            if (useRF) valueToUse *= RF_Period;
            if (valueToUse == 0) valueToUse = meanOffset;
            outFile << "0 " << column << " " << valueToUse + tagm_tdc_time_offsets[index-1] - meanOffset << endl;
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    valueToUse = selectedTAGMOffset->GetBinContent(index);
                    if (useRF) valueToUse *= RF_Period;
                    if (valueToUse == 0) valueToUse = meanOffset;
                    outFile << row << " " << column << " " << valueToUse + tagm_tdc_time_offsets[index-1] - meanOffset << endl;
                }
            }
        }
        outFile.close();
        outFile.open(prefix + "tagm_base_time.txt", ios::out);
        outFile << tagm_t_base_fadc - meanOffset << " " << tagm_t_base_tdc - meanOffset << endl;
        outFile.close();

    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGH - SC Target Time");
    if (useRF) thisHist = Get2DHistogram("HLDetectorTiming", "TRACKING", "TAGH - RFBunch Time");
    if(thisHist != NULL){
        outFile.open(prefix + "tagh_tdc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file
        outFile.open(prefix + "tagh_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file

        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        TH1D * selectedTAGHOffset = new TH1D("selectedTAGHOffset", "Selected TAGH Offset; ID; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        TH1I * TAGHOffsetDistribution = new TH1I("TAGHOffsetDistribution", "TAGH Offset; TAGH Offset [ns]; Entries", 500, -250, 250);
        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            //chose the correct number of bins based on the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 2; //ns (Full Width)
            int binWindow = int(timeWindow / nsPerBin);

            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + binWindow) <= projY->GetNbinsX() ? (j + binWindow) : projY->GetNbinsX();
                double sum = 0; 
                double nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++){
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin){
                        if (nEntries > maxEntries){
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    }
                }
            }

            if(useRF) {
                int beamBucket = int((maxMean / RF_Period) + 0.5); // +0.5 to handle rounding correctly
                selectedTAGHOffset->SetBinContent(i, beamBucket);
                TAGHOffsetDistribution->Fill(beamBucket);
            }
            else{
                selectedTAGHOffset->SetBinContent(i, maxMean);
            }
            /*
               outFile.open("tagh_tdc_timing_offsets.txt", ios::out | ios::app);
               outFile << i << " " << maxMean + tagh_tdc_time_offsets[i] << endl;
               outFile.close();
               outFile.open("tagh_adc_timing_offsets.txt", ios::out | ios::app);
               outFile << i << " " << maxMean + tagh_fadc_time_offsets[i] << endl;
               outFile.close();
               */
        }

        // Fit 1D histogram. If value is far from the fit use the fitted value
        // Two behaviors above and below microscope
        // This isn't working well, so removing...
        /*
           TFitResultPtr fr1 = selectedTAGHOffset->Fit("pol2", "SQ", "", 0.5, 131.5);
           TFitResultPtr fr2 = selectedTAGHOffset->Fit("pol2", "SQ", "", 182.5, 274.5);        

           for (int i = 1 ; i <= nBinsX; i++){
           double fitResult = 0.0;
           if (i < 150){
           double x0 = fr1->Parameter(0);
           double x1 = fr1->Parameter(1);
           double x2 = fr1->Parameter(2);
           fitResult = x0 + i*x1 + i*i*x2;
           }
           else{
           double x0 = fr2->Parameter(0);
           double x1 = fr2->Parameter(1);
           double x2 = fr2->Parameter(2);
           fitResult = x0 + i*x1 + i*i*x2;
           }

           double outlierCut = 7;
           double valueToUse = selectedTAGHOffset->GetBinContent(i);
           if (fabs(selectedTAGHOffset->GetBinContent(i) - fitResult) > outlierCut && valueToUse != 0.0){
           valueToUse = fitResult;
           }

           selectedTAGHOffset->SetBinContent(i, valueToUse);
           if(valueToUse != 0) TAGHOffsetDistribution->Fill(valueToUse);
           }
           */
        double meanOffset = TAGHOffsetDistribution->GetMean();
        if (useRF) meanOffset *= RF_Period;
        for (int i = 1 ; i <= nBinsX; i++){
            valueToUse = selectedTAGHOffset->GetBinContent(i);
            if (useRF) valueToUse *= RF_Period;
            if (valueToUse == 0) valueToUse = meanOffset;
            outFile.open(prefix + "tagh_tdc_timing_offsets.txt", ios::out | ios::app);
            outFile << i << " " << valueToUse + tagh_tdc_time_offsets[i-1] - meanOffset << endl;
            outFile.close();
            outFile.open(prefix + "tagh_adc_timing_offsets.txt", ios::out | ios::app);
            outFile << i << " " << valueToUse + tagh_fadc_time_offsets[i-1] - meanOffset << endl;
            outFile.close();
        }

        outFile.open(prefix + "tagh_base_time.txt", ios::out);
        outFile << tagh_t_base_fadc - meanOffset << " " << tagh_t_base_tdc - meanOffset << endl;
        outFile.close();
    }

    // We can use the RF time to calibrate the SC time (Experimental for now)
    double meanSCOffset = 0.0; // In case we change the time of the SC, we need this in this scope
    if(useRF){
        TH1F * selectedSCSectorOffset = new TH1F("selectedSCSectorOffset", "Selected TDC-RF offset;Sector; Time", 30, 0.5, 30.5);
        TH1F * selectedSCSectorOffsetDistribution = new TH1F("selectedSCSectorOffsetDistribution", "Selected TDC-RF offset;Time;Entries", 100, -3.0, 3.0);
        TF1* f = new TF1("f","pol0(0)+gaus(1)", -3.0, 3.0);
        for (int sector = 1; sector <= 30; sector++){
            TH1I *scRFHist = Get1DHistogram("HLDetectorTiming", "SC_Target_RF_Compare", Form("Sector %.2i", sector));
            if (scRFHist == NULL) continue;
            //Do the fit
            TFitResultPtr fr = scRFHist->Fit("pol0", "SQ", "", -2, 2);
            double p0 = fr->Parameter(0);

            f->FixParameter(0,p0);
            f->SetParLimits(2, -2, 2);
            f->SetParLimits(3, 0, 2);
            f->SetParameter(1, 10);
            f->SetParameter(2, scRFHist->GetBinCenter(scRFHist->GetMaximumBin()));
            f->SetParameter(3, 0);

            fr = scRFHist->Fit(f, "SQ", "", -2, 2);
            double SCOffset = fr->Parameter(2);
            selectedSCSectorOffset->SetBinContent(sector, SCOffset);
            selectedSCSectorOffsetDistribution->Fill(SCOffset);
        }
        // Now write out the offsets
        meanSCOffset = selectedSCSectorOffsetDistribution->GetMean();
        outFile.open(prefix + "sc_tdc_timing_offsets.txt");
        for (int sector = 1; sector <= 30; sector++){
            outFile << sc_tdc_time_offsets[sector-1] + selectedSCSectorOffset->GetBinContent(sector) - meanSCOffset << endl;
        }
        outFile.close();
        outFile.open(prefix + "sc_adc_timing_offsets.txt");
        for (int sector = 1; sector <= 30; sector++){
            outFile << sc_fadc_time_offsets[sector-1] + selectedSCSectorOffset->GetBinContent(sector) - meanSCOffset << endl;
        }
        outFile.close();
        
        outFile.open(prefix + "sc_base_time.txt");
        outFile << sc_t_base_fadc - meanSCOffset << " " << sc_t_base_tdc - meanSCOffset << endl;
        outFile.close();
    }

    TH1I *this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "TOF - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 1.5, maximum + 1.5);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "tof_base_time.txt");
        outFile << tof_t_base_fadc - mean - meanSCOffset<< " " << tof_t_base_tdc - mean - meanSCOffset<< endl;
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "BCAL - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 5, maximum + 5);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "bcal_base_time.txt");
        outFile << bcal_t_base_fadc - mean - meanSCOffset << " " << bcal_t_base_tdc - mean - meanSCOffset << endl; // TDC info not used
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "FCAL - SC Target Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 5, maximum + 5);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "fcal_base_time.txt");
        outFile << fcal_t_base - mean - meanSCOffset<< endl; 
        outFile.close();
    }

    this1DHist = Get1DHistogram("HLDetectorTiming", "TRACKING", "Earliest CDC Time Minus Matched SC Time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 15, maximum + 10);
        float mean = fr->Parameter(1);
        outFile.open(prefix + "cdc_base_time.txt");
        outFile << cdc_t_base - mean - meanSCOffset << endl;
        outFile.close();
    }
    thisFile->Write();
    return;
    }
