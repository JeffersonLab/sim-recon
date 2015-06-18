TFile * thisFile;

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

int GetCCDBIndexTAGM(int column, int row){
    int CCDBIndex = column + row;
    if (column > 9) CCDBIndex += 5;
    if (column > 27) CCDBIndex += 5;
    if (column > 81) CCDBIndex += 5;
    if (column > 99) CCDBIndex += 5;

    return CCDBIndex;
}

Double_t FitFunctionLeft(Double_t *x, Double_t *par)
{
    Float_t xx =x[0];
    Double_t f = par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1]) / par[2] , 2) ) )+ par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1] - 4.008) / par[2] , 2) ) );
    return f;
}

Double_t FitFunctionRight(Double_t *x, Double_t *par)
{
    Float_t xx =x[0];
    Double_t f = par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1]) / par[2] , 2) ) )+ par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1] + 4.008) / par[2] , 2) ) );
    return f;
}

void ExtractTDCADCTiming(int runNumber){

    TString fileName = Form ("Run%i/TDCADCTiming.root", runNumber);
    TString prefix = Form ("Run%i/constants/TDCADCTiming/",runNumber);

    thisFile = TFile::Open( fileName , "UPDATE");
    if (thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }
    ofstream outFile;

    //All the existing constants have been zero'd out already, But we should make copies in case one of the methods below fails

    outFile.open(prefix + "bcal_base_time.txt");
    outFile << "0.0 0.0" << endl;
    outFile.close();

    outFile.open(prefix + "cdc_base_time.txt");
    outFile << "0.0" << endl;
    outFile.close();

    outFile.open(prefix + "fdc_base_time.txt");
    outFile << "0.0 0.0" << endl;
    outFile.close();

    outFile.open(prefix + "fcal_base_time.txt");
    outFile << "0.0" << endl;
    outFile.close();

    outFile.open(prefix + "sc_base_time.txt");
    outFile << "0.0 0.0" << endl;
    outFile.close();

    outFile.open(prefix + "tagh_base_time.txt");
    outFile << "0.0 0.0" << endl;
    outFile.close();

    outFile.open(prefix + "tagm_base_time.txt");
    outFile << "0.0 0.0" << endl;
    outFile.close();

    outFile.open(prefix + "tof_base_time.txt");
    outFile << "0.0 0.0" << endl;
    outFile.close();

    outFile.open(prefix + "bcal_adc_timing_offsets.txt");
    for (int i = 1; i <= 1536; i++){
        outFile << "0.0" << endl;
    }
    outFile.close();

    outFile.open(prefix + "bcal_tdc_timing_offsets.txt");
    for (int i = 1; i <= 1152; i++){
        outFile << "0.0" << endl;
    }
    outFile.close();

    outFile.open(prefix + "fcal_adc_timing_offsets.txt");
    for (int i = 1; i <= 2800; i++){
        outFile << "0.0" << endl;
    }
    outFile.close();

    outFile.open(prefix + "sc_adc_timing_offsets.txt");
    for (int i = 1; i <= 30; i++){
        outFile << "0.0" << endl;
    }
    outFile.close();

    outFile.open(prefix + "sc_tdc_timing_offsets.txt");
    for (int i = 1; i <= 30; i++){
        outFile << "0.0" << endl;
    }
    outFile.close();

    outFile.open(prefix + "tagm_adc_timing_offsets.txt");
    for (int i = 1; i <= 102; i++){
        outFile << "0 " << i << " 0.0" << endl;
        if ( i == 7 || i == 25 || i == 79 || i == 97){
            for(int j = 1; j <= 5; j++){
                outFile << j << " " << i << " 0.0" << endl;
            }
        }
    }
    outFile.close();

    outFile.open(prefix + "tagm_tdc_timing_offsets.txt");
    for (int i = 1; i <= 102; i++){
        outFile << "0 " << i << " 0.0" << endl;
        if ( i == 7 || i == 25 || i == 79 || i == 97){
            for(int j = 1; j <= 5; j++){
                outFile << j << " " << i << " 0.0" << endl;
            }
        }
    }
    outFile.close();

    outFile.open(prefix + "tagh_adc_timing_offsets.txt");
    for (int i = 1; i <= 274; i++){
        outFile << i << " 0.0" << endl;
    }
    outFile.close();

    outFile.open(prefix + "tagh_tdc_timing_offsets.txt");
    for (int i = 1; i <= 274; i++){
        outFile << i << " 0.0" << endl;
    }
    outFile.close();

    outFile.open(prefix + "tof_adc_timing_offsets.txt");
    for (int i = 1; i <= 176; i++){
        outFile << "0.0" << endl;
    }
    outFile.close();

    //Move the base times just for the tracking

    float CDC_ADC_Offset = 0.0;
    TH1I * this1DHist = Get1DHistogram("HLDetectorTiming", "CDC", "CDCHit time");
    if(this1DHist != NULL){
        Int_t firstBin = this1DHist->FindLastBinAbove( 1 , 1); // Find first bin with content above 1 in the histogram
        for (int i = 0; i <= 25; i++){
            if ((firstBin + i) > 0) this1DHist->SetBinContent((firstBin + i), 0);
        }
        //Fit a gaussian to the left of the main peak
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TF1 *f = new TF1("f", "gaus");
        f->SetParameters(100, maximum, 20);
        f->FixParameter(1 , maximum);
        TFitResultPtr fr = this1DHist->Fit(f, "S", "", maximum - 40, maximum + 10); // Cant fix value at end of range
        double mean = fr->Parameter(1);
        float sigma = fr->Parameter(2);
        CDC_ADC_Offset = mean - sigma;
        delete f;
    }

    CDC_ADC_Offset *= -1;
    outFile.open(prefix + "cdc_base_time.txt");
    outFile << CDC_ADC_Offset << endl;
    outFile.close();

    float FDC_ADC_Offset = 0.0, FDC_TDC_Offset = 0.0;
    this1DHist = Get1DHistogram("HLDetectorTiming", "FDC", "FDCHit time");
    if(this1DHist != NULL){
        Int_t firstBin = this1DHist->FindLastBinAbove( 1 , 1); // Find first bin with content above 1 in the histogram
        for (int i = 0; i <= 25; i++){
            if ((firstBin + i) > 0) this1DHist->SetBinContent((firstBin + i), 0);
        }
        //Fit a gaussian to the left of the main peak
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TF1 *f = new TF1("f", "gaus");
        f->SetParameters(100, maximum, 20);
        f->FixParameter(1 , maximum);
        TFitResultPtr fr = this1DHist->Fit(f, "S", "", maximum - 40, maximum + 10); // Cant fix value at end of range
        double mean = fr->Parameter(1);
        float sigma = fr->Parameter(2);
        FDC_ADC_Offset = mean - sigma;
        delete f;
    }

    FDC_ADC_Offset *= -1;
    outFile.open(prefix + "fdc_base_time.txt");
    outFile << FDC_ADC_Offset << " " << FDC_TDC_Offset << endl;
    outFile.close();

    // Now that we have the file open, do all of the fits and write the output
    // Fit all plots with expected funtional form, output files for CCDB input
    // Do a finer alignement of the TDC and ADC's in the detectors that have both
    int minHits = 7;

    // In order to calibrate the SC in one step, we need to work with the base times and with the TDC/ADC offsets at the same time
    // Sort of complicates things but saves considerable time.

    float sc_tdc_base_time = 0.0;

    this1DHist = Get1DHistogram("HLDetectorTiming", "SC", "SCHit TDC time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 50, maximum + 50);
        double mean = fr->Parameter(1);
        sc_tdc_base_time = mean;
    }

    TH2I *thisHist = Get2DHistogram("HLDetectorTiming", "SC", "SCHit TDC_ADC Difference");
    TH1D * selected_SC_TDCADCOffset = NULL;
    TH1I * SCOffsetDistribution = NULL;
    if(thisHist != NULL){

        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();

        selected_SC_TDCADCOffset = new TH1D("selected_SC_TDCADCOffset", "Selected SC TDC/ADC Offset; CCDB Index; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        SCOffsetDistribution = new TH1I("SCOffsetDistribution", "SC ADC Offset; ADC Offset [ns]; Entries", 1000, -500, 500);

        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 2; //ns (Full Width)
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
            //outFile.open(prefix + "sc_adc_timing_offsets.txt", ios::out | ios::app);
            //outFile << -1 * maxMean << endl;
            //outFile.close();
            // Some histograms for tracking the magnitude of the shifts
            selected_SC_TDCADCOffset->SetBinContent(i, -1*maxMean);
            SCOffsetDistribution->Fill(-1*maxMean);
        }
    }

    // Two Output files to write for the SC
    float meanDiff = 0.0;
    if (selected_SC_TDCADCOffset != NULL){
        meanDiff = SCOffsetDistribution->GetMean();
        outFile.open(prefix + "sc_adc_timing_offsets.txt", ios::out);
        for( int iSC = 1; iSC <= 30; iSC++){ 
            float SC_ADC_Offset = selected_SC_TDCADCOffset->GetBinContent(iSC);
            outFile << SC_ADC_Offset - meanDiff << endl;
        }
        outFile.close();
    }

    outFile.open(prefix + "sc_base_time.txt", ios::out);
    outFile << -1*(sc_tdc_base_time + meanDiff) << " " << -1*sc_tdc_base_time << endl;
    outFile.close();

    thisHist = Get2DHistogram("HLDetectorTiming", "TOF", "TOFHit TDC_ADC Difference");
    if(thisHist != NULL){
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();

        TH1D * selected_TOF_TDCADCOffset = new TH1D("selected_TOF_TDCADCOffset", "Selected TOF TDC/ADC Offset; CCDB Index; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        TH1I * TOFOffsetDistribution = new TH1I("TOFOffsetDistribution", "TOF ADC Offset; ADC Offset [ns]; Entries", 1000, -500, 500);

        outFile.open(prefix + "tof_adc_timing_offsets.txt");
        outFile.close(); // Clear Output File

        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            //chose the correct number of bins based on the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 0.5; //ns (Full Width) -- Pretty kiler resolution
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
            // Some histograms for tracking the magnitude of the shifts
            selected_TOF_TDCADCOffset->SetBinContent(i, -1*maxMean);
            if (maxMean != 0) TOFOffsetDistribution->Fill(-1*maxMean);
        }

        double meanOffset = TOFOffsetDistribution->GetMean();
        outFile.open(prefix + "tof_adc_timing_offsets.txt", ios::out);
        for (int i = 1 ; i <= nBinsX; i++){
            double valueToUse = selected_TOF_TDCADCOffset->GetBinContent(i);
            if (valueToUse == 0) valueToUse = meanOffset;
            outFile << valueToUse - meanOffset << endl;
        }
        outFile.close();
        outFile.open(prefix + "tof_base_time.txt", ios::out);
        outFile << -1 * meanOffset << " 0.0" << endl;
        outFile.close();

    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TAGM", "TAGMHit TDC_ADC Difference");
    if(thisHist != NULL){
        float tdcRFOffset[122] = {}; // 122 possible offsets
        char name[100];
        //First loop through the TDC RF times to see which ones we can correct
        float period = 4.008161;
        TF1 * fitLeft = new TF1("fa",FitFunctionLeft, -1*period / 2, period / 2, 3);
        TF1 * fitRight = new TF1("fa",FitFunctionRight, -1*period / 2, period / 2, 3);
        for (unsigned int column = 1; column <= 102; column++){
            for (unsigned int row = 0; row <= 5; row++){
                sprintf(name, "Column %.3i Row %.1i", column, row);
                TH1I *rf_hist = Get1DHistogram("HLDetectorTiming", "TAGM_TDC_RF_Compare", name, false);
                if (rf_hist != NULL){
                    // do the fit and store the result
                    // Some fancy footwork to fit this periodic function
                    cout << "Fitting TAGM " << name << endl;
                    TF1 *fa;
                    if (rf_hist->GetBinCenter(rf_hist->GetMaximumBin()) < 0.0) fa = fitLeft;
                    else fa = fitRight;
                    fa->SetParLimits(1, -1*period / 2, period / 2);
                    fa->SetParameter(1,rf_hist->GetBinCenter(rf_hist->GetMaximumBin()));
                    fa->SetParLimits(2, 0.175, 1.2);
                    fa->SetParameter(2, 0.4);
                    TFitResultPtr r = rf_hist->Fit(fa, "S", "", -1*period / 2, period / 2);
                    Int_t status = r;
                    if ( status == 0){ // Fit succeeded
                        Double_t offset  = r->Parameter(1);
                        tdcRFOffset[ GetCCDBIndexTAGM(column, row) - 1 ] = offset;
                    }
                }
            }
        }

        //int nbins = meanHist->GetNbinsX();
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();

        TH1D * selectedTAGMOffset = new TH1D("selectedTAGMOffset", "Selected TAGM Offset; Column; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        TH1I * TAGMOffsetDistribution = new TH1I("TAGMOffsetDistribution", "TAGM ADC Offset; ADC Offset [ns]; Entries", 1000, -500, 500);

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
            selectedTAGMOffset->SetBinContent(i, -1 * maxMean);
            if (maxMean != 0) TAGMOffsetDistribution->Fill(-1 * maxMean);
        }
        //Put the average offset in the base times
        double meanOffset = TAGMOffsetDistribution->GetMean();
        outFile.open(prefix + "tagm_adc_timing_offsets.txt", ios::out);
        //for (int i = 1 ; i <= nBinsX; i++){
        // Loop over rows
        for (unsigned int column = 1; column <= 102; column++){
            int index = GetCCDBIndexTAGM(column, 0);
            double valueToUse = selectedTAGMOffset->GetBinContent(index);
            if (valueToUse == 0) valueToUse = meanOffset;
            outFile << "0 " << column << " " << valueToUse - meanOffset + tdcRFOffset[index-1] << endl;
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    valueToUse = selectedTAGMOffset->GetBinContent(index);
                    double valueToUse = selectedTAGMOffset->GetBinContent(index);
                    if (valueToUse == 0) valueToUse = meanOffset; 
                    outFile << row << " " << column << " " << valueToUse - meanOffset + tdcRFOffset[index-1] << endl;
                }
            }
        }
        outFile.close();

        outFile.open(prefix + "tagm_tdc_timing_offsets.txt", ios::out);
        //for (int i = 1 ; i <= nBinsX; i++){
        // Loop over rows
        for (unsigned int column = 1; column <= 102; column++){
            int index = GetCCDBIndexTAGM(column, 0);
            outFile << "0 " << column << " " << tdcRFOffset[index-1] << endl;
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    outFile << row << " " << column << " " << tdcRFOffset[index-1] << endl;
                }
            }
        }
        outFile.close();

        outFile.open(prefix + "tagm_base_time.txt", ios::out);
        outFile << -1 * meanOffset << " 0.0" << endl;
        outFile.close();
    }

    thisHist = Get2DHistogram("HLDetectorTiming", "TAGH", "TAGHHit TDC_ADC Difference");
    if(thisHist != NULL){
        // First take a look at the RF offsets
        float tdcRFOffsetTAGH[274] = {}; // 274 possible offsets
        char name[100];
        //First loop through the TDC RF times to see which ones we can correct
        float period = 4.008161;
        TF1 * fitLeft = new TF1("fa",FitFunctionLeft, -1*period / 2, period / 2, 3);
        TF1 * fitRight = new TF1("fa",FitFunctionRight, -1*period / 2, period / 2, 3);
        for (unsigned int counter_id = 1; counter_id <= 274; counter_id++){
            sprintf(name, "Counter ID %.3i", counter_id);
            cout << "Fitting TAGH " << name << endl;
            TH1I *rf_hist = Get1DHistogram("HLDetectorTiming", "TAGH_TDC_RF_Compare", name, false);
            if (rf_hist != NULL){
                // do the fit and store the result
                // Some fancy footwork to fit this periodic function
                TF1 *fa;
                if (rf_hist->GetBinCenter(rf_hist->GetMaximumBin()) < 0.0) fa = fitLeft;
                else fa = fitRight;
                fa->SetParLimits(1, -1*period / 2, period / 2);
                fa->SetParameter(1,rf_hist->GetBinCenter(rf_hist->GetMaximumBin()));
                fa->SetParLimits(2, 0.175, 1.2);
                fa->SetParameter(2, 0.4);
                TFitResultPtr r = rf_hist->Fit(fa, "S", "", -1*period / 2, period / 2);
                Int_t status = r;
                if ( status == 0){ // Fit succeeded
                    Double_t offset  = r->Parameter(1);
                    tdcRFOffsetTAGH[counter_id - 1] = offset;
                }
            }
        }
        outFile.open(prefix + "tagh_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file

        //int nbins = meanHist->GetNbinsX();
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();

        TH1D * selectedTAGHOffset = new TH1D("selectedTAGHOffset", "Selected TAGH Offset; Column; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        TH1I * TAGHOffsetDistribution = new TH1I("TAGHOffsetDistribution", "TAGH ADC Offset; ADC Offset [ns]; Entries", 1000, -500, 500);

        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 2; //ns (Full Width)
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
            selectedTAGHOffset->SetBinContent(i, -1 * maxMean);
            if(maxMean != 0) TAGHOffsetDistribution->Fill(-1 * maxMean);
        }

        double meanOffset = TAGHOffsetDistribution->GetMean();
        outFile.open(prefix + "tagh_adc_timing_offsets.txt");
        for (int i = 1; i <= nBinsX; i++){
            double mean = selectedTAGHOffset->GetBinContent(i);
            if (mean == 0) mean = meanOffset;
            outFile << i << " " << mean - meanOffset + tdcRFOffsetTAGH[i-1] << endl;
        }
        outFile.close();

        outFile.open(prefix + "tagh_tdc_timing_offsets.txt");
        for (int i = 1; i <= nBinsX; i++){
            outFile << i << " " << tdcRFOffsetTAGH[i-1] << endl;
        }
        outFile.close();

        outFile.open(prefix + "tagh_base_time.txt", ios::out);
        outFile << -1 * meanOffset << " 0.0" << endl;
        outFile.close();
    }

    // Now for the BCAL. First do the upstream, then the downstream. 
    // The entries in the CCDB are interleaved, so we need to store them temporarily.
    // A histogram will do.
    // We will use the averaging method and not FitSlices

    outFile.open(prefix + "bcal_tdc_timing_offsets.txt", ios::out | ios::trunc);
    outFile.close(); // clear file

    TH1D * selectedBCALOffset = new TH1D("selectedBCALOffset", "Selected BCAL TDC Offset; Column; Offset [ns]", 1152, 0.5, 1152 + 0.5);
    TH1I * BCALOffsetDistribution = new TH1I("BCALOffsetDistribution", "BCAL TDC Offset; TDC Offset [ns]; Entries", 500, -250, 250); 

    thisHist = Get2DHistogram("HLDetectorTiming", "BCAL", "BCALHit Upstream Per Channel TDC-ADC Hit Time");
    if(thisHist != NULL){
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 2; //ns (Full Width)
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
            selectedBCALOffset->SetBinContent(2*i - 1, maxMean);
            BCALOffsetDistribution->Fill(maxMean);
        }
    }

    thisHist = Get2DHistogram("HLDetectorTiming", "BCAL", "BCALHit Downstream Per Channel TDC-ADC Hit Time");
    if(thisHist != NULL){
        int nBinsX = thisHist->GetNbinsX();
        int nBinsY = thisHist->GetNbinsY();
        for (int i = 1 ; i <= nBinsX; i++){
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over the histogram
            float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            float timeWindow = 2; //ns (Full Width)
            int binWindow = int(timeWindow / nsPerBin);
            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++){
                int minBin = j;
                int maxBin = (j + binWindow) <= projY->GetNbinsX() ? (j + binWindow) : projY->GetNbinsX(); // Width is about 10ns (1ns per bin)
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

            selectedBCALOffset->SetBinContent(2*i, maxMean);
            BCALOffsetDistribution->Fill(maxMean);
        }
    }

    // We want to split things up so the bulk of the offset goes into the base time offset
    double meanTDCOffset = BCALOffsetDistribution->GetMean();

    // One is added, the other subtracted @_@
    outFile.open(prefix + "bcal_base_time.txt", ios::out);
    outFile << "0.0 " << -1*meanTDCOffset << endl;
    outFile.close();

    for (int i = 1 ; i <= selectedBCALOffset->GetNbinsX(); i++){
        double valueToUse = selectedBCALOffset->GetBinContent(i);
        if (valueToUse == 0) valueToUse = meanTDCOffset;
        outFile.open(prefix + "bcal_tdc_timing_offsets.txt", ios::out | ios::app);
        outFile << valueToUse - meanTDCOffset << endl;
        outFile.close();
    }

    thisFile->Write();
    return;
    }
