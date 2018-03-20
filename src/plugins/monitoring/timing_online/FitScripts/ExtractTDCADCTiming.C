namespace ExtractTDCADCTimingNS {
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
    thisFile->GetObject(fullName, histogram);
    if (histogram == 0){
        cout << "Unable to find histogram " << fullName.Data() << endl;
        return NULL;
    }
    return histogram;
}
};

void GetCCDBConstants(TString path, Int_t run, TString variation, vector<double>& container, Int_t column = 1){
    char command[256];
    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return;
    //get the remaining lines
    double entry;
    int counter = 0;
    while(fgets(buff, sizeof(buff), inputPipe) != NULL){
        istringstream locConstantsStream(buff);
        while (locConstantsStream >> entry){
            counter++;
            if (counter % column == 0) container.push_back(entry);
        }
    }
    //Close the pipe
    gSystem->ClosePipe(inputPipe);
}
//Overload this function to handle the base time offsets
void GetCCDBConstants1(TString path, Int_t run, TString variation, double& constant1){
    char command[256];
    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return;
    //get the line containing the values
    while(fgets(buff, sizeof(buff), inputPipe) != NULL){
        istringstream locConstantsStream(buff);
	locConstantsStream >> constant1;
    }
    //Close the pipe
    gSystem->ClosePipe(inputPipe);
}

void GetCCDBConstants2(TString path, Int_t run, TString variation, double& constant1, double& constant2){
    char command[256];
    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return;
    //get the line containing the values
    while(fgets(buff, sizeof(buff), inputPipe) != NULL){
        istringstream locConstantsStream(buff);
	locConstantsStream >> constant1 >> constant2;
    }
    //Close the pipe
    gSystem->ClosePipe(inputPipe);
}

int GetCCDBIndexTAGM(int column, int row){
    int CCDBIndex = column + row;
    if (column > 9) CCDBIndex += 5;
    if (column > 27) CCDBIndex += 5;
    if (column > 81) CCDBIndex += 5;
    if (column > 99) CCDBIndex += 5;

    return CCDBIndex;
}

int GetF1TDCslotTAGH(int id) {
    double N = 32.0; // channels per slot
    if (id >= 132 && id <= 172) throw("TAGH: unknown id in [132,172]");
    int HVid = (id <= 131) ? id : (id - 274 + 233);
    return int((HVid-1)/N) + 1;
}

Double_t FitFunctionLeft(Double_t *x, Double_t *par)
{
    Float_t xx =x[0];
    Double_t f = par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1]) / par[2] , 2) ) )+ par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1] - par[3]) / par[2] , 2) ) );
    return f;
}

Double_t FitFunctionRight(Double_t *x, Double_t *par)
{
    Float_t xx =x[0];
    Double_t f = par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1]) / par[2] , 2) ) )+ par[0]*TMath::Exp(-0.5 * (TMath::Power((xx - par[1] + par[3]) / par[2] , 2) ) );
    return f;
}

void ExtractTDCADCTiming(TString fileName = "hd_root.root", int runNumber = 10390, TString variation = "default", TString prefix = ""){

    // set "prefix" in case you want to think about thowing the text files into another directory
    cout << "Performing TDC/ADC timing fits for File: " << fileName.Data() << " Run: " << runNumber << " Variation: " << variation.Data() << endl;

    ExtractTDCADCTimingNS::thisFile = TFile::Open( fileName , "UPDATE");
    if (ExtractTDCADCTimingNS::thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }
    ofstream outFile;

    // We need to grab the existing values from the CCDB so that we know what was used in the creation of the ROOT file
    cout << "Grabbing CCDB constants..." << endl;

    // First are the base times for all of the detectors
    double cdc_base_time,      fcal_base_time;
    double sc_base_time_tdc,   sc_base_time_adc;
    double fdc_base_time_tdc,  fdc_base_time_adc;
    double bcal_base_time_tdc, bcal_base_time_adc;
    double tagm_base_time_tdc, tagm_base_time_adc;
    double tagh_base_time_tdc, tagh_base_time_adc;
    double tof_base_time_tdc,  tof_base_time_adc;

    double beam_period;

    GetCCDBConstants1("/CDC/base_time_offset" ,runNumber, variation, cdc_base_time);
    GetCCDBConstants1("/FCAL/base_time_offset",runNumber, variation, fcal_base_time);
    GetCCDBConstants1("/PHOTON_BEAM/RF/beam_period",runNumber, variation, beam_period);
    GetCCDBConstants2("/FDC/base_time_offset" ,runNumber, variation, fdc_base_time_adc, fdc_base_time_tdc);
    GetCCDBConstants2("/BCAL/base_time_offset" ,runNumber, variation, bcal_base_time_adc, bcal_base_time_tdc);
    GetCCDBConstants2("/PHOTON_BEAM/microscope/base_time_offset" ,runNumber, variation, tagm_base_time_adc, tagm_base_time_tdc);
    GetCCDBConstants2("/PHOTON_BEAM/hodoscope/base_time_offset" ,runNumber, variation, tagh_base_time_adc, tagh_base_time_tdc);
    GetCCDBConstants2("/START_COUNTER/base_time_offset" ,runNumber, variation, sc_base_time_adc, sc_base_time_tdc);
    GetCCDBConstants2("/TOF/base_time_offset" ,runNumber, variation, tof_base_time_adc, tof_base_time_tdc);

    cout << "CDC base times = " << cdc_base_time << endl;
    cout << "FCAL base times = " << fcal_base_time << endl;
    cout << "FDC base times = " << fdc_base_time_adc << ", " << fdc_base_time_tdc << endl;
    cout << "BCAL base times = " << bcal_base_time_adc << ", " << bcal_base_time_tdc << endl;
    cout << "SC base times = " << sc_base_time_adc << ", " << sc_base_time_tdc << endl;
    cout << "TOF base times = " << tof_base_time_adc << ", " << tof_base_time_tdc << endl;
    cout << "TAGH base times = " << tagh_base_time_adc << ", " << tagh_base_time_tdc << endl;
    cout << "TAGM base times = " << tagm_base_time_adc << ", " << tagm_base_time_tdc << endl;

    cout << "beam_period = " << beam_period << endl;

    // Then the channel by channel ADC and TDC times for those that need the calibration
    vector<double> bcal_tdc_offsets;
    vector<double> fcal_adc_offsets;
    vector<double> sc_adc_offsets, sc_tdc_offsets;
    vector<double> tagm_adc_offsets, tagm_tdc_offsets;
    vector<double> tagh_adc_offsets, tagh_tdc_offsets, tagh_counter_quality;
    vector<double> tof_adc_offsets, tof_tdc_offsets;
    GetCCDBConstants("/BCAL/TDC_offsets"    ,runNumber, variation, bcal_tdc_offsets);
    GetCCDBConstants("/FCAL/timing_offsets" ,runNumber, variation, fcal_adc_offsets);
    GetCCDBConstants("/START_COUNTER/adc_timing_offsets" ,runNumber, variation, sc_adc_offsets);
    GetCCDBConstants("/START_COUNTER/tdc_timing_offsets" ,runNumber, variation, sc_tdc_offsets);
    GetCCDBConstants("/PHOTON_BEAM/microscope/fadc_time_offsets" ,runNumber, variation, tagm_adc_offsets,3);// Interested in 3rd column
    GetCCDBConstants("/PHOTON_BEAM/microscope/tdc_time_offsets"  ,runNumber, variation, tagm_tdc_offsets,3);
    GetCCDBConstants("/PHOTON_BEAM/hodoscope/fadc_time_offsets"  ,runNumber, variation, tagh_adc_offsets,2);// Interested in 2nd column
    GetCCDBConstants("/PHOTON_BEAM/hodoscope/tdc_time_offsets"   ,runNumber, variation, tagh_tdc_offsets,2);
    GetCCDBConstants("/PHOTON_BEAM/hodoscope/counter_quality"    ,runNumber, variation, tagh_counter_quality,2);
    GetCCDBConstants("/TOF/adc_timing_offsets",runNumber, variation, tof_adc_offsets);
    GetCCDBConstants("/TOF/timing_offsets",runNumber, variation, tof_tdc_offsets);

    cout << "Done grabbing CCDB constants...Entering fits..." << endl;

    //Move the base times just for the tracking

    float CDC_ADC_Offset = 0.0;
    TH1I * this1DHist = ExtractTDCADCTimingNS::Get1DHistogram("HLDetectorTiming", "CDC", "CDCHit time");
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
        TFitResultPtr fr = this1DHist->Fit(f, "S", "", maximum - 20, maximum + 20); // Cant fix value at end of range
        double mean = fr->Parameter(1);
        float sigma = fr->Parameter(2);
        CDC_ADC_Offset = mean - sigma;
        delete f;
    }

    CDC_ADC_Offset *= -1;
    outFile.open(prefix + "cdc_base_time.txt");
    outFile << cdc_base_time + CDC_ADC_Offset << endl;
    outFile.close();

    float FDC_ADC_Offset = 0.0, FDC_TDC_Offset = 0.0;
    this1DHist = ExtractTDCADCTimingNS::Get1DHistogram("HLDetectorTiming", "FDC", "FDCHit Cathode time");
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
        TFitResultPtr fr = this1DHist->Fit(f, "S", "", maximum - 30, maximum + 20); // Cant fix value at end of range
        double mean = fr->Parameter(1);
        float sigma = fr->Parameter(2);
        FDC_ADC_Offset = mean;
        delete f;
    }
    this1DHist = ExtractTDCADCTimingNS::Get1DHistogram("HLDetectorTiming", "FDC", "FDCHit Wire time");
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
        TFitResultPtr fr = this1DHist->Fit(f, "S", "", maximum - 30, maximum + 20); // Cant fix value at end of range
        double mean = fr->Parameter(1);
        float sigma = fr->Parameter(2);
        FDC_TDC_Offset = mean;
        delete f;
    }

    FDC_ADC_Offset *= -1; FDC_TDC_Offset *= -1;
    outFile.open(prefix + "fdc_base_time.txt");
    outFile << fdc_base_time_adc + FDC_ADC_Offset << " " << fdc_base_time_tdc + FDC_TDC_Offset << endl;
    outFile.close();

    // Now that we have the file open, do all of the fits and write the output
    // Fit all plots with expected funtional form, output files for CCDB input
    // Do a finer alignement of the TDC and ADC's in the detectors that have both
    int minHits = 7;

    // In order to calibrate the SC in one step, we need to work with the base times and with the TDC/ADC offsets at the same time
    // Sort of complicates things but saves considerable time.

    float sc_tdc_base_time = 0.0;

    this1DHist = ExtractTDCADCTimingNS::Get1DHistogram("HLDetectorTiming", "SC", "SCHit TDC time");
    if(this1DHist != NULL){
        //Gaussian
        Double_t maximum = this1DHist->GetBinCenter(this1DHist->GetMaximumBin());
        TFitResultPtr fr = this1DHist->Fit("gaus", "S", "", maximum - 5, maximum + 5);
        double mean = fr->Parameter(1);
        sc_tdc_base_time = mean;
    }

    TH2I *thisHist = ExtractTDCADCTimingNS::Get2DHistogram("HLDetectorTiming", "SC", "SCHit TDC_ADC Difference");
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
            outFile << SC_ADC_Offset + sc_adc_offsets[iSC-1] - meanDiff << endl;
        }
        outFile.close();
    }

    outFile.open(prefix + "sc_base_time.txt", ios::out);
    outFile << sc_base_time_adc - (sc_tdc_base_time + meanDiff) << " " << sc_base_time_tdc - sc_tdc_base_time << endl;
    outFile.close();

    thisHist = ExtractTDCADCTimingNS::Get2DHistogram("HLDetectorTiming", "TOF", "TOFHit TDC_ADC Difference");
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
            outFile << tof_adc_offsets[i-1] + valueToUse - meanOffset << endl;
        }
        outFile.close();
        outFile.open(prefix + "tof_base_time.txt", ios::out);
        outFile << tof_base_time_adc - meanOffset << " " << tof_base_time_tdc << endl;
        outFile.close();

    }

    thisHist = ExtractTDCADCTimingNS::Get2DHistogram("HLDetectorTiming", "TAGM", "TAGMHit TDC_ADC Difference");
    if(thisHist != NULL){
        float tdcRFOffset[122] = {}; // 122 possible offsets
        char name[100];
        //First loop through the TDC RF times to see which ones we can correct
        TF1 * fitLeft = new TF1("fa",FitFunctionLeft, -1*beam_period / 2, beam_period / 2, 4);
        TF1 * fitRight = new TF1("fa",FitFunctionRight, -1*beam_period / 2, beam_period / 2, 4);
        for (unsigned int column = 1; column <= 102; column++){
            for (unsigned int row = 0; row <= 5; row++){
                sprintf(name, "Column %.3i Row %.1i", column, row);
                TH1I *rf_hist = ExtractTDCADCTimingNS::Get1DHistogram("HLDetectorTiming", "TAGM_TDC_RF_Compare", name, false);
                if (rf_hist != NULL){
                    // do the fit and store the result
                    // Some fancy footwork to fit this periodic function
                    cout << "Fitting TAGM " << name << endl;
                    TF1 *fa;
                    if (rf_hist->GetBinCenter(rf_hist->GetMaximumBin()) < 0.0) fa = fitLeft;
                    else fa = fitRight;
                    fa->SetParLimits(1, -1*beam_period / 2, beam_period / 2);
                    fa->SetParameter(1,rf_hist->GetBinCenter(rf_hist->GetMaximumBin()));
                    fa->SetParLimits(2, 0.175, 1.2);
                    fa->SetParameter(2, 0.4);
                    fa->FixParameter(3, beam_period);
                    TFitResultPtr r = rf_hist->Fit(fa, "S", "", -1*beam_period / 2, beam_period / 2);
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
            outFile << "0 " << column << " " << tagm_adc_offsets[index-1] + valueToUse - meanOffset + tdcRFOffset[index-1] << endl;
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    valueToUse = selectedTAGMOffset->GetBinContent(index);
                    double valueToUse = selectedTAGMOffset->GetBinContent(index);
                    if (valueToUse == 0) valueToUse = meanOffset; 
                    outFile << row << " " << column << " " << tagm_adc_offsets[index-1] + valueToUse - meanOffset + tdcRFOffset[index-1] << endl;
                }
            }
        }
        outFile.close();

        outFile.open(prefix + "tagm_tdc_timing_offsets.txt", ios::out);
        //for (int i = 1 ; i <= nBinsX; i++){
        // Loop over rows
        for (unsigned int column = 1; column <= 102; column++){
            int index = GetCCDBIndexTAGM(column, 0);
            outFile << "0 " << column << " " << tagm_tdc_offsets[index-1] + tdcRFOffset[index-1] << endl;
            if (column == 9 || column == 27 || column == 81 || column == 99){
                for (unsigned int row = 1; row <= 5; row++){
                    index = GetCCDBIndexTAGM(column, row);
                    outFile << row << " " << column << " " << tagm_tdc_offsets[index-1] + tdcRFOffset[index-1] << endl;
                }
            }
        }
        outFile.close();

        outFile.open(prefix + "tagm_base_time.txt", ios::out);
        outFile << tagm_base_time_adc - meanOffset << " " << tagm_base_time_tdc << endl;
        outFile.close();
    }

    thisHist = ExtractTDCADCTimingNS::Get2DHistogram("HLDetectorTiming", "TAGH", "TAGHHit TDC_ADC Difference");
    if (thisHist != NULL) {

        // Setup histograms for determining the most probable change in offset for each F1TDC slot
        // This is needed to account for the occasional uniform shift in offsets of the 32 counters in a slot
        const int NtdcSlots = 8;
        TH1I * tdcDist[NtdcSlots];
        TH1I * tdcadcDist[NtdcSlots];
        for (int i = 1; i <= NtdcSlots; i++) {
            stringstream ss; ss << i;
            TString s = ss.str();
            double range = 1000.0; double width = 0.1;
            int Nbins = range/width;
            double low = -0.5*range - 0.5*width;
            double high = 0.5*range - 0.5*width;
            tdcDist[i-1] = new TH1I("TAGHTDCOffsetDistribution_"+s, "TAGH TDC Offset (slot "+s+"); TDC Offset [ns]; Entries", Nbins, low, high);
            tdcadcDist[i-1] = new TH1I("TAGHTDCADCOffsetDistribution_"+s, "TAGH TDC/ADC Offset (slot "+s+"); TDC/ADC Offset [ns]; Entries", Nbins, low, high);
        }

        // First take a look at the RF offsets
        double tdcRFOffsetTAGH[274] = {}; // 274 possible offsets
        char name[100];
        // Loop through the TDC RF times to see which ones we can correct
        TF1 * fitLeft = new TF1("fa",FitFunctionLeft, -1*beam_period / 2, beam_period / 2, 4);
        TF1 * fitRight = new TF1("fa",FitFunctionRight, -1*beam_period / 2, beam_period / 2, 4);

        for (unsigned int counter_id = 1; counter_id <= 274; counter_id++) {
            sprintf(name, "Counter ID %.3i", counter_id);
            cout << "Fitting TAGH " << name << endl;
            TH1I *rf_hist = ExtractTDCADCTimingNS::Get1DHistogram("HLDetectorTiming", "TAGH_TDC_RF_Compare", name, false);
            if (rf_hist != NULL) {
                // do the fit and store the result
                // Some fancy footwork to fit this periodic function
                TF1 *fa;
                if (rf_hist->GetBinCenter(rf_hist->GetMaximumBin()) < 0.0) fa = fitLeft;
                else fa = fitRight;
                fa->SetParLimits(1, -1*beam_period / 2, beam_period / 2);
                fa->SetParameter(1,rf_hist->GetBinCenter(rf_hist->GetMaximumBin()));
                fa->SetParLimits(2, 0.175, 1.2);
                fa->SetParameter(2, 0.4);
                fa->FixParameter(3, beam_period);
                TFitResultPtr r = rf_hist->Fit(fa, "S", "", -1*beam_period / 2, beam_period / 2);
                Int_t status = r;
                if (status == 0) { // Fit succeeded
                    Double_t offset  = r->Parameter(1);
                    tdcRFOffsetTAGH[counter_id - 1] = offset;
                    if (tagh_counter_quality[counter_id-1] == 0.0) continue;
                    int tdc_slot = GetF1TDCslotTAGH(counter_id);
                    tdcDist[tdc_slot - 1]->Fill(offset);
                }
            }
        }
        // Most probable change in TDC/RF offset per F1TDC slot
        double mpDeltaTDC[NtdcSlots];
        for (int i = 1; i <= NtdcSlots; i++) {
            int mpBin = tdcDist[i-1]->GetMaximumBin();
            mpDeltaTDC[i-1] = (mpBin > 0) ? tdcDist[i-1]->GetBinCenter(mpBin) : 0.0;
        }

        outFile.open(prefix + "tagh_adc_timing_offsets.txt", ios::out | ios::trunc);
        outFile.close(); // clear file

        int nBinsX = thisHist->GetNbinsX();

        TH1D * selectedTAGHOffset = new TH1D("selectedTAGHOffset", "Selected TAGH Offset; Column; Offset [ns]", nBinsX, 0.5, nBinsX + 0.5);
        for (int i = 1 ; i <= nBinsX; i++) {
            TH1D *projY = thisHist->ProjectionY("temp", i, i);
            // Scan over histogram to find mean offset in timeWindow with largest integral
            double nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
            double timeWindow = 2.0; //ns (Full Width)
            int binWindow = int(timeWindow / nsPerBin);
            double maxEntries = 0;
            double maxMean = 0;
            for (int j = 1 ; j <= projY->GetNbinsX();j++) {
                int minBin = j;
                int maxBin = (j + binWindow) <= projY->GetNbinsX() ? (j + binWindow) : projY->GetNbinsX();
                double sum = 0, nEntries = 0;
                for (int bin = minBin; bin <= maxBin; bin++) {
                    sum += projY->GetBinContent(bin) * projY->GetBinCenter(bin);
                    nEntries += projY->GetBinContent(bin);
                    if (bin == maxBin) {
                        if (nEntries > maxEntries) {
                            maxMean = sum / nEntries;
                            maxEntries = nEntries;
                        }
                    }
                }
            }
            selectedTAGHOffset->SetBinContent(i, maxMean);
            if (tagh_counter_quality[i-1] == 0.0) continue;
            int tdc_slot = GetF1TDCslotTAGH(i);
            if (maxEntries != 0.0) tdcadcDist[tdc_slot - 1]->Fill(maxMean);
        }
        // Most probable change in TDC/ADC offset per F1TDC slot
        double mpDeltaADC[NtdcSlots];
        for (int i = 1; i <= NtdcSlots; i++) {
            int mpBin = tdcadcDist[i-1]->GetMaximumBin();
            double mpDeltaTDCADC = (mpBin > 0) ? tdcadcDist[i-1]->GetBinCenter(mpBin) : 0.0;
            mpDeltaADC[i-1] = mpDeltaTDC[i-1] - mpDeltaTDCADC;
        }

        outFile.open(prefix + "tagh_adc_timing_offsets.txt");
        double limit = 2.5; // ns
        double ccdb_sum = 0.0;
        for (int i = 1; i <= nBinsX; i++) ccdb_sum += tagh_adc_offsets[i-1];
        double c1_adcOffset = 0.0;
        for (int i = 1; i <= nBinsX; i++) {
            if (tagh_counter_quality[i-1] == 0.0) {
                outFile << i << " " << 0 << endl;
                continue;
            }
            int tdc_slot = GetF1TDCslotTAGH(i);
            double ccdb = tagh_adc_offsets[i-1];
            double delta = tdcRFOffsetTAGH[i-1] - selectedTAGHOffset->GetBinContent(i);
            if (ccdb_sum > 0.0 && fabs(delta - mpDeltaADC[tdc_slot-1]) > limit) {
                delta = mpDeltaADC[tdc_slot-1];
            }
            double offset = ccdb + delta;
            if (i == 1) c1_adcOffset = offset;
            offset -= c1_adcOffset;
            outFile << i << " " << offset << endl;
        }
        outFile.close();

        outFile.open(prefix + "tagh_tdc_timing_offsets.txt");
        ccdb_sum = 0.0;
        for (int i = 1; i <= nBinsX; i++) ccdb_sum += tagh_tdc_offsets[i-1];
        double c1_tdcOffset = 0.0;
        for (int i = 1; i <= nBinsX; i++) {
            if (tagh_counter_quality[i-1] == 0.0) {
                outFile << i << " " << 0 << endl;
                continue;
            }
            int tdc_slot = GetF1TDCslotTAGH(i);
            double ccdb = tagh_tdc_offsets[i-1];
            double delta = tdcRFOffsetTAGH[i-1];
            if (ccdb_sum > 0.0 && fabs(delta - mpDeltaTDC[tdc_slot-1]) > limit) {
                delta = mpDeltaTDC[tdc_slot-1];
            }
            double offset = ccdb + delta;
            if (i == 1) c1_tdcOffset = offset;
            offset -= c1_tdcOffset;
            outFile << i << " " << offset << endl;
        }
        outFile.close();

        outFile.open(prefix + "tagh_base_time.txt", ios::out);
        outFile << tagh_base_time_adc - c1_adcOffset << " " << tagh_base_time_tdc - c1_tdcOffset << endl;
        outFile.close();
    }

    // Now for the BCAL. First do the upstream, then the downstream. 
    // The entries in the CCDB are interleaved, so we need to store them temporarily.
    // A histogram will do.
    // We will use the averaging method and not FitSlices

    outFile.open(prefix + "bcal_tdc_timing_offsets.txt", ios::out | ios::trunc);
    outFile.close(); // clear file

    TH1D * selectedBCALOffset = new TH1D("selectedBCALOffset", "Selected BCAL TDC Offset; Column; Offset [ns]", 1152, 0.5, 1152 + 0.5);
    TH1I * BCALOffsetDistribution = new TH1I("BCALOffsetDistribution", "BCAL TDC Offset; TDC Offset [ns]; Entries", 500, -500, 500); 

    thisHist = ExtractTDCADCTimingNS::Get2DHistogram("HLDetectorTiming", "BCAL", "BCALHit Upstream Per Channel TDC-ADC Hit Time");
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
            BCALOffsetDistribution->Fill(maxMean+bcal_tdc_offsets[2*i -2]);
        }
    }

    thisHist = ExtractTDCADCTimingNS::Get2DHistogram("HLDetectorTiming", "BCAL", "BCALHit Downstream Per Channel TDC-ADC Hit Time");
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
            BCALOffsetDistribution->Fill(maxMean+bcal_tdc_offsets[2*i -1]);
        }
    }

    // We want to split things up so the bulk of the offset goes into the base time offset
    double meanTDCOffset = BCALOffsetDistribution->GetMean();

    // One is added, the other subtracted @_@
    outFile.open(prefix + "bcal_base_time.txt", ios::out);
    outFile << bcal_base_time_adc << " " << bcal_base_time_tdc - meanTDCOffset << endl;
    outFile.close();

    for (int i = 1 ; i <= selectedBCALOffset->GetNbinsX(); i++){
        double valueToUse = selectedBCALOffset->GetBinContent(i);
        if (valueToUse == 0) valueToUse = meanTDCOffset;
        outFile.open(prefix + "bcal_tdc_timing_offsets.txt", ios::out | ios::app);
        outFile << bcal_tdc_offsets[i-1] + valueToUse - meanTDCOffset << endl;
        outFile.close();
    }

    ExtractTDCADCTimingNS::thisFile->Write();
    return;
    }
