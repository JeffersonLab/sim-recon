TFile * thisFile;

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

void GetCCDBConstants(TString path, Int_t run, TString variation, vector<double>& container, Int_t column = 1){
    char command[256];
    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return 0;
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
void GetCCDBConstants(TString path, Int_t run, TString variation, double& constant1, double& constant2 = NULL){
    char command[256];
    sprintf(command, "ccdb dump %s:%i:%s", path.Data(), run, variation.Data());
    FILE* inputPipe = gSystem->OpenPipe(command, "r");
    if(inputPipe == NULL)
        return 0;
    //get the first (comment) line
    char buff[1024];
    if(fgets(buff, sizeof(buff), inputPipe) == NULL)
        return 0;
    //get the line containing the values
    while(fgets(buff, sizeof(buff), inputPipe) != NULL){
        istringstream locConstantsStream(buff);
        if(constant2 != NULL){
            locConstantsStream >> constant1 >> constant2;
        }
        else{
            locConstantsStream >> constant1;
        }
    }
    //Close the pipe
    gSystem->ClosePipe(inputPipe);
}

void ExtractFCALTiming(TString fileName = "hd_root.root", int runNumber = 2931, TString variation = "default", TString prefix = "")
{
    // set "prefix" in case you want to think about thowing the text files into another directory
    cout << "Performing TDC/ADC timing fits for File: " << fileName.Data() << " Run: " << runNumber << " Variation: " << variation.Data() << endl;

    // Open input ROOT file
    thisFile = TFile::Open( fileName , "UPDATE");
    if (thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }

    // Open output ROOT file
    TFile *outputROOTFile = new TFile("FCAL_Timing_Results.root","RECREATE");
    ofstream outFile;

    // We need to grab the existing values from the CCDB so that we know what was used in the creation of the ROOT file
    cout << "Grabbing CCDB constants..." << endl;
    double fcal_base_time;
    GetCCDBConstants("/FCAL/base_time_offset",runNumber, variation, fcal_base_time);
    vector<double> fcal_adc_offsets;
    GetCCDBConstants("/FCAL/timing_offsets" ,runNumber, variation, fcal_adc_offsets);

    cout << "Performing the fits..." << endl;
    TH2I *thisHistogram = Get2DHistogram("FCAL_Timing", "", "Target Time minus RF Time Vs. Channel Number");
    if (thisHistogram == NULL) return;

    int nBinsX = thisHistogram->GetNbinsX();
    TH1D * selectedFCALOffset = new TH1D("selectedFCALOffset", "Selected FCAL Offset; Column; Offset [ns]", nBinsX, -0.5, nBinsX - 0.5);

    for (int i = 1 ; i <= nBinsX; i++){ 
        TH1D *projY = thisHistogram->ProjectionY("temp", i, i);
        // Scan over the histogram
        //chose the correct number of bins based on the histogram
        //float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
        float timeWindow = 5; //ns (Full window)

        //Fit a Gaussian around the maximum bin 
        double maxBinLocation = projY->GetBinCenter(projY->GetMaximumBin());
        TFitResultPtr r = projY->Fit("gaus", "S", "", maxBinLocation - timeWindow / 2, maxBinLocation + timeWindow / 2);

        Int_t status = r;
        if ( status == 0){ // Fit succeeded
            double mean = r->Parameter(1);

            selectedFCALOffset->SetBinContent(i,mean);
        }
    }

    // Run through the second method for low statistics. If the values are found to be much different than the fit value, prefer these values...
    for (int i = 1 ; i <= nBinsX; i++){
        TH1D *projY = thisHistogram->ProjectionY("temp", i, i);
        // Scan over the histogram
        //chose the correct number of bins based on the histogram
        float nsPerBin = (projY->GetBinCenter(projY->GetNbinsX()) - projY->GetBinCenter(1)) / projY->GetNbinsX();
        float timeWindow = 5; //ns (Full window)

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
        double fitValue = selectedFCALOffset->GetBinContent(i);
        float threshold = 10; // If difference is greater than this number, take the low stats number
        if (fabs(maxMean - fitValue) > 3){
            selectedFCALOffset->SetBinContent(i,maxMean);
        }
    }

    // Add the old values to the new values, get the mean, and write the output
    vector<double> newValues;
    double sum = 0.0;
    for (int i = 1 ; i <= nBinsX; i++){
        double newVal = selectedFCALOffset->GetBinContent(i) + fcal_adc_offsets[i-1];
        sum += newVal;
        newValues.push_back(newVal);
    }
    // Average to go to base times
    double average = sum / nBinsX;
    cout << "average = " << average << endl;

    outFile.open("fcal_base_time.txt");
    outFile << (fcal_base_time - average) << endl;
    outFile.close();

    outFile.open("fcal_adc_offsets.txt");
    for (int i = 0 ; i < nBinsX; i++){
        outFile << (newValues[i] - average) << endl;
    }
    outFile.close();

    return;
}
