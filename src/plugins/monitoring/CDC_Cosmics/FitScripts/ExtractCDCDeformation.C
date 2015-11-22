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
void ExtractCDCDeformation(TString filename = "hd_root.root"){

    // Open our input and output file
    thisFile = TFile::Open(filename);
    TFile *outputFile = TFile::Open("CDCDeformation_Results.root", "RECREATE");

    // Check to make sure it is open
    if (thisFile == 0) {
        cout << "Unable to open file " << fileName.Data() << "...Exiting" << endl;
        return;
    }

    // This stream will be for outputting the results in a format suitable for the CCDB
    // Will wait to open until needed
    //ofstream textFile;
    //textFile.open("TimewalkBCAL.txt");

    // We want to display the direction of the shift as well as the magnitude in the "CDC view"
    // Let's make it happen
    int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};
    int Nstraws[28] = {42, 42, 54, 54, 66, 66, 80, 80, 93, 93, 106, 106, 123, 123, 135, 135, 146, 146, 158, 158, 170, 170, 182, 182, 197, 197, 209, 209};
    double radius[28] = {10.72134, 12.08024, 13.7795, 15.14602, 18.71726, 20.2438, 22.01672, 23.50008, 25.15616, 26.61158, 28.33624, 29.77388, 31.3817, 32.75838, 34.43478, 35.81146, 38.28542, 39.7002, 41.31564, 42.73042, 44.34078, 45.75302, 47.36084, 48.77054, 50.37582, 51.76012, 53.36286, 54.74716};
    double phi[28] = {0, 0.074707844, 0.038166294, 0.096247609, 0.05966371, 0.012001551, 0.040721951, 0.001334527, 0.014963808, 0.048683644, 0.002092645, 0.031681749, 0.040719354, 0.015197341, 0.006786058, 0.030005892, 0.019704045, -0.001782064, -0.001306618, 0.018592421, 0.003686784, 0.022132975, 0.019600866, 0.002343723, 0.021301449, 0.005348855, 0.005997358, 0.021018761};

    TH2D * Amplitude_view[29];
    TH2D * Direction_view[29];

    for(unsigned int iring=0; iring<28; iring++){
        double r_start = radius[iring] - 0.8;
        double r_end = radius[iring] + 0.8;
        double phi_start = phi[iring]; 
        double phi_end = phi_start + TMath::TwoPi();

        char hname[256];
        sprintf(hname, "Amplitude_view_ring[%d]", iring+1);
        Amplitude_view[iring+1] = new TH2D(hname, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
        sprintf(hname, "Direction_view_ring[%d]", iring+1);
        Direction_view[iring+1] = new TH2D(hname, "", Nstraws[iring], phi_start, phi_end, 1, r_start, r_end);
    }

    //Fit function for 
    TF1 *f1 = new TF1("f1", "[0] + [1] * TMath::Sin(x + [2])", -3.14, 3.14);
    f1->SetParLimits(0, 0.5, 1.0);
    f1->SetParLimits(1, 0.0, 0.35);
    //f1->SetParLimits(2, -3.14, 3.14);
    f1->SetParameters(0.78, 0.0, 0.0);

    // Make some histograms to get the distributions of the fit parameters
    TH1I *h1_c0 = new TH1I("h1_c0", "Distribution of Constant", 100, 0.5, 1.0);
    TH1I *h1_c1 = new TH1I("h1_c1", "Distribution of Amplitude", 100, 0.0, 0.35);
    TH1I *h1_c2 = new TH1I("h1_c2", "Distribution of Phase", 100, -10, 10);
    TH2I *h2_c0_c1 = new TH2I("h2_c0_c1", "c_{1} Vs. c_{0}; c_{0}; c_{1}", 100, 0.5, 1.0, 100, 0, 0.35);
    TH2I *h2_c0_c2 = new TH2I("h2_c0_c2", "c_{2} Vs. c_{0}; c_{0}; c_{2}", 100, 0.5, 1.0, 100, -10, 10);
    TH2I *h2_c1_c2 = new TH2I("h2_c1_c2", "c_{2} Vs. c_{1}; c_{1}; c_{2}", 100, 0.0, 0.35, 100, -10, 10);


    // Now we want to loop through all available module/layer/sector and try to make a fit of each one
    int ring = 1, straw = 1;
    while (ring <= 28){
        cout << "Entering Fit " << endl;
        char folder[100];
        sprintf(folder, "Ring %.2i", ring);
        char strawname[100];
        sprintf(strawname,"Straw %.3i Predicted Drift Distance Vs phi_DOCA", straw);
        TH2I *thisStrawHistogram = Get2DHistogram("CDC_Cosmic_Per_Straw",folder,strawname);

        if (thisStrawHistogram != NULL) {

            // Now to do our fits. This time we know there are 16 bins.
            double percentile95[16], percentile97[16], percentile99[16]; // Location of 95, 97,and 99th percentile bins
            double binCenter[16];
            char name[100];
            sprintf(name,"Ring %.2i Straw %.3i", ring, straw);
            TH1D *extractedPoints = new TH1D(name, name, 16, -3.14, 3.14);
            for (int i = 1; i <= thisStrawHistogram->GetNbinsX() ; i++){
                TH1D *projY = thisStrawHistogram->ProjectionY(" ", i, i);
                binCenter[i-1] = thisStrawHistogram->GetXaxis()->GetBinCenter(i);
                int nbins = projY->GetNbinsX();
                //Get the total nubmer of entries
                int nEntries = projY->GetEntries();
                if (nEntries == 0) continue;
                double errorFraction = TMath::Sqrt(nEntries) / nEntries;
                double perc95 = 0.95*nEntries, perc97 = 0.97 * nEntries, perc99 = 0.99 * nEntries;
                //Accumulate from the beginning to get total, mark 95, 97, 99% location
                int total = 0;
                for (int j = 0; j <= nbins; j++){
                    total += projY->GetBinContent(j);
                    if (total > perc99) percentile99[i-1] = projY->GetBinCenter(j);
                    else if (total > perc97) {
                        percentile97[i-1] = projY->GetBinCenter(j);
                        extractedPoints->SetBinContent(i, projY->GetBinCenter(j));
                        extractedPoints->SetBinError(i, errorFraction * projY->GetBinCenter(j));
                    }
                    else if (total > perc95) percentile95[i-1] = projY->GetBinCenter(j);
                }
            }
            // Now make a graph with the values found
            //TGraph *graph = new TGraph(16, binCenter, percentile97);
            //char name[100];
            //sprintf(name,"Ring %.2i Straw %.3i", ring, straw);
            //graph->SetName(name);
            //graph->SetTitle(name); 
            //graph->Write();
            //outputFile->Add(graph);
            //TFitResultPtr fr = graph->Fit(f1, "SR+");
            f1->SetParameters(0.78, 0.0, 0.0);
            TFitResultPtr fr = extractedPoints->Fit(f1, "SR");
            Int_t fitStatus = fr;
            if (fitStatus == 0){
                double c0 = fr->Parameter(0);
                double c1 = fr->Parameter(1);
                double c2 = fr->Parameter(2);
                // Move c2 to fit on our range
                while (c2 > TMath::Pi()) c2 -= 2 * TMath::Pi();
                while (c2 < -1* TMath::Pi()) c2 += 2 * TMath::Pi();
                h1_c0->Fill(c0); h1_c1->Fill(c1); h1_c2->Fill(c2);
                h2_c0_c1->Fill(c0,c1); h2_c0_c2->Fill(c0,c2); h2_c1_c2->Fill(c1,c2);
                Amplitude_view[ring]->SetBinContent(straw, c1);
                Direction_view[ring]->SetBinContent(straw, c2 + TMath::Pi() / 2); // Maximum of sin 
            }
            else {
                cout << "WARNING: Fit Status "<< fitStatus << " for ring " << ring << " straw " << straw << endl;
            }

        }

        // On to the next one
        straw++;
        if(straw > Nstraws[ring-1]){
            straw = 1;
            ring++;
        } 
    }
    // Draw CDC view plots
    TCanvas *c_Amplitude = new TCanvas("c_Amplitude","c_Amplitude", 800, 800);
    TH2D *axes = new TH2D("axes", "Amplitude of Sinusoid", 1, -65.0, 65.0, 1, -65.0, 65.0);
    axes->SetBinContent(1,1,-1.0/0.0); // Don't draw background color
    double minScale = 0.0, maxScale = 0.3;
    axes->SetStats(0);
    axes->Fill(100,100); // without this, the color ramp is not drawn
    axes->GetZaxis()->SetRangeUser(minScale, maxScale);
    axes->Draw("colz");
    for(unsigned int iring=1; iring<=28; iring++){
        Amplitude_view[iring]->GetZaxis()->SetRangeUser(minScale, maxScale);
        Amplitude_view[iring]->SetStats(0);
        Amplitude_view[iring]->Draw("same col pol");  // draw remaining histos without overwriting color palette
    }

    TCanvas *c_Amplitude = new TCanvas("c_Amplitude","c_Amplitude", 800, 800);
    TH2D *axes2 = new TH2D("axes", "Amplitude of Sinusoid", 1, -65.0, 65.0, 1, -65.0, 65.0);
    axes2->SetBinContent(1,1,-1.0/0.0); // Don't draw background color
    minScale = -3.14; maxScale = 3.14;
    axes2->SetStats(0);
    axes2->Fill(100,100); // without this, the color ramp is not drawn
    axes2->GetZaxis()->SetRangeUser(minScale, maxScale);

    axes2->Draw("colz");
    for(unsigned int iring=1; iring<=28; iring++){
        Direction_view[iring]->GetZaxis()->SetRangeUser(minScale, maxScale);
        Direction_view[iring]->SetStats(0);
        Direction_view[iring]->Draw("same col pol");  // draw remaining histos without overwriting color palette
    }

    cout << "Closing Files..." << endl;
    outputFile->Write();
    thisFile->Close();
    return;
}



