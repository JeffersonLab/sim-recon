TObjArray* Get_RFDeltaTHists(string locSystem);
void Draw_Array(string locSystem, TObjArray* locTotalArray);

int RFMacro_InternalDeltaTs(void)
{
	gDirectory->cd("/"); //return to file base directory

	//Goto RF Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("DeltaT_All");
	TObjArray* locFDCArray = Get_RFDeltaTHists("FDC");
	TObjArray* locPSCArray = Get_RFDeltaTHists("PSC");
	TObjArray* locTAGHArray = Get_RFDeltaTHists("TAGH");
	TObjArray* locTOFArray = Get_RFDeltaTHists("TOF");

	//Draw
	Draw_Array("FDC", locFDCArray);
	Draw_Array("PSC", locPSCArray);
	Draw_Array("TAGH", locTAGHArray);
	Draw_Array("TOF", locTOFArray);
}

TObjArray* Get_RFDeltaTHists(string locSystem)
{
	int locFirstHit = 0, locSecondHit = 0;
	TObjArray* locTotalArray = new TObjArray();
	TObjArray* locRowArray = NULL;
	while(true)
	{
		ostringstream locHistName;
		locHistName << locSystem << "RF_DeltaT_" << locFirstHit << "_" << locSecondHit;
		TObject* locObject = gDirectory->Get(locHistName.str().c_str());
		if(locObject == NULL)
		{
			if(locSecondHit == locFirstHit)
				break;
			++locFirstHit;
			locSecondHit = locFirstHit;
			continue;
		}

		if(locSecondHit == locFirstHit)
		{
			locRowArray = new TObjArray();
			for(int loc_i = 0; loc_i < locSecondHit; ++loc_i)
				locRowArray->AddLast(NULL);
			locTotalArray->AddLast(locRowArray);
		}
		locRowArray->AddLast(locObject);
		++locSecondHit;
	}
	return locTotalArray;
}

void Draw_Array(string locSystem, TObjArray* locTotalArray)
{
	int locNumHits = locTotalArray->GetEntriesFast();

	//Means
	string locHistName = string("DeltaTMeans_") + locSystem + string("RF");
	string locHistTitle = string("#Deltat Means (ns), ") + locSystem + string(" RF;First Hit;Second Hit");
	TH2F* locMeanHist = new TH2F(locHistName.c_str(), locHistTitle.c_str(), locNumHits, 0.5, locNumHits + 0.5, locNumHits, 0.5, locNumHits + 0.5);

	//Std Dev
	locHistName = string("DeltaTSigmas_") + locSystem + string("RF");
	locHistTitle = string("#Deltat #sigma/#sqrt{2} (ns), ") + locSystem + string(" RF;First Hit;Second Hit");
	TH2F* locStdDevHist = new TH2F(locHistName.c_str(), locHistTitle.c_str(), locNumHits, 0.5, locNumHits + 0.5, locNumHits, 0.5, locNumHits + 0.5);

	//Skewness
	locHistName = string("DeltaTSkewness_") + locSystem + string("RF");
	locHistTitle = string("#Deltat Skewness, ") + locSystem + string(" RF;First Hit;Second Hit");
	TH2F* locSkewnessHist = new TH2F(locHistName.c_str(), locHistTitle.c_str(), locNumHits, 0.5, locNumHits + 0.5, locNumHits, 0.5, locNumHits + 0.5);

	//Kurtosis
	locHistName = string("DeltaTKurtosis_") + locSystem + string("RF");
	locHistTitle = string("#Deltat Kurtosis, ") + locSystem + string(" RF;First Hit;Second Hit");
	TH2F* locKurtosisHist = new TH2F(locHistName.c_str(), locHistTitle.c_str(), locNumHits, 0.5, locNumHits + 0.5, locNumHits, 0.5, locNumHits + 0.5);

	//Get values, build mega-canvas
	string locCanvasName = string("DeltaTs_") + locSystem + string("RFCanvas");
	int locPadHeight = 200;
	gStyle->SetOptStat(0000);
	gStyle->SetNumberContours(20);

	TCanvas* locCanvas = new TCanvas(locCanvasName.c_str(), locCanvasName.c_str(), locNumHits*locPadHeight*16.0/9.0, locNumHits*locPadHeight);
	locCanvas->Divide(locNumHits, locNumHits);
	for(int loc_i = 0; loc_i < locNumHits; ++loc_i) //first hit: columns
	{
		TObjArray* locRowArray = (TObjArray*)locTotalArray->At(loc_i);
		for(int loc_j = 0; loc_j < locNumHits; ++loc_j) //second hit: rows
		{
			TH1* locDeltaTHist = (TH1*)locRowArray->At(loc_j);
			if((locDeltaTHist == NULL) || (loc_i == loc_j))
			{
				locMeanHist->SetBinContent(loc_i + 1, loc_j + 1, -1.0/0.0);
				locStdDevHist->SetBinContent(loc_i + 1, loc_j + 1, -1.0/0.0);
				locSkewnessHist->SetBinContent(loc_i + 1, loc_j + 1, -1.0/0.0);
				locKurtosisHist->SetBinContent(loc_i + 1, loc_j + 1, -1.0/0.0);
				if(locDeltaTHist == NULL)
					continue;
			}

			double locPeakLocation = locDeltaTHist->GetBinCenter(locDeltaTHist->GetMaximumBin());
			locDeltaTHist->GetXaxis()->SetRangeUser(locPeakLocation - 10.0, locPeakLocation + 10.0);
			double locMean = locDeltaTHist->GetMean();
			if(loc_i != loc_j)
			{
				locMeanHist->SetBinContent(loc_i + 1, loc_j + 1, locMean);
				locStdDevHist->SetBinContent(loc_i + 1, loc_j + 1, locDeltaTHist->GetStdDev()/sqrt(2.0));
				locSkewnessHist->SetBinContent(loc_i + 1, loc_j + 1, locDeltaTHist->GetSkewness());
				locKurtosisHist->SetBinContent(loc_i + 1, loc_j + 1, locDeltaTHist->GetKurtosis());
			}

			locDeltaTHist->GetXaxis()->SetRangeUser(locMean - 0.5, locMean + 0.5);
			//i is
			int locCanvasIndex = (locNumHits - loc_j - 1)*locNumHits + loc_i + 1;
			locCanvas->cd(locCanvasIndex);
			locDeltaTHist->Draw();
		}
	}

	locCanvasName = string("DeltaTMeans_") + locSystem + string("RFCanvas");
	new TCanvas(locCanvasName.c_str(), locCanvasName.c_str(), 1600, 900);
	locMeanHist->Draw("COLZ");

	locCanvasName = string("DeltaTSigmas_") + locSystem + string("RFCanvas");
	new TCanvas(locCanvasName.c_str(), locCanvasName.c_str(), 1600, 900);
	locStdDevHist->Draw("COLZ");

	locCanvasName = string("DeltaTSkewness_") + locSystem + string("RFCanvas");
	new TCanvas(locCanvasName.c_str(), locCanvasName.c_str(), 1600, 900);
	locSkewnessHist->Draw("COLZ");

	locCanvasName = string("DeltaTKurtosis_") + locSystem + string("RFCanvas");
	new TCanvas(locCanvasName.c_str(), locCanvasName.c_str(), 1600, 900);
	locKurtosisHist->Draw("COLZ");

    return 1;
}
