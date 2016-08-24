int RFMacro_ROCTITimes(void)
{
	//Each tdc-tick of the TI-time corresponds to 4 ns, so max deviation from average allowed is 3.9
		//3.9 instead of 0: is average of all other ROCs. If a ROC is out of time, it will throw off the average when analyzing the other ROCs
	double locMaxDeviation = 3.9;

	gDirectory->cd("/"); //return to file base directory

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("ROCTIs");

	TList* locListOfKeys = gDirectory->GetListOfKeys();
	TKey* locPreviousKey = NULL;

	// loop over all keys in this directory
	for(int loc_i = 0; loc_i < locListOfKeys->GetEntries(); ++loc_i)
	{
		TKey* locKey = (TKey*)locListOfKeys->At(loc_i);
		if(locKey == NULL)
			continue;
		string locKeyName = locKey->GetName();

		//use only the first cycle number for each key
		if(locPreviousKey != NULL)
		{
			string locPreviousKeyName = locPreviousKey->GetName();
			if(locPreviousKeyName == locKeyName)
				continue;
		}
		locPreviousKey = locKey;

		TH1I* locHist = (TH1I*)gDirectory->Get(locKeyName.c_str());
		double locMean = locHist->GetMean();
		if(fabs(locMean) < locMaxDeviation)
			continue;

		//extract the roc # from the hist name
		size_t locROCIndex = locKeyName.rfind("ROC") + 3;
		string locROCNumString = locKeyName.substr(locROCIndex);
		cout << "ROC " << locROCNumString << " TI-time disagrees with the other ROCs." << endl;

		string locCanvasName = locKeyName + string("_Canvas");
		TCanvas* locCanvas = new TCanvas(locCanvasName.c_str(), locCanvasName.c_str(), 1200, 800);
		locHist->Draw();
	}

    return 1;
}
