// hnamepath: /Independent/Hist_DetectorMatching/WireBased/SC/SCTrackDeltaPhiVsTheta
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCTrackDeltaPhiVsTheta
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddleVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddleVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_BarrelRegion_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_BarrelRegion_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_NoseRegion_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_NoseRegion_NoHit

{
	double locMinNumCountsForRatio = 20.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("WireBased/SC");
	TH2I* locHist_SC_TrackDeltaPhiVsTheta_WireBased = (TH2I*)gDirectory->Get("SCTrackDeltaPhiVsTheta");
	gDirectory->cd("../../TimeBased/SC");
	TH2I* locHist_SC_TrackDeltaPhiVsTheta_TimeBased = (TH2I*)gDirectory->Get("SCTrackDeltaPhiVsTheta");
	TH2I* locHist_PVsTheta_HasHit_SC = (TH2I*)gDirectory->Get("PVsTheta_HasHit");
	TH2I* locHist_PVsTheta_NoHit_SC = (TH2I*)gDirectory->Get("PVsTheta_NoHit");
	TH2I* locHist_SCPaddleVsTheta_HasHit_SC = (TH2I*)gDirectory->Get("SCPaddleVsTheta_HasHit");
	TH2I* locHist_SCPaddleVsTheta_NoHit_SC = (TH2I*)gDirectory->Get("SCPaddleVsTheta_NoHit");
	TH1I* locHist_SCPaddle_BarrelRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_HasHit");
	TH1I* locHist_SCPaddle_BarrelRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_NoHit");
	TH1I* locHist_SCPaddle_NoseRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_HasHit");
	TH1I* locHist_SCPaddle_NoseRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_NoHit");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_SC", "Matching_SC", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SC_TrackDeltaPhiVsTheta_WireBased != NULL)
	{
		locHist_SC_TrackDeltaPhiVsTheta_WireBased->GetXaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_WireBased->GetYaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_WireBased->GetXaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_WireBased->GetYaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_WireBased->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SC_TrackDeltaPhiVsTheta_TimeBased != NULL)
	{
		locHist_SC_TrackDeltaPhiVsTheta_TimeBased->GetXaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_TimeBased->GetYaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_TimeBased->GetXaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_TimeBased->GetYaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsTheta_TimeBased->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_PVsTheta_HasHit_SC != NULL) && (locHist_PVsTheta_NoHit_SC != NULL))
	{
		locHist_PVsTheta_HasHit_SC->Rebin2D(4, 5); //280x250 -> 70x50
		locHist_PVsTheta_NoHit_SC->Rebin2D(4, 5); //280x250 -> 70x50

		TH2I* locFoundHist = locHist_PVsTheta_HasHit_SC;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit_SC;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / SC Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

		TH2D* locAcceptanceHist = new TH2D(locHistName.c_str(), locHistTitle.c_str(), locFoundHist->GetNbinsX(), locFoundHist->GetXaxis()->GetXmin(), locFoundHist->GetXaxis()->GetXmax(), locFoundHist->GetNbinsY(), locFoundHist->GetYaxis()->GetXmin(), locFoundHist->GetYaxis()->GetXmax());
		for(Int_t loc_m = 1; loc_m <= locFoundHist->GetNbinsX(); ++loc_m)
		{
			for(Int_t loc_j = 1; loc_j <= locFoundHist->GetNbinsY(); ++loc_j)
			{
				double locNumMissing = locMissingHist->GetBinContent(loc_m, loc_j);
				double locNumFound = locFoundHist->GetBinContent(loc_m, loc_j);
				double locTotal = locNumMissing + locNumFound;
				if(!(locTotal >= locMinNumCountsForRatio))
				{
					locAcceptanceHist->SetBinContent(loc_m, loc_j, 0.0);
					locAcceptanceHist->SetBinError(loc_m, loc_j, 1.0/0.0);
					continue;
				}

				double locAcceptance = locNumFound/locTotal;
				if(!(locAcceptance > 0.0))
					locAcceptance = 0.00001; //so that it shows up on the histogram
				locAcceptanceHist->SetBinContent(loc_m, loc_j, locAcceptance);
				double locNumFoundError = sqrt(locNumFound*(1.0 - locAcceptance));

				double locAcceptanceError = locNumFoundError/locTotal;
				locAcceptanceHist->SetBinError(loc_m, loc_j, locAcceptanceError);
			}
		}
		locAcceptanceHist->SetEntries(locMissingHist->GetEntries() + locFoundHist->GetEntries());
		locAcceptanceHist->SetStats(kFALSE);

		locAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		locAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		locAcceptanceHist->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_SCPaddleVsTheta_HasHit_SC != NULL) && (locHist_SCPaddleVsTheta_NoHit_SC != NULL))
	{
		locHist_SCPaddleVsTheta_HasHit_SC->Rebin2D(4, 1); //280x30 -> 70x30
		locHist_SCPaddleVsTheta_NoHit_SC->Rebin2D(4, 1); //280x30 -> 70x30

		TH2I* locFoundHist = locHist_SCPaddleVsTheta_HasHit_SC;
		TH2I* locMissingHist = locHist_SCPaddleVsTheta_NoHit_SC;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / SC Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

		TH2D* locAcceptanceHist = new TH2D(locHistName.c_str(), locHistTitle.c_str(), locFoundHist->GetNbinsX(), locFoundHist->GetXaxis()->GetXmin(), locFoundHist->GetXaxis()->GetXmax(), locFoundHist->GetNbinsY(), locFoundHist->GetYaxis()->GetXmin(), locFoundHist->GetYaxis()->GetXmax());
		for(Int_t loc_m = 1; loc_m <= locFoundHist->GetNbinsX(); ++loc_m)
		{
			for(Int_t loc_j = 1; loc_j <= locFoundHist->GetNbinsY(); ++loc_j)
			{
				double locNumMissing = locMissingHist->GetBinContent(loc_m, loc_j);
				double locNumFound = locFoundHist->GetBinContent(loc_m, loc_j);
				double locTotal = locNumMissing + locNumFound;
				if(!(locTotal >= locMinNumCountsForRatio))
				{
					locAcceptanceHist->SetBinContent(loc_m, loc_j, 0.0);
					locAcceptanceHist->SetBinError(loc_m, loc_j, 1.0/0.0);
					continue;
				}

				double locAcceptance = locNumFound/locTotal;
				if(!(locAcceptance > 0.0))
					locAcceptance = 0.00001; //so that it shows up on the histogram
				locAcceptanceHist->SetBinContent(loc_m, loc_j, locAcceptance);
				double locNumFoundError = sqrt(locNumFound*(1.0 - locAcceptance));

				double locAcceptanceError = locNumFoundError/locTotal;
				locAcceptanceHist->SetBinError(loc_m, loc_j, locAcceptanceError);
			}
		}
		locAcceptanceHist->SetEntries(locMissingHist->GetEntries() + locFoundHist->GetEntries());
		locAcceptanceHist->SetStats(kFALSE);

		locAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		locAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		locAcceptanceHist->Draw("COLZ");
	}

	TH1I* locHist_SCPaddle_BarrelRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_HasHit");
	TH1I* locHist_SCPaddle_BarrelRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_NoHit");
	TH1I* locHist_SCPaddle_NoseRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_HasHit");
	TH1I* locHist_SCPaddle_NoseRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_NoHit");
	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_SCPaddle_BarrelRegion_HasHit != NULL) && (locHist_SCPaddle_BarrelRegion_NoHit != NULL) && (locHist_SCPaddle_NoseRegion_HasHit != NULL) && (locHist_SCPaddle_NoseRegion_NoHit != NULL))
	{
		int locNumCounters = 30;
		vector<double> locNoseXAxisVector, locNoseEffVector, locNoseEffUncertaintyVector;
		vector<double> locBarrelXAxisVector, locBarrelEffVector, locBarrelEffUncertaintyVector;

		for(int loc_i = 0; loc_i < locNumCounters; ++loc_i)
		{
			//Barrel
			double locNumFound_Barrel = locHist_SCPaddle_BarrelRegion_HasHit->GetBinContent(loc_i + 1);
			double locNumMissing_Barrel = locHist_SCPaddle_BarrelRegion_NoHit->GetBinContent(loc_i + 1);
			double locTotal_Barrel = locNumMissing_Barrel + locNumFound_Barrel;
			if(locTotal_Barrel >= locMinNumCountsForRatio)
			{
				double locAcceptance = locNumFound_Barrel/locTotal_Barrel;
				locBarrelEffVector.push_back(locAcceptance);
				double locNumFoundError = sqrt(locNumFound*(1.0 - locAcceptance));

				double locAcceptanceError = locNumFoundError/locTotal;
				locBarrelEffUncertaintyVector.push_back(locAcceptanceError);

				locBarrelXAxisVector.push_back(loc_i + 1);
			}

			//Nose & Bend
			double locNumFound_Nose = locHist_SCPaddle_NoseRegion_HasHit->GetBinContent(loc_i + 1);
			double locNumMissing_Nose = locHist_SCPaddle_NoseRegion_NoHit->GetBinContent(loc_i + 1);
			double locTotal_Nose = locNumMissing_Nose + locNumFound_Nose;
			if(locTotal_Nose >= locMinNumCountsForRatio)
			{
				double locAcceptance = locNumFound_Nose/locTotal_Nose;
				locNoseEffVector.push_back(locAcceptance);
				double locNumFoundError = sqrt(locNumFound*(1.0 - locAcceptance));

				double locAcceptanceError = locNumFoundError/locTotal;
				locNoseEffUncertaintyVector.push_back(locAcceptanceError);

				locNoseXAxisVector.push_back(loc_i + 1);
			}
		}

		double* locNoseXAxisArray = new double[locNoseEffVector.size()];
		for(size_t loc_f = 0; loc_f < locNoseEffVector.size(); ++loc_f)
			locNoseXAxisArray[loc_f] = locNoseXAxisVector[loc_f];

		double* locNoseXUncertaintyArray = new double[locNoseEffVector.size()];
		for(size_t loc_f = 0; loc_f < locNoseEffVector.size(); ++loc_f)
			locNoseXUncertaintyArray[loc_f] = 0.0;

		double* locNoseEffArray = new double[locNoseEffVector.size()];
		for(size_t loc_f = 0; loc_f < locNoseEffVector.size(); ++loc_f)
			locNoseEffArray[loc_f] = locNoseEffVector[loc_f];

		double* locNoseEffUncertaintyArray = new double[locNoseEffVector.size()];
		for(size_t loc_f = 0; loc_f < locNoseEffVector.size(); ++loc_f)
			locNoseEffUncertaintyArray[loc_f] = locNoseEffUncertaintyVector[loc_f];



		double* locBarrelXAxisArray = new double[locBarrelEffVector.size()];
		for(size_t loc_f = 0; loc_f < locBarrelEffVector.size(); ++loc_f)
			locBarrelXAxisArray[loc_f] = locBarrelXAxisVector[loc_f];

		double* locBarrelXUncertaintyArray = new double[locBarrelEffVector.size()];
		for(size_t loc_f = 0; loc_f < locBarrelEffVector.size(); ++loc_f)
			locBarrelXUncertaintyArray[loc_f] = 0.0;

		double* locBarrelEffArray = new double[locBarrelEffVector.size()];
		for(size_t loc_f = 0; loc_f < locBarrelEffVector.size(); ++loc_f)
			locBarrelEffArray[loc_f] = locBarrelEffVector[loc_f];

		double* locBarrelEffUncertaintyArray = new double[locBarrelEffVector.size()];
		for(size_t loc_f = 0; loc_f < locBarrelEffVector.size(); ++loc_f)
			locBarrelEffUncertaintyArray[loc_f] = locBarrelEffUncertaintyVector[loc_f];

		TGraphErrors* locNoseGraph = new TGraphErrors(locNoseEffVector.size(), locNoseXAxisArray, locNoseEffArray, locNoseXUncertaintyArray, locNoseEffUncertaintyArray);
		locNoseGraph->SetLineColor(kBlack);
		locNoseGraph->SetMarkerColor(kBlack);
		locNoseGraph->SetMarkerStyle(20);
		locNoseGraph->SetMarkerSize(1);

		TGraphErrors* locBarrelGraph = new TGraphErrors(locBarrelEffVector.size(), locBarrelXAxisArray, locBarrelEffArray, locBarrelXUncertaintyArray, locBarrelEffUncertaintyArray);
		locBarrelGraph->SetLineColor(kBlue);
		locBarrelGraph->SetMarkerColor(kBlue);
		locBarrelGraph->SetMarkerStyle(22);
		locBarrelGraph->SetMarkerSize(1);

		string locGraphName = "SC_EffGraphs";
		string locGraphTitle = string("Track / SC Match Rate;") + string(locHist_SCPaddle_BarrelRegion_HasHit->GetXaxis()->GetTitle());
		TMultiGraph* locMultiGraph = new TMultiGraph(locGraphName.c_str(), locGraphTitle.c_str());
		locMultiGraph->Add(locNoseGraph);
		locMultiGraph->Add(locBarrelGraph);
		locMultiGraph->Draw("ap");

		//make legend!
		TLegend *locLegend = new TLegend(0.75, 0.88, 0.99, 0.999); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locBarrelGraph, "Barrel", "EP");
		locLegend->AddEntry(locNoseGraph, "Nose/Bend", "EP");
		locLegend->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_SCPaddleVsTheta_HasHit_SC != NULL) && (locHist_SCPaddleVsTheta_NoHit_SC != NULL))
	{
		TH1D* locHist_SCPaddle_HasHit = locHist_SCPaddleVsTheta_HasHit_SC->ProjectionY("SCPaddle_HasHit");
		TH1D* locHist_SCPaddle_NoHit = locHist_SCPaddleVsTheta_NoHit_SC->ProjectionY("SCPaddle_NoHit");

		TH1D* loc1DFoundHist = locHist_SCPaddle_HasHit;
		TH1D* loc1DMissingHist = locHist_SCPaddle_NoHit;
		string locHistName = string(loc1DFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / SC Match Rate;") + string(loc1DFoundHist->GetXaxis()->GetTitle()) + string(";") + string(loc1DFoundHist->GetYaxis()->GetTitle());

		TH1D* loc1DAcceptanceHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), loc1DFoundHist->GetNbinsX(), loc1DFoundHist->GetXaxis()->GetXmin(), loc1DFoundHist->GetXaxis()->GetXmax());
		for(Int_t loc_m = 1; loc_m <= loc1DFoundHist->GetNbinsX(); ++loc_m)
		{
			double locNumMissing = loc1DMissingHist->GetBinContent(loc_m);
			double locNumFound = loc1DFoundHist->GetBinContent(loc_m);
			double locTotal = locNumMissing + locNumFound;
			if(!(locTotal >= locMinNumCountsForRatio))
			{
				loc1DAcceptanceHist->SetBinContent(loc_m, 0.0);
				loc1DAcceptanceHist->SetBinError(loc_m, 1.0/0.0);
				continue;
			}

			double locAcceptance = locNumFound/locTotal;
			if(!(locAcceptance > 0.0))
				locAcceptance = 0.00001; //so that it shows up on the histogram
			loc1DAcceptanceHist->SetBinContent(loc_m, locAcceptance);
			double locNumFoundError = sqrt(locNumFound*(1.0 - locAcceptance));

			double locAcceptanceError = locNumFoundError/locTotal;
			loc1DAcceptanceHist->SetBinError(loc_m, locAcceptanceError);
		}
		loc1DAcceptanceHist->SetEntries(loc1DMissingHist->GetEntries() + loc1DFoundHist->GetEntries());
		loc1DAcceptanceHist->SetStats(kFALSE);

		loc1DAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		loc1DAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		loc1DAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		loc1DAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		loc1DAcceptanceHist->Draw("E1");
	}
}

