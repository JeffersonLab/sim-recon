// hnamepath: /Independent/Hist_DetectorMatching/WireBased/SC/SCTrackDeltaPhiVsZ
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCTrackDeltaPhiVsZ
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddleVsZ_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddleVsZ_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_BarrelRegion_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_BarrelRegion_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_NoseRegion_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddle_NoseRegion_NoHit

{
	double locMinNumCountsForRatio = 50.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("WireBased/SC");
	TH2I* locHist_SC_TrackDeltaPhiVsZ_WireBased = (TH2I*)gDirectory->Get("SCTrackDeltaPhiVsZ");
	gDirectory->cd("../../TimeBased/SC");
	TH2I* locHist_SC_TrackDeltaPhiVsZ_TimeBased = (TH2I*)gDirectory->Get("SCTrackDeltaPhiVsZ");
	TH2I* locHist_PVsTheta_HasHit_SC = (TH2I*)gDirectory->Get("PVsTheta_HasHit");
	TH2I* locHist_PVsTheta_NoHit_SC = (TH2I*)gDirectory->Get("PVsTheta_NoHit");
	TH2I* locHist_SCPaddleVsZ_HasHit_SC = (TH2I*)gDirectory->Get("SCPaddleVsZ_HasHit");
	TH2I* locHist_SCPaddleVsZ_NoHit_SC = (TH2I*)gDirectory->Get("SCPaddleVsZ_NoHit");
	TH1I* locHist_SCPaddle_BarrelRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_HasHit");
	TH1I* locHist_SCPaddle_BarrelRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_NoHit");
	TH1I* locHist_SCPaddle_NoseRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_HasHit");
	TH1I* locHist_SCPaddle_NoseRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_NoHit");

	//Get original pad margins
	double locLeftPadMargin = gStyle->GetPadLeftMargin();
	double locRightPadMargin = gStyle->GetPadRightMargin();
	double locTopPadMargin = gStyle->GetPadTopMargin();
	double locBottomPadMargin = gStyle->GetPadBottomMargin();

	//Set new pad margins
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.15);
//	gStyle->SetPadTopMargin(locTopPadMargin);
//	gStyle->SetPadBottomMargin(locBottomPadMargin);

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
	if(locHist_SC_TrackDeltaPhiVsZ_WireBased != NULL)
	{
		locHist_SC_TrackDeltaPhiVsZ_WireBased->GetYaxis()->SetTitleOffset(1.3);
		locHist_SC_TrackDeltaPhiVsZ_WireBased->GetXaxis()->SetRangeUser(35.0, 105.0);
		locHist_SC_TrackDeltaPhiVsZ_WireBased->GetXaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_WireBased->GetYaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_WireBased->GetXaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_WireBased->GetYaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_WireBased->Draw("COLZ");
		TF1* locFunc_High = new TF1("SC_PhiCut_High_WB", "[0] + [1]*exp([2]*(x - [3]))", 30.0, 110.0);
		locFunc_High->SetParameters(10.0, 0.5, 0.1, 60.0);
		locFunc_High->Draw("SAME");
		TF1* locFunc_Low = new TF1("SC_PhiCut_Low_WB", "-1.0*([0] + [1]*exp([2]*(x - [3])))", 30.0, 110.0);
		locFunc_Low->SetParameters(10.0, 0.5, 0.1, 60.0);
		locFunc_Low->Draw("SAME");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SC_TrackDeltaPhiVsZ_TimeBased != NULL)
	{
		locHist_SC_TrackDeltaPhiVsZ_TimeBased->GetYaxis()->SetTitleOffset(1.3);
		locHist_SC_TrackDeltaPhiVsZ_TimeBased->GetXaxis()->SetRangeUser(35.0, 105.0);
		locHist_SC_TrackDeltaPhiVsZ_TimeBased->GetXaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_TimeBased->GetYaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_TimeBased->GetXaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_TimeBased->GetYaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsZ_TimeBased->Draw("COLZ");
		TF1* locFunc_High = new TF1("SC_PhiCut_High_TB", "[0] + [1]*exp([2]*(x - [3]))", 30.0, 110.0);
		locFunc_High->SetParameters(8.0, 0.5, 0.1, 60.0);
		locFunc_High->Draw("SAME");
		TF1* locFunc_Low = new TF1("SC_PhiCut_Low_TB", "-1.0*([0] + [1]*exp([2]*(x - [3])))", 30.0, 110.0);
		locFunc_Low->SetParameters(8.0, 0.5, 0.1, 60.0);
		locFunc_Low->Draw("SAME");
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
		locAcceptanceHist->GetYaxis()->SetTitleOffset(1.3);
		locAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		locAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		locAcceptanceHist->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_SCPaddleVsZ_HasHit_SC != NULL) && (locHist_SCPaddleVsZ_NoHit_SC != NULL))
	{
		locHist_SCPaddleVsZ_HasHit_SC->Rebin2D(4, 1); //280x30 -> 70x30
		locHist_SCPaddleVsZ_NoHit_SC->Rebin2D(4, 1); //280x30 -> 70x30

		TH2I* locFoundHist = locHist_SCPaddleVsZ_HasHit_SC;
		TH2I* locMissingHist = locHist_SCPaddleVsZ_NoHit_SC;
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
		locAcceptanceHist->GetYaxis()->SetTitleOffset(1.3);
		locAcceptanceHist->GetXaxis()->SetRangeUser(35.0, 105.0);
		locAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		locAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		locAcceptanceHist->Draw("COLZ");
	}
    
	locHist_SCPaddle_BarrelRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_HasHit");
	locHist_SCPaddle_BarrelRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_BarrelRegion_NoHit");
	locHist_SCPaddle_NoseRegion_HasHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_HasHit");
	locHist_SCPaddle_NoseRegion_NoHit = (TH1I*)gDirectory->Get("SCPaddle_NoseRegion_NoHit");
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
				double locNumFoundError = sqrt(locNumFound_Barrel*(1.0 - locAcceptance));

				double locAcceptanceError = locNumFoundError/locTotal_Barrel;
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
				double locNumFoundError = sqrt(locNumFound_Nose*(1.0 - locAcceptance));

				double locAcceptanceError = locNumFoundError/locTotal_Nose;
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
		TLegend *locLegend = new TLegend(0.75, 0.17, 0.99, 0.32); //botleft x/y, topright x/y
		locLegend->SetHeader("Legend");
		locLegend->AddEntry(locBarrelGraph, "Barrel", "EP");
		locLegend->AddEntry(locNoseGraph, "Nose/Bend", "EP");
		locLegend->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_SCPaddleVsZ_HasHit_SC != NULL) && (locHist_SCPaddleVsZ_NoHit_SC != NULL))
	{
		string locHistName = string(locHist_SCPaddleVsZ_HasHit_SC->GetName()) + string("_ProjX");
		TH1* locFoundHist = locHist_SCPaddleVsZ_HasHit_SC->ProjectionX(locHistName.c_str());

		locHistName = string(locHist_SCPaddleVsZ_NoHit_SC->GetName()) + string("_ProjX");
		TH1* locMissingHist = locHist_SCPaddleVsZ_NoHit_SC->ProjectionX(locHistName.c_str());

		locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / SC Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle());

		TH1D* locAcceptanceHist = new TH1D(locHistName.c_str(), locHistTitle.c_str(), locFoundHist->GetNbinsX(), locFoundHist->GetXaxis()->GetXmin(), locFoundHist->GetXaxis()->GetXmax());
		for(Int_t loc_m = 1; loc_m <= locFoundHist->GetNbinsX(); ++loc_m)
		{
			double locNumMissing = locMissingHist->GetBinContent(loc_m);
			double locNumFound = locFoundHist->GetBinContent(loc_m);
			double locTotal = locNumMissing + locNumFound;
			if(!(locTotal >= locMinNumCountsForRatio))
			{
				locAcceptanceHist->SetBinContent(loc_m, 0.0);
				locAcceptanceHist->SetBinError(loc_m, 0.0);
				continue;
			}

			double locAcceptance = locNumFound/locTotal;
			locAcceptanceHist->SetBinContent(loc_m, locAcceptance);
			double locNumFoundError = sqrt(locNumFound*(1.0 - locAcceptance));

			double locAcceptanceError = locNumFoundError/locTotal;
			locAcceptanceHist->SetBinError(loc_m, locAcceptanceError);
		}
		locAcceptanceHist->SetEntries(locMissingHist->GetEntries() + locFoundHist->GetEntries());
		locAcceptanceHist->SetStats(kFALSE);

		locAcceptanceHist->GetXaxis()->SetRangeUser(35.0, 105.0);
		locAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		locAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		locAcceptanceHist->Draw("E1");
	}

	//Reset original pad margins
	gStyle->SetPadLeftMargin(locLeftPadMargin);
	gStyle->SetPadRightMargin(locRightPadMargin);
	gStyle->SetPadTopMargin(locTopPadMargin);
	gStyle->SetPadBottomMargin(locBottomPadMargin);
}

