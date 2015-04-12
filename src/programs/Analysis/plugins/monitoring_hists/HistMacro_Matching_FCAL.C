// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/FCAL/FCALTrackDistanceVsP
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/FCAL/FCALTrackDistanceVsTheta
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/FCAL/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/FCAL/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/FCAL/TrackFCALRowVsColumn_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/FCAL/TrackFCALRowVsColumn_NoHit

{
	double locMinNumCountsForRatio = 20.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("TimeBased/FCAL");
	TH2I* locHist_FCAL_TrackDistanceVsP = (TH2I*)gDirectory->Get("FCALTrackDistanceVsP");
	TH2I* locHist_FCAL_TrackDistanceVsTheta = (TH2I*)gDirectory->Get("FCALTrackDistanceVsTheta");
	TH2I* locHist_PVsTheta_HasHit_FCAL = (TH2I*)gDirectory->Get("PVsTheta_HasHit");
	TH2I* locHist_PVsTheta_NoHit_FCAL = (TH2I*)gDirectory->Get("PVsTheta_NoHit");
	TH2I* locHist_TrackFCALRowVsColumn_HasHit_FCAL = (TH2I*)gDirectory->Get("TrackFCALRowVsColumn_HasHit");
	TH2I* locHist_TrackFCALRowVsColumn_NoHit_FCAL = (TH2I*)gDirectory->Get("TrackFCALRowVsColumn_NoHit");

	//FCAL, by element (1 plot)
	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_FCAL", "Matching_FCAL", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCAL_TrackDistanceVsP != NULL)
	{
		locHist_FCAL_TrackDistanceVsP->Rebin2D(2, 2);
		locHist_FCAL_TrackDistanceVsP->GetXaxis()->SetTitleSize(0.05);
		locHist_FCAL_TrackDistanceVsP->GetYaxis()->SetTitleSize(0.05);
		locHist_FCAL_TrackDistanceVsP->GetXaxis()->SetLabelSize(0.05);
		locHist_FCAL_TrackDistanceVsP->GetYaxis()->SetLabelSize(0.05);
		locHist_FCAL_TrackDistanceVsP->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FCAL_TrackDistanceVsTheta != NULL)
	{
		locHist_FCAL_TrackDistanceVsTheta->Rebin2D(2, 2);
		locHist_FCAL_TrackDistanceVsTheta->GetXaxis()->SetTitleSize(0.05);
		locHist_FCAL_TrackDistanceVsTheta->GetYaxis()->SetTitleSize(0.05);
		locHist_FCAL_TrackDistanceVsTheta->GetXaxis()->SetLabelSize(0.05);
		locHist_FCAL_TrackDistanceVsTheta->GetYaxis()->SetLabelSize(0.05);
		locHist_FCAL_TrackDistanceVsTheta->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_PVsTheta_HasHit_FCAL != NULL) && (locHist_PVsTheta_NoHit_FCAL != NULL))
	{
		locHist_PVsTheta_HasHit_FCAL->Rebin2D(8, 5); //280x250 -> 35x50
		locHist_PVsTheta_NoHit_FCAL->Rebin2D(8, 5); //280x250 -> 35x50

		TH2I* locFoundHist = locHist_PVsTheta_HasHit_FCAL;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit_FCAL;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / FCAL Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());
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
	if((locHist_TrackFCALRowVsColumn_HasHit_FCAL != NULL) && (locHist_TrackFCALRowVsColumn_NoHit_FCAL != NULL))
	{
		TH2I* locFoundHist = locHist_TrackFCALRowVsColumn_HasHit_FCAL;
		TH2I* locMissingHist = locHist_TrackFCALRowVsColumn_NoHit_FCAL;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / FCAL Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());
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
}

