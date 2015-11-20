// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDeltaXVsVerticalPaddle
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TOFTrackDeltaYVsHorizontalPaddle
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TrackTOF2DPaddles_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPoint/TrackTOF2DPaddles_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPaddle/TrackYVsVerticalPaddle_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPaddle/TrackYVsVerticalPaddle_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPaddle/HorizontalPaddleVsTrackX_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOFPaddle/HorizontalPaddleVsTrackX_NoHit
{
	double locMinNumCountsForRatio = 20.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("TimeBased/TOFPaddle");
	TH2I* locHist_TrackYVsVerticalPaddle_HasHit = (TH2I*)gDirectory->Get("TrackYVsVerticalPaddle_HasHit");
	TH2I* locHist_TrackYVsVerticalPaddle_NoHit = (TH2I*)gDirectory->Get("TrackYVsVerticalPaddle_NoHit");
	TH2I* locHist_HorizontalPaddleVsTrackX_HasHit = (TH2I*)gDirectory->Get("HorizontalPaddleVsTrackX_HasHit");
	TH2I* locHist_HorizontalPaddleVsTrackX_NoHit = (TH2I*)gDirectory->Get("HorizontalPaddleVsTrackX_NoHit");

	gDirectory->cd("../TOFPoint");
	TH2I* locHist_TrackDeltaXVsVerticalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaXVsVerticalPaddle");
	TH2I* locHist_TrackDeltaYVsHorizontalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaYVsHorizontalPaddle");
	TH2I* locHist_PVsTheta_HasHit = (TH2I*)gDirectory->Get("PVsTheta_HasHit");
	TH2I* locHist_PVsTheta_NoHit = (TH2I*)gDirectory->Get("PVsTheta_NoHit");
	TH2I* locHist_TrackTOF2DPaddles_HasHit = (TH2I*)gDirectory->Get("TrackTOF2DPaddles_HasHit");
	TH2I* locHist_TrackTOF2DPaddles_NoHit = (TH2I*)gDirectory->Get("TrackTOF2DPaddles_NoHit");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_p1", "Matching_TOF", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackDeltaXVsVerticalPaddle != NULL)
	{
		locHist_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaXVsVerticalPaddle->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_TrackYVsVerticalPaddle_HasHit != NULL) && (locHist_TrackYVsVerticalPaddle_NoHit != NULL))
	{
		locHist_TrackYVsVerticalPaddle_HasHit->Rebin2D(1, 5); //44x260 -> 44x52
		locHist_TrackYVsVerticalPaddle_NoHit->Rebin2D(1, 5); //44x260 -> 44x52

		TH2I* locFoundHist = locHist_TrackYVsVerticalPaddle_HasHit;
		TH2I* locMissingHist = locHist_TrackYVsVerticalPaddle_NoHit;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / TOF Paddle Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

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
        locAcceptanceHist->GetZaxis()->SetRangeUser(0.6, 1.);
		locAcceptanceHist->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_HorizontalPaddleVsTrackX_HasHit != NULL) && (locHist_HorizontalPaddleVsTrackX_NoHit != NULL))
	{
		locHist_HorizontalPaddleVsTrackX_HasHit->Rebin2D(5, 1); //260x44 -> 52x44
		locHist_HorizontalPaddleVsTrackX_NoHit->Rebin2D(5, 1); //260x44 -> 52x44

		TH2I* locFoundHist = locHist_HorizontalPaddleVsTrackX_HasHit;
		TH2I* locMissingHist = locHist_HorizontalPaddleVsTrackX_NoHit;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / TOF Paddle Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

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
        locAcceptanceHist->GetZaxis()->SetRangeUser(0.6, 1.);
		locAcceptanceHist->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TrackDeltaYVsHorizontalPaddle != NULL)
	{
		locHist_TrackDeltaYVsHorizontalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TrackDeltaYVsHorizontalPaddle->Draw("COLZ");
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_PVsTheta_HasHit != NULL) && (locHist_PVsTheta_NoHit != NULL))
	{
		locHist_PVsTheta_HasHit->Rebin2D(8, 5); //280x400 -> 35x80
		locHist_PVsTheta_NoHit->Rebin2D(8, 5); //280x400 -> 35x80

		TH2I* locFoundHist = locHist_PVsTheta_HasHit;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / TOF Point Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

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
        locAcceptanceHist->GetZaxis()->SetRangeUser(0.6, 1.);
		locAcceptanceHist->Draw("COLZ");
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_TrackTOF2DPaddles_HasHit != NULL) && (locHist_TrackTOF2DPaddles_NoHit != NULL))
	{
		TH2I* locFoundHist = locHist_TrackTOF2DPaddles_HasHit;
		TH2I* locMissingHist = locHist_TrackTOF2DPaddles_NoHit;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / TOF Point Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

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
        locAcceptanceHist->GetZaxis()->SetRangeUser(0.6, 1.);
		locAcceptanceHist->Draw("COLZ");
	}
}
