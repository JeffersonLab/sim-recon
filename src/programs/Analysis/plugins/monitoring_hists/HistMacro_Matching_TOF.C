// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOF/TOFTrackDeltaXVsVerticalPaddle
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOF/TOFTrackDeltaYVsHorizontalPaddle
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOF/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOF/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOF/TrackTOF2DPaddles_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/TOF/TrackTOF2DPaddles_NoHit

{
	double locMinNumCountsForRatio = 5.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("TimeBased/TOF");
	TH2I* locHist_TOF_TrackDeltaXVsVerticalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaXVsVerticalPaddle");
	TH2I* locHist_TOF_TrackDeltaYVsHorizontalPaddle = (TH2I*)gDirectory->Get("TOFTrackDeltaYVsHorizontalPaddle");
	TH2I* locHist_PVsTheta_HasHit_TOF = (TH2I*)gDirectory->Get("PVsTheta_HasHit");
	TH2I* locHist_PVsTheta_NoHit_TOF = (TH2I*)gDirectory->Get("PVsTheta_NoHit");
	TH2I* locHist_TrackTOF2DPaddles_HasHit_TOF = (TH2I*)gDirectory->Get("TrackTOF2DPaddles_HasHit");
	TH2I* locHist_TrackTOF2DPaddles_NoHit_TOF = (TH2I*)gDirectory->Get("TrackTOF2DPaddles_NoHit");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_p1", "Matching_TOF", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOF_TrackDeltaXVsVerticalPaddle != NULL)
	{
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOF_TrackDeltaYVsHorizontalPaddle != NULL)
	{
		locHist_TOF_TrackDeltaYVsHorizontalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TOF_TrackDeltaYVsHorizontalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TOF_TrackDeltaYVsHorizontalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TOF_TrackDeltaYVsHorizontalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TOF_TrackDeltaYVsHorizontalPaddle->Draw("COLZ");
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_PVsTheta_HasHit_TOF != NULL) && (locHist_PVsTheta_NoHit_TOF != NULL))
	{
		locHist_PVsTheta_HasHit_TOF->Rebin2D(8, 5); //280x400 -> 35x80
		locHist_PVsTheta_NoHit_TOF->Rebin2D(8, 5); //280x400 -> 35x80

		TH2I* locFoundHist = locHist_PVsTheta_HasHit_TOF;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit_TOF;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / TOF Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

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
	if((locHist_TrackTOF2DPaddles_HasHit_TOF != NULL) && (locHist_TrackTOF2DPaddles_NoHit_TOF != NULL))
	{
		TH2I* locFoundHist = locHist_TrackTOF2DPaddles_HasHit_TOF;
		TH2I* locMissingHist = locHist_TrackTOF2DPaddles_NoHit_TOF;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / TOF Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());

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
