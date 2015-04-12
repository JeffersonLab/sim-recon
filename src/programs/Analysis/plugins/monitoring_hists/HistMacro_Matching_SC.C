// hnamepath: /Independent/Hist_DetectorMatching/WireBased/SC/SCTrackDeltaPhiVsTheta
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCTrackDeltaPhiVsTheta
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddleVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/SC/SCPaddleVsTheta_NoHit

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

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_SC", "Matching_SC", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

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
}

