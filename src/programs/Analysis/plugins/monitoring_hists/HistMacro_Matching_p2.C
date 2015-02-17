// hnamepath: /Independent/Hist_DetectorStudies/Matching/BCAL_DeltaPhiVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matching/BCAL_DeltaZVsTheta
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_WireBased_GoodTrackFOM_HasHit_BCAL
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_WireBased_GoodTrackFOM_NoHit_BCAL
// hnamepath: /Independent/Hist_DetectorStudies/Matching/SC_TrackDeltaPhiVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_WireBased_GoodTrackFOM_HasHit_ST
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_WireBased_GoodTrackFOM_NoHit_ST

{
	double locMinNumCountsForRatio = 5.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Matching");
	TH2I* locHist_BCAL_DeltaPhiVsP = (TH2I*)gDirectory->Get("BCAL_DeltaPhiVsP");
	TH2I* locHist_BCAL_DeltaZVsTheta = (TH2I*)gDirectory->Get("BCAL_DeltaZVsTheta");
	TH2I* locHist_PVsTheta_HasHit_BCAL = (TH2I*)gDirectory->Get("PVsTheta_WireBased_GoodTrackFOM_HasHit_BCAL");
	TH2I* locHist_PVsTheta_NoHit_BCAL = (TH2I*)gDirectory->Get("PVsTheta_WireBased_GoodTrackFOM_NoHit_BCAL");

	TH2I* locHist_SC_TrackDeltaPhiVsP = (TH2I*)gDirectory->Get("SC_TrackDeltaPhiVsP");
	TH2I* locHist_PVsTheta_HasHit_SC = (TH2I*)gDirectory->Get("PVsTheta_WireBased_GoodTrackFOM_HasHit_ST");
	TH2I* locHist_PVsTheta_NoHit_SC = (TH2I*)gDirectory->Get("PVsTheta_WireBased_GoodTrackFOM_NoHit_ST");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_p2", "Matching_p2", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCAL_DeltaPhiVsP != NULL)
	{
		locHist_BCAL_DeltaPhiVsP->Rebin2D(2, 2);
		locHist_BCAL_DeltaPhiVsP->GetXaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaPhiVsP->GetYaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaPhiVsP->GetXaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaPhiVsP->GetYaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaPhiVsP->Draw("COLZ");
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCAL_DeltaZVsTheta != NULL)
	{
		locHist_BCAL_DeltaZVsTheta->Rebin2D(2, 2);
		locHist_BCAL_DeltaZVsTheta->GetXaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaZVsTheta->GetYaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaZVsTheta->GetXaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaZVsTheta->GetYaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaZVsTheta->Draw("COLZ");
		gPad->SetLogz();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_PVsTheta_HasHit_BCAL != NULL) && (locHist_PVsTheta_NoHit_BCAL != NULL))
	{
		locHist_PVsTheta_HasHit_BCAL->Rebin2D(4, 10); //280x400 -> 70x40
		locHist_PVsTheta_NoHit_BCAL->Rebin2D(4, 10); //280x400 -> 70x40

		TH2I* locFoundHist = locHist_PVsTheta_HasHit_BCAL;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit_BCAL;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		TH2D* locAcceptanceHist = new TH2D(locHistName.c_str(), "Track / BCAL Match Rate;#theta#circ;p (GeV/c)", locFoundHist->GetNbinsX(), locFoundHist->GetXaxis()->GetXmin(), locFoundHist->GetXaxis()->GetXmax(), locFoundHist->GetNbinsY(), locFoundHist->GetYaxis()->GetXmin(), locFoundHist->GetYaxis()->GetXmax());
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

		locAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		locAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		locAcceptanceHist->Draw("COLZ");
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_SC_TrackDeltaPhiVsP != NULL)
	{
		locHist_SC_TrackDeltaPhiVsP->Rebin2D(2, 2);
		locHist_SC_TrackDeltaPhiVsP->GetXaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsP->GetYaxis()->SetTitleSize(0.05);
		locHist_SC_TrackDeltaPhiVsP->GetXaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsP->GetYaxis()->SetLabelSize(0.05);
		locHist_SC_TrackDeltaPhiVsP->Draw("COLZ");
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_PVsTheta_HasHit_SC != NULL) && (locHist_PVsTheta_NoHit_SC != NULL))
	{
		locHist_PVsTheta_HasHit_SC->Rebin2D(4, 5); //280x400 -> 70x80
		locHist_PVsTheta_NoHit_SC->Rebin2D(4, 5); //280x400 -> 70x80

		TH2I* locFoundHist = locHist_PVsTheta_HasHit_SC;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit_SC;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		TH2D* locAcceptanceHist = new TH2D(locHistName.c_str(), "Track / SC Match Rate;#theta#circ;p (GeV/c)", locFoundHist->GetNbinsX(), locFoundHist->GetXaxis()->GetXmin(), locFoundHist->GetXaxis()->GetXmax(), locFoundHist->GetNbinsY(), locFoundHist->GetYaxis()->GetXmin(), locFoundHist->GetYaxis()->GetXmax());
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

		locAcceptanceHist->GetXaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetYaxis()->SetTitleSize(0.05);
		locAcceptanceHist->GetXaxis()->SetLabelSize(0.05);
		locAcceptanceHist->GetYaxis()->SetLabelSize(0.05);
		locAcceptanceHist->Draw("COLZ");
	}
}

