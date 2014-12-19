// hnamepath: /Independent/Hist_DetectorStudies/Matching/FCAL_TrackDistanceVsP
// hnamepath: /Independent/Hist_DetectorStudies/Matching/FCAL_TrackDistanceVsTheta
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_TimeBased_GoodTrackFOM_HasHit_FCAL
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_TimeBased_GoodTrackFOM_NoHit_FCAL
// hnamepath: /Independent/Hist_DetectorStudies/Matching/TOF_TrackDeltaXVsVerticalPaddle
// hnamepath: /Independent/Hist_DetectorStudies/Matching/TOF_TrackDeltaYVsHorizontalPaddle
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_TimeBased_GoodTrackFOM_HasHit_TOF
// hnamepath: /Independent/Hist_DetectorStudies/Matching/PVsTheta_TimeBased_GoodTrackFOM_NoHit_TOF

{
	double locMinNumCountsForRatio = 5.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorStudies");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("Matching");
	TH2I* locHist_FCAL_TrackDistanceVsP = (TH2I*)gDirectory->Get("FCAL_TrackDistanceVsP");
	TH2I* locHist_FCAL_TrackDistanceVsTheta = (TH2I*)gDirectory->Get("FCAL_TrackDistanceVsTheta");
	TH2I* locHist_PVsTheta_HasHit_FCAL = (TH2I*)gDirectory->Get("PVsTheta_TimeBased_GoodTrackFOM_HasHit_FCAL");
	TH2I* locHist_PVsTheta_NoHit_FCAL = (TH2I*)gDirectory->Get("PVsTheta_TimeBased_GoodTrackFOM_NoHit_FCAL");

	TH2I* locHist_TOF_TrackDeltaXVsVerticalPaddle = (TH2I*)gDirectory->Get("TOF_TrackDeltaXVsVerticalPaddle");
	TH2I* locHist_TOF_TrackDeltaYVsHorizontalPaddle = (TH2I*)gDirectory->Get("TOF_TrackDeltaYVsHorizontalPaddle");
	TH2I* locHist_PVsTheta_HasHit_TOF = (TH2I*)gDirectory->Get("PVsTheta_TimeBased_GoodTrackFOM_HasHit_TOF");
	TH2I* locHist_PVsTheta_NoHit_TOF = (TH2I*)gDirectory->Get("PVsTheta_TimeBased_GoodTrackFOM_NoHit_TOF");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_p1", "Matching_p1", 1200, 800);
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

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
		locHist_PVsTheta_HasHit_FCAL->Rebin2D(8, 5); //280x400 -> 35x80
		locHist_PVsTheta_NoHit_FCAL->Rebin2D(8, 5); //280x400 -> 35x80

		TH2I* locFoundHist = locHist_PVsTheta_HasHit_FCAL;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit_FCAL;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		TH2D* locAcceptanceHist = new TH2D(locHistName.c_str(), "Track / FCAL Match Rate;#theta#circ;p (GeV/c)", locFoundHist->GetNbinsX(), locFoundHist->GetXaxis()->GetXmin(), locFoundHist->GetXaxis()->GetXmax(), locFoundHist->GetNbinsY(), locFoundHist->GetYaxis()->GetXmin(), locFoundHist->GetYaxis()->GetXmax());
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
	if(locHist_TOF_TrackDeltaXVsVerticalPaddle != NULL)
	{
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetTitleSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetTitleSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetXaxis()->SetLabelSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->GetYaxis()->SetLabelSize(0.05);
		locHist_TOF_TrackDeltaXVsVerticalPaddle->Draw("COLZ");
	}

	locCanvas->cd(5);
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

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_PVsTheta_HasHit_TOF != NULL) && (locHist_PVsTheta_NoHit_TOF != NULL))
	{
		locHist_PVsTheta_HasHit_TOF->Rebin2D(8, 5); //280x400 -> 35x80
		locHist_PVsTheta_NoHit_TOF->Rebin2D(8, 5); //280x400 -> 35x80

		TH2I* locFoundHist = locHist_PVsTheta_HasHit_TOF;
		TH2I* locMissingHist = locHist_PVsTheta_NoHit_TOF;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		TH2D* locAcceptanceHist = new TH2D(locHistName.c_str(), "Track / TOF Match Rate;#theta#circ;p (GeV/c)", locFoundHist->GetNbinsX(), locFoundHist->GetXaxis()->GetXmin(), locFoundHist->GetXaxis()->GetXmax(), locFoundHist->GetNbinsY(), locFoundHist->GetYaxis()->GetXmin(), locFoundHist->GetYaxis()->GetXmax());
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

