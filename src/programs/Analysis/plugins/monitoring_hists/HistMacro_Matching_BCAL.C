// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/BCALDeltaPhiVsP
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/BCALDeltaZVsTheta
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/TrackBCALModuleVsZ_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/TrackBCALModuleVsZ_NoHit

{
	double locMinNumCountsForRatio = 20.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("TimeBased/BCAL");
	TH2I* locHist_BCAL_DeltaPhiVsP = (TH2I*)gDirectory->Get("BCALDeltaPhiVsP");
	TH2I* locHist_BCAL_DeltaZVsTheta = (TH2I*)gDirectory->Get("BCALDeltaZVsTheta");
	TH2I* locHist_PVsTheta_HasHit_BCAL = (TH2I*)gDirectory->Get("PVsTheta_HasHit");
	TH2I* locHist_PVsTheta_NoHit_BCAL = (TH2I*)gDirectory->Get("PVsTheta_NoHit");
	TH2I* locHist_TrackBCALModuleVsZ_HasHit_BCAL = (TH2I*)gDirectory->Get("TrackBCALModuleVsZ_HasHit");
	TH2I* locHist_TrackBCALModuleVsZ_NoHit_BCAL = (TH2I*)gDirectory->Get("TrackBCALModuleVsZ_NoHit");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("Matching_BCAL", "Matching_BCAL", 1200, 800);
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
		string locHistTitle = string("Track / BCAL Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());
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
	if((locHist_TrackBCALModuleVsZ_HasHit_BCAL != NULL) && (locHist_TrackBCALModuleVsZ_NoHit_BCAL != NULL))
	{
		locHist_TrackBCALModuleVsZ_HasHit_BCAL->Rebin2D(5, 1); //450x48 -> 90x48
		locHist_TrackBCALModuleVsZ_NoHit_BCAL->Rebin2D(5, 1); //450x48 -> 90x48

		TH2I* locFoundHist = locHist_TrackBCALModuleVsZ_HasHit_BCAL;
		TH2I* locMissingHist = locHist_TrackBCALModuleVsZ_NoHit_BCAL;
		string locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / BCAL Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle()) + string(";") + string(locFoundHist->GetYaxis()->GetTitle());
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

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_TrackBCALModuleVsZ_HasHit_BCAL != NULL) && (locHist_TrackBCALModuleVsZ_NoHit_BCAL != NULL))
	{
		TH1D* locHist_TrackBCALModule_HasHit = locHist_TrackBCALModuleVsZ_HasHit_BCAL->ProjectionY("SCPaddle_HasHit");
		TH1D* locHist_TrackBCALModule_NoHit = locHist_TrackBCALModuleVsZ_NoHit_BCAL->ProjectionY("SCPaddle_NoHit");

		TH1D* loc1DFoundHist = locHist_TrackBCALModule_HasHit;
		TH1D* loc1DMissingHist = locHist_TrackBCALModule_NoHit;
		string locHistName = string(loc1DFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / BCAL Match Rate;") + string(loc1DFoundHist->GetXaxis()->GetTitle()) + string(";") + string(loc1DFoundHist->GetYaxis()->GetTitle());

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

