// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/BCALDeltaPhiVsP
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/BCALDeltaZVsZ
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/PVsTheta_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/PVsTheta_NoHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/TrackBCALModuleVsZ_HasHit
// hnamepath: /Independent/Hist_DetectorMatching/TimeBased/BCAL/TrackBCALModuleVsZ_NoHit

{
	double locMinNumCountsForRatio = 50.0;

	//Goto Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("Hist_DetectorMatching");
	if(!locDirectory)
		return;
	locDirectory->cd();

	//Get Histograms
	gDirectory->cd("TimeBased/BCAL");
	TH2I* locHist_BCAL_DeltaPhiVsP = (TH2I*)gDirectory->Get("BCALDeltaPhiVsP");
	TH2I* locHist_BCAL_DeltaZVsZ = (TH2I*)gDirectory->Get("BCALDeltaZVsZ");
	TH2I* locHist_PVsTheta_HasHit_BCAL = (TH2I*)gDirectory->Get("PVsTheta_HasHit");
	TH2I* locHist_PVsTheta_NoHit_BCAL = (TH2I*)gDirectory->Get("PVsTheta_NoHit");
	TH2I* locHist_TrackBCALModuleVsZ_HasHit_BCAL = (TH2I*)gDirectory->Get("TrackBCALModuleVsZ_HasHit");
	TH2I* locHist_TrackBCALModuleVsZ_NoHit_BCAL = (TH2I*)gDirectory->Get("TrackBCALModuleVsZ_NoHit");

	//Get original pad margins
	double locLeftPadMargin = gStyle->GetPadLeftMargin();
	double locRightPadMargin = gStyle->GetPadRightMargin();
	double locTopPadMargin = gStyle->GetPadTopMargin();
	double locBottomPadMargin = gStyle->GetPadBottomMargin();

	//Set new pad margins
	gStyle->SetPadLeftMargin(0.15);
	gStyle->SetPadRightMargin(0.12);
//	gStyle->SetPadTopMargin(locTopPadMargin);
//	gStyle->SetPadBottomMargin(locBottomPadMargin);

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
		locHist_BCAL_DeltaPhiVsP->GetYaxis()->SetTitleOffset(1.3);
		locHist_BCAL_DeltaPhiVsP->GetXaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaPhiVsP->GetYaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaPhiVsP->GetXaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaPhiVsP->GetYaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaPhiVsP->Draw("COLZ");
		TF1* locFunc_High = new TF1("BCAL_PhiCut_High", "[0] + [1]*exp(-1.0*[2]*x)", 0.0, 4.0);
		locFunc_High->SetParameters(3.0, 12.0, 0.8);
		locFunc_High->Draw("SAME");
		TF1* locFunc_Low = new TF1("BCAL_PhiCut_Low", "-1.0*([0] + [1]*exp(-1.0*[2]*x))", 0.0, 4.0);
		locFunc_Low->SetParameters(3.0, 12.0, 0.8);
		locFunc_Low->Draw("SAME");
		gPad->SetLogz();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_BCAL_DeltaZVsZ != NULL)
	{
		locHist_BCAL_DeltaZVsZ->Rebin2D(2, 2);
		locHist_BCAL_DeltaZVsZ->GetYaxis()->SetTitleOffset(1.3);
		locHist_BCAL_DeltaZVsZ->GetXaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaZVsZ->GetYaxis()->SetTitleSize(0.05);
		locHist_BCAL_DeltaZVsZ->GetXaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaZVsZ->GetYaxis()->SetLabelSize(0.05);
		locHist_BCAL_DeltaZVsZ->Draw("COLZ");
		TF1* locFunc_High = new TF1("BCAL_ZCut_High", "30.0", 0.0, 450.0);
		locFunc_High->Draw("SAME");
		TF1* locFunc_Low = new TF1("BCAL_ZCut_Low", "-30.0", 0.0, 450.0);
		locFunc_Low->Draw("SAME");
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
		locAcceptanceHist->GetYaxis()->SetTitleOffset(1.3);
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
		TH1D* locHist_TrackBCALModule_HasHit = locHist_TrackBCALModuleVsZ_HasHit_BCAL->ProjectionY("BCALModule_HasHit");
		TH1D* locHist_TrackBCALModule_NoHit = locHist_TrackBCALModuleVsZ_NoHit_BCAL->ProjectionY("BCALModule_NoHit");

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

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if((locHist_TrackBCALModuleVsZ_HasHit_BCAL != NULL) && (locHist_TrackBCALModuleVsZ_NoHit_BCAL != NULL))
	{
		string locHistName = string(locHist_TrackBCALModuleVsZ_HasHit_BCAL->GetName()) + string("_ProjX");
		TH1* locFoundHist = locHist_TrackBCALModuleVsZ_HasHit_BCAL->ProjectionX(locHistName.c_str());

		locHistName = string(locHist_TrackBCALModuleVsZ_NoHit_BCAL->GetName()) + string("_ProjX");
		TH1* locMissingHist = locHist_TrackBCALModuleVsZ_NoHit_BCAL->ProjectionX(locHistName.c_str());

		locHistName = string(locFoundHist->GetName()) + string("_Acceptance");
		string locHistTitle = string("Track / BCAL Match Rate;") + string(locFoundHist->GetXaxis()->GetTitle());
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

