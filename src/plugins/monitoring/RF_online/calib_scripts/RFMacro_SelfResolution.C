// hnamepath: /RF/DeltaT_RF_Itself/FDCRF_SelfDeltaT
// hnamepath: /RF/DeltaT_RF_Itself/TOFRF_SelfDeltaT
// hnamepath: /RF/DeltaT_RF_Itself/TAGHRF_SelfDeltaT
// hnamepath: /RF/DeltaT_RF_Itself/PSCRF_SelfDeltaT

int RFMacro_SelfResolution(void)
{
	gStyle->SetOptStat(1111);
	gDirectory->cd("/"); //return to file base directory

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("DeltaT_RF_Itself");
	TH1I* locHist_FDCRF_SelfDeltaT = (TH1I*)gDirectory->Get("FDCRF_SelfDeltaT");
	TH1I* locHist_TOFRF_SelfDeltaT = (TH1I*)gDirectory->Get("TOFRF_SelfDeltaT");
	TH1I* locHist_TAGHRF_SelfDeltaT = (TH1I*)gDirectory->Get("TAGHRF_SelfDeltaT");
	TH1I* locHist_PSCRF_SelfDeltaT = (TH1I*)gDirectory->Get("PSCRF_SelfDeltaT");

	//Time resolutions
	double locTimeResolutionSq_TOF = locHist_TOFRF_SelfDeltaT->GetStdDev() / sqrt(2.0);
	locTimeResolutionSq_TOF *= locTimeResolutionSq_TOF;
	double locTimeResolutionSq_FDC = locHist_FDCRF_SelfDeltaT->GetStdDev() / sqrt(2.0);
	locTimeResolutionSq_FDC *= locTimeResolutionSq_FDC;
	double locTimeResolutionSq_PSC = locHist_PSCRF_SelfDeltaT->GetStdDev() / sqrt(2.0);
	locTimeResolutionSq_PSC *= locTimeResolutionSq_PSC;
	double locTimeResolutionSq_TAGH = locHist_TAGHRF_SelfDeltaT->GetStdDev() / sqrt(2.0);
	locTimeResolutionSq_TAGH *= locTimeResolutionSq_TAGH;

	//Print Coarse offsets to screen:
	cout << "Time-resolution-squared for TOF/TAGH/PSC/FDC are: " << locTimeResolutionSq_TOF << ", " << locTimeResolutionSq_TAGH << ", " << locTimeResolutionSq_PSC << ", " << locTimeResolutionSq_FDC << endl;

	//Print Coarse offsets to file:
	ofstream locOutputFileStream;
	locOutputFileStream.open("rf_time_resolution_sq.txt");
	locOutputFileStream << std::setprecision(6) << locTimeResolutionSq_TOF << " " << locTimeResolutionSq_TAGH << " " << locTimeResolutionSq_PSC << " " << locTimeResolutionSq_FDC << endl;
	locOutputFileStream.close();

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_SelfResolution", "RF_SelfResolution", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(2, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TOFRF_SelfDeltaT != NULL)
	{
		TH1I* locHist = locHist_TOFRF_SelfDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_TAGHRF_SelfDeltaT != NULL)
	{
		TH1I* locHist = locHist_TAGHRF_SelfDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_PSCRF_SelfDeltaT != NULL)
	{
		TH1I* locHist = locHist_PSCRF_SelfDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_FDCRF_SelfDeltaT != NULL)
	{
		TH1I* locHist = locHist_FDCRF_SelfDeltaT;
		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

    return 1;
}
