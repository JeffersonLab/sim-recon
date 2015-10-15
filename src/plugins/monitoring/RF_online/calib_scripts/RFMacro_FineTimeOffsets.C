// hnamepath: /RF/AverageDeltaT_RF_OtherRFs/RFDeltaT_FDC_TOF
// hnamepath: /RF/AverageDeltaT_RF_OtherRFs/RFDeltaT_FDC_TAGH
// hnamepath: /RF/AverageDeltaT_RF_OtherRFs/RFDeltaT_FDC_PSC
// hnamepath: /RF/AverageDeltaT_RF_OtherRFs/RFDeltaT_PSC_TAGH
// hnamepath: /RF/AverageDeltaT_RF_OtherRFs/RFDeltaT_PSC_TOF
// hnamepath: /RF/AverageDeltaT_RF_OtherRFs/RFDeltaT_TAGH_TOF

double gRFSignalPeriod = 1000.0/499.0;
int gRebinF1sAmount = 25;
int gRebinTOFAmount = 2;

//Commit constants with:
//ccdb add PHOTON_BEAM/RF/time_offset -r <run_min>-<run_max> rf_fine_time_offsets.txt #"fine time offsets"
//ccdb add PHOTON_BEAM/RF/time_offset_var -r <run_min>-<run_max> rf_time_offset_vars.txt #"time offset variances"

Double_t Periodic_Gaussian_Func(Double_t* locXArray, Double_t* locParamArray)
{
	Double_t locValue = locParamArray[0]*TMath::Gaus(locXArray[0], locParamArray[1], locParamArray[2]);
	if(locParamArray[1] < 0.0)
		locValue += locParamArray[0]*TMath::Gaus(locXArray[0], locParamArray[1] + gRFSignalPeriod, locParamArray[2]);
	else
		locValue += locParamArray[0]*TMath::Gaus(locXArray[0], locParamArray[1] - gRFSignalPeriod, locParamArray[2]);

	return locValue;
}

TF1* Create_FitFunc(TH1I* locHist)
{
	Int_t locMaxBin = locHist->GetMaximumBin();
	double locMean = locHist->GetBinCenter(locMaxBin);

	string locFuncName = string(locHist->GetName()) + string("_Func");
	TF1 *locFunc = new TF1(locFuncName.c_str(), Periodic_Gaussian_Func, -0.5*gRFSignalPeriod, 0.5*gRFSignalPeriod, 3);
	locFunc->SetParameters(locHist->GetBinContent(locMaxBin), locMean, 0.1);
	locFunc->SetParNames("Gaussian Height", "Gaussian #mu", "Gaussian #sigma");

	return locFunc;
}

int main(int locRunNumber)
{
	//INPUT RUN NUMBER MUST BE A RUN NUMBER IN THE RUN RANGE YOU ARE TRYING TO COMMIT CONSTANTS TO
	TDirectory *locTopDirectory = gDirectory;

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("AverageDeltaT_RF_OtherRFs");
	TH1I* locHist_RFDeltaT_FDC_TOF = (TH1I*)gDirectory->Get("RFDeltaT_FDC_TOF");
	TH1I* locHist_RFDeltaT_FDC_TAGH = (TH1I*)gDirectory->Get("RFDeltaT_FDC_TAGH");
	TH1I* locHist_RFDeltaT_FDC_PSC = (TH1I*)gDirectory->Get("RFDeltaT_FDC_PSC");
	TH1I* locHist_RFDeltaT_PSC_TAGH = (TH1I*)gDirectory->Get("RFDeltaT_PSC_TAGH");
	TH1I* locHist_RFDeltaT_PSC_TOF = (TH1I*)gDirectory->Get("RFDeltaT_PSC_TOF");
	TH1I* locHist_RFDeltaT_TAGH_TOF = (TH1I*)gDirectory->Get("RFDeltaT_TAGH_TOF");

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_FineTimeOffsets", "RF_FineTimeOffsets", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	TF1* locTOFFitFunc_TAGH = NULL;
	if(locHist_RFDeltaT_TAGH_TOF != NULL)
	{
		TH1I* locHist = locHist_RFDeltaT_TAGH_TOF;
		locHist->Rebin(gRebinTOFAmount);
		locTOFFitFunc_TAGH = Create_FitFunc(locHist);
		locHist->Fit(locTOFFitFunc_TAGH, "QR");

		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	TF1* locTOFFitFunc_PSC = NULL;
	if(locHist_RFDeltaT_PSC_TOF != NULL)
	{
		TH1I* locHist = locHist_RFDeltaT_PSC_TOF;
		locHist->Rebin(gRebinTOFAmount);
		locTOFFitFunc_PSC = Create_FitFunc(locHist);
		locHist->Fit(locTOFFitFunc_PSC, "QR");

		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	TF1* locTOFFitFunc_FDC = NULL;
	if(locHist_RFDeltaT_FDC_TOF != NULL)
	{
		TH1I* locHist = locHist_RFDeltaT_FDC_TOF;
		locTOFFitFunc_FDC = Create_FitFunc(locHist);
		locHist->Rebin(gRebinTOFAmount);
		locHist->Fit(locTOFFitFunc_FDC, "QR");

		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_TAGH != NULL)
	{
		TH1I* locHist = locHist_RFDeltaT_FDC_TAGH;
		locHist->Rebin(gRebinF1sAmount);
		TF1* locFunc = Create_FitFunc(locHist);
		locHist->Fit(locFunc, "QR");

		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_PSC != NULL)
	{
		TH1I* locHist = locHist_RFDeltaT_FDC_PSC;
		locHist->Rebin(gRebinF1sAmount);
		TF1* locFunc = Create_FitFunc(locHist);
		locHist->Fit(locFunc, "QR");

		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_PSC_TAGH != NULL)
	{
		TH1I* locHist = locHist_RFDeltaT_PSC_TAGH;
		locHist->Rebin(gRebinF1sAmount);
		TF1* locFunc = Create_FitFunc(locHist);
		locHist->Fit(locFunc, "QR");

		locHist->GetXaxis()->SetTitleSize(0.05);
		locHist->GetYaxis()->SetTitleSize(0.05);
		locHist->GetXaxis()->SetLabelSize(0.05);
		locHist->GetYaxis()->SetLabelSize(0.05);
		locHist->Draw();
	}

	//Fine time offsets (with respect to TOF)
	double locTimeOffset_TOF = 0.0;

	double locTimeOffset_FDC = locTOFFitFunc_FDC->GetParameter(1);
	if(!(locHist_RFDeltaT_FDC_TOF->GetEntries() > 0.0))
		locTimeOffset_FDC = 0.0;

	double locTimeOffset_PSC = locTOFFitFunc_PSC->GetParameter(1);
	if(!(locHist_RFDeltaT_PSC_TOF->GetEntries() > 0.0))
		locTimeOffset_PSC = 0.0;

	double locTimeOffset_TAGH = locTOFFitFunc_TAGH->GetParameter(1);
	if(!(locHist_RFDeltaT_TAGH_TOF->GetEntries() > 0.0))
		locTimeOffset_TAGH = 0.0;

	//Fine time offset variances
	double locTimeOffsetVariance_TOF = 0.0;
	double locTimeOffsetVariance_FDC = locTOFFitFunc_FDC->GetParError(1)*locTOFFitFunc_FDC->GetParError(1);
	double locTimeOffsetVariance_PSC = locTOFFitFunc_PSC->GetParError(1)*locTOFFitFunc_PSC->GetParError(1);
	double locTimeOffsetVariance_TAGH = locTOFFitFunc_TAGH->GetParError(1)*locTOFFitFunc_TAGH->GetParError(1);

	//Print Fine offsets to screen:
	cout << "Fine offsets for TOF/TAGH/PSC/FDC are: 0.0, " << locTimeOffset_TAGH << ", " << locTimeOffset_PSC << ", " << locTimeOffset_FDC << endl;

	//Pipe the current constants into this macro
		//NOTE: This dumps the "LATEST" values. If you need something else, modify this script.
	ostringstream locCommandStream;
	locCommandStream << "ccdb dump PHOTON_BEAM/RF/time_offset -r " << locRunNumber;
	FILE* locInputFile = gSystem->OpenPipe(locCommandStream.str().c_str(), "r");
	if(locInputFile == NULL)
		return 0;

	//get the first (comment) line
	char buff[1024]; // I HATE char buffers
	if(fgets(buff, sizeof(buff), locInputFile) == NULL)
		return 0;

	//get the second (data) line
	if(fgets(buff, sizeof(buff), locInputFile) == NULL)
		return 0;
	istringstream locConstantsStream(buff);

	//Close the pipe
	gSystem->ClosePipe(locInputFile);

	//extract the numbers
	double locCoarseTimeOffset_TOF = 0.0, locCoarseTimeOffset_TAGH = 0.0, locCoarseTimeOffset_PSC = 0.0, locCoarseTimeOffset_FDC = 0.0;
	if(!(locConstantsStream >> locCoarseTimeOffset_TOF))
		return 0;
	if(!(locConstantsStream >> locCoarseTimeOffset_TAGH))
		return 0;
	if(!(locConstantsStream >> locCoarseTimeOffset_PSC))
		return 0;
	if(!(locConstantsStream >> locCoarseTimeOffset_FDC))
		return 0;

	//Print old Coarse offsets to screen:
	cout << "Old coarse offsets for TOF/TAGH/PSC/FDC are: " << locCoarseTimeOffset_TOF << ", " << locCoarseTimeOffset_TAGH << ", " << locCoarseTimeOffset_PSC << ", " << locCoarseTimeOffset_FDC << endl;

	//Compute new time offsets
	locTimeOffset_TOF += locCoarseTimeOffset_TOF;
	locTimeOffset_TAGH += locCoarseTimeOffset_TAGH;
	locTimeOffset_PSC += locCoarseTimeOffset_PSC;
	locTimeOffset_FDC += locCoarseTimeOffset_FDC;

	//Print Final offsets to screen:
	cout << "Final offsets for TOF/TAGH/PSC/FDC are: 0.0, " << locTimeOffset_TAGH << ", " << locTimeOffset_PSC << ", " << locTimeOffset_FDC << endl;

	//Print Final variances to screen:
	cout << "Final offset variances for TOF/TAGH/PSC/FDC are: 0.0, " << locTimeOffsetVariance_TAGH << ", " << locTimeOffsetVariance_PSC << ", " << locTimeOffsetVariance_FDC << endl;

	//Print Fine offsets to file:
	ofstream locOutputFileStream;
	locOutputFileStream.open("rf_fine_time_offsets.txt");
	locOutputFileStream << "0.0 " << std::setprecision(8) << locTimeOffset_TAGH << " " << locTimeOffset_PSC << " " << locTimeOffset_FDC << endl;
	locOutputFileStream.close();

	//Print Fine offset variances to file:
	ofstream locOutputFileStream2;
	locOutputFileStream2.open("rf_time_offset_vars.txt");
	locOutputFileStream2 << "0.0 " << std::setprecision(8) << locTimeOffsetVariance_TAGH << " " << locTimeOffsetVariance_PSC << " " << locTimeOffsetVariance_FDC << endl;
	locOutputFileStream2.close();
}
