// hnamepath: /RF/AbsoluteDeltaT_RF_OtherRFs/RFDeltaT_FDC_TOF
// hnamepath: /RF/AbsoluteDeltaT_RF_OtherRFs/RFDeltaT_FDC_TAGH
// hnamepath: /RF/AbsoluteDeltaT_RF_OtherRFs/RFDeltaT_FDC_PSC
// hnamepath: /RF/AbsoluteDeltaT_RF_OtherRFs/RFDeltaT_PSC_TAGH
// hnamepath: /RF/AbsoluteDeltaT_RF_OtherRFs/RFDeltaT_PSC_TOF
// hnamepath: /RF/AbsoluteDeltaT_RF_OtherRFs/RFDeltaT_TAGH_TOF

double Calc_Mean(TH1I* locHist)
{
        double locTotal = 0.0;
        Int_t locEntries = 0;
        for(Int_t loc_i = 1; loc_i <= locHist->GetNbinsX(); ++loc_i)
        {
                Int_t locThisBinEntries = locHist->GetBinContent(loc_i);
                locEntries += locThisBinEntries;
                locTotal += double(locThisBinEntries)*locHist->GetBinCenter(loc_i);
        }
        if(locEntries == 0)
        	return 0.0;
        return (locTotal / double(locEntries));
}

int RFMacro_CoarseTimeOffsets(int locRunNumber, string locVariation = "default")
{
	//INPUT RUN NUMBER MUST BE A RUN NUMBER IN THE RUN RANGE YOU ARE TRYING TO COMMIT CONSTANTS TO
	gDirectory->cd("/"); //return to file base directory

	//Goto Beam Path
	TDirectory *locDirectory = (TDirectory*)gDirectory->FindObjectAny("RF");
	if(!locDirectory)
		return 0;
	locDirectory->cd();

	//Get RF DeltaT Histograms
	gDirectory->cd("AbsoluteDeltaT_RF_OtherRFs");
	TH1I* locHist_RFDeltaT_FDC_TOF = (TH1I*)gDirectory->Get("RFDeltaT_FDC_TOF");
	TH1I* locHist_RFDeltaT_FDC_TAGH = (TH1I*)gDirectory->Get("RFDeltaT_FDC_TAGH");
	TH1I* locHist_RFDeltaT_FDC_PSC = (TH1I*)gDirectory->Get("RFDeltaT_FDC_PSC");
	TH1I* locHist_RFDeltaT_PSC_TAGH = (TH1I*)gDirectory->Get("RFDeltaT_PSC_TAGH");
	TH1I* locHist_RFDeltaT_PSC_TOF = (TH1I*)gDirectory->Get("RFDeltaT_PSC_TOF");
	TH1I* locHist_RFDeltaT_TAGH_TOF = (TH1I*)gDirectory->Get("RFDeltaT_TAGH_TOF");

	//New coarse time offsets (with respect to TOF)
	double locTimeOffset_TOF = 0.0;
	double locTimeOffset_FDC = Calc_Mean(locHist_RFDeltaT_FDC_TOF);
	double locTimeOffset_PSC = Calc_Mean(locHist_RFDeltaT_PSC_TOF);
	double locTimeOffset_TAGH = Calc_Mean(locHist_RFDeltaT_TAGH_TOF);

	//Print extracted coarse offsets to screen:
	cout << "Extracted offsets for TOF/TAGH/PSC/FDC are: 0.0, " << locTimeOffset_TAGH << ", " << locTimeOffset_PSC << ", " << locTimeOffset_FDC << endl;

	//Pipe the current constants into this macro
		//NOTE: This dumps the "LATEST" values. If you need something else, modify this script.
	ostringstream locCommandStream;
	locCommandStream << "ccdb -v " << locVariation << " dump PHOTON_BEAM/RF/time_offset -r " << locRunNumber;
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
	double locOldTimeOffset_TOF = 0.0, locOldTimeOffset_TAGH = 0.0, locOldTimeOffset_PSC = 0.0, locOldTimeOffset_FDC = 0.0;
	if(!(locConstantsStream >> locOldTimeOffset_TOF))
		return 0;
	if(!(locConstantsStream >> locOldTimeOffset_TAGH))
		return 0;
	if(!(locConstantsStream >> locOldTimeOffset_PSC))
		return 0;
	if(!(locConstantsStream >> locOldTimeOffset_FDC))
		return 0;

	//Print old time offsets to screen:
	cout << "Old time offsets for TOF/TAGH/PSC/FDC are: " << locOldTimeOffset_TOF << ", " << locOldTimeOffset_TAGH << ", " << locOldTimeOffset_PSC << ", " << locOldTimeOffset_FDC << endl;

	//Compute new time offsets
	locTimeOffset_TOF += locOldTimeOffset_TOF;
	locTimeOffset_TAGH += locOldTimeOffset_TAGH;
	locTimeOffset_PSC += locOldTimeOffset_PSC;
	locTimeOffset_FDC += locOldTimeOffset_FDC;

	//Print new coarse time offsets to screen:
	cout << "New coarse time offsets for TOF/TAGH/PSC/FDC are: " << locTimeOffset_TOF << ", " << locTimeOffset_TAGH << ", " << locTimeOffset_PSC << ", " << locTimeOffset_FDC << endl;

	//Print new coarse offsets to file:
	ofstream locOutputFileStream;
	locOutputFileStream.open("rf_coarse_time_offsets.txt");
	locOutputFileStream << "0.0 " << std::setprecision(8) << locTimeOffset_TAGH << " " << locTimeOffset_PSC << " " << locTimeOffset_FDC << endl;
	locOutputFileStream.close();

	//Get/Make Canvas
	TCanvas *locCanvas = NULL;
	if(TVirtualPad::Pad() == NULL)
		locCanvas = new TCanvas("RF_CoarseTimeOffsets", "RF_CoarseTimeOffsets", 1200, 800); //for testing
	else
		locCanvas = gPad->GetCanvas();
	locCanvas->Divide(3, 2);

	//Draw
	locCanvas->cd(1);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_TAGH_TOF != NULL)
	{
		locHist_RFDeltaT_TAGH_TOF->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_TAGH_TOF->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_TAGH_TOF->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_TAGH_TOF->Draw();
	}

	locCanvas->cd(2);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_PSC_TOF != NULL)
	{
		locHist_RFDeltaT_PSC_TOF->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TOF->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TOF->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_PSC_TOF->Draw();
	}

	locCanvas->cd(3);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_TOF != NULL)
	{
		locHist_RFDeltaT_FDC_TOF->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TOF->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TOF->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_TOF->Draw();
	}

	locCanvas->cd(4);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_TAGH != NULL)
	{
		locHist_RFDeltaT_FDC_TAGH->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TAGH->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_TAGH->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_TAGH->Draw();
	}

	locCanvas->cd(5);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_FDC_PSC != NULL)
	{
		locHist_RFDeltaT_FDC_PSC->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_PSC->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_FDC_PSC->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_FDC_PSC->Draw();
	}

	locCanvas->cd(6);
	gPad->SetTicks();
	gPad->SetGrid();
	if(locHist_RFDeltaT_PSC_TAGH != NULL)
	{
		locHist_RFDeltaT_PSC_TAGH->GetXaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TAGH->GetYaxis()->SetTitleSize(0.05);
		locHist_RFDeltaT_PSC_TAGH->GetYaxis()->SetLabelSize(0.05);
		locHist_RFDeltaT_PSC_TAGH->Draw();
	}
}
