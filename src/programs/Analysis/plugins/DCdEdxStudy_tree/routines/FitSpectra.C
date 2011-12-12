int gRoutineSetupFlag = 0;
int gDrawFlag = 2; //0 = none, 1 = all, 2 = all but 1d-fits, 3 = only final fit

void Setup_Routines(){
	if(gRoutineSetupFlag == 1)
		return;

	gSystem->Load("$(DECAYDRAW_PATH)/depends/$(OS_NAME)/libDecayDraw.so");
	gROOT->LoadMacro("$(DECAYDRAW_PATH)/DecayDrawHistograms.C");
	gROOT->LoadMacro("$(DECAYDRAW_PATH)/DecayDrawUtilities.C");
	Setup_dqfFitter();

	DecayDrawInstructions *locDrawInstructions;
	//default 1d mass fit & counts (widescreen)
	locDrawInstructions = new DecayDrawInstructions();
	locDrawInstructions->Set_OptFit(102);
	locDrawInstructions->Set_NumPadsPerCanvas(25);
	locDrawInstructions->Set_CanvasWidth(1400);
	Set_ddhTemplateDrawInstructions(locDrawInstructions, 0);
	//fit result graphs
	locDrawInstructions = new DecayDrawInstructions();
	locDrawInstructions->Set_DrawOption("ap");
	locDrawInstructions->Set_OptFit(102);
	locDrawInstructions->Set_HistOrGraphColor(1);
	locDrawInstructions->Set_MarkerStyle(20);
	locDrawInstructions->Set_MarkerSize(1);
	Set_ddhTemplateDrawInstructions(locDrawInstructions, 1);

	gRoutineSetupFlag = 1;
}

void Draw_ConfidenceLevels(const string& locFileName, int locCDCFlag){
	Setup_Routines();

	TH1F* locHist;
	TFile *locFile = new TFile(locFileName.c_str(), "READ");
	unsigned int locNumHitsBins = (locCDCFlag == 1) ? 6 : 4;
	unsigned int locNumBetaGammaBins = 6;
	ostringstream locHistName;

	DecayDrawInstructions *locDrawInstructions;
	locDrawInstructions = new DecayDrawInstructions();
	locDrawInstructions->Set_OptFit(102);
	locDrawInstructions->Set_NumPadsPerCanvas(25);
	locDrawInstructions->Set_CanvasWidth(1400);
	locDrawInstructions->Set_HistOrGraphColor(kTeal);
	locDrawInstructions->Set_LogYFlag(1);
/*
	if(locCDCFlag == 1)
		locHist = (TH1F*)locFile->Get("dSelectorHist_ConfidenceLevel_CDC");
	else
		locHist = (TH1F*)locFile->Get("dSelectorHist_ConfidenceLevel_FDC");
	DrawAndSave_ddhHistOrGraph<TH1F>(locHist, locDrawInstructions, 0);

	for(unsigned int loc_i = 0; loc_i < locNumHitsBins; loc_i++){
		locHistName.str("");
		locHistName << "dSelectorHist_ConfidenceLevel_";
		if(locCDCFlag == 1)
			locHistName << "CDC_" << 2*loc_i + 4 << "Hits";
		else
			locHistName << "FDC_" << 3*loc_i + 3 << "Hits";
		locHist = (TH1F*)locFile->Get(locHistName.str().c_str());
		DrawAndSave_ddhHistOrGraph<TH1F>(locHist, locDrawInstructions, 0);
	}
*/
	for(unsigned int loc_i = 0; loc_i < locNumBetaGammaBins; loc_i++){
		locHistName.str("");
		locHistName << "dSelectorHist_ConfidenceLevel_";
		if(locCDCFlag == 1)
			locHistName << "CDC_BetaGammaBin" << loc_i + 1;
		else
			locHistName << "FDC_BetaGammaBin" << loc_i + 1;
		locHist = (TH1F*)locFile->Get(locHistName.str().c_str());
		DrawAndSave_ddhHistOrGraph<TH1F>(locHist, locDrawInstructions, 0);
	}
}

void Fit_dEdxUncertainty_Mean(const string& locFileName, int locParticleFlag, int locCDCFlag, int locDrawFlag = 1){
	Setup_Routines();
	gDrawFlag = locDrawFlag;

	TFile *locFile = new TFile(locFileName.c_str(), "READ");
	TH2F *locHist;
	TF1 *locFunc;

	float locMass, locMaximumMomentum, locMinBetaGamma, locMaxBeta, locMaxBetaGamma;
	unsigned int locNumBinsToMergeX = 1, locNumBinsToMergeY = 1;
	int locFitFlag;
	TGraphErrors *locGraphErrors_Mean, *locGraphErrors_Sigma;
	TObjArray *locdEdxArray_2D, *locdEdxArray_Means, *locdEdxArray_Sigmas;
	DecayDrawInstructions* locDrawInstructions;
	ostringstream locHistName, locFuncName;
	double locFitRangeMin, locFitRangeMax;
	vector<double> locInputParams;
	double locdEdxAxisMin_Mean, locdEdxAxisMax_Mean, locdEdxAxisMin_Sigma, locdEdxAxisMax_Sigma, locdEdxAxisMax_2D;

	if(locParticleFlag == 1){ //proton
		locMass = 0.93827203;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.3;
		locNumBinsToMergeX = 4;
		locNumBinsToMergeY = 2;
		locFitFlag = 14;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.0;
		locdEdxAxisMax_Mean = 15.0;
		locdEdxAxisMin_Sigma = 0.0;
		locdEdxAxisMax_Sigma = 4.0;
		locFitRangeMin = 0.3;
		locFitRangeMax = 2.1;

		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.2;  locInputParams[3] = 0.0;
	}
	if(locParticleFlag == 2){ //pi- or pi+
		locMass = 0.1395700;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.0;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.3;
		locdEdxAxisMax_Mean = 1.8;
		locdEdxAxisMin_Sigma = 0.1;
		locdEdxAxisMax_Sigma = 0.65;
		locNumBinsToMergeX = 32;
		locNumBinsToMergeY = 1;
		locFitRangeMin = 1.0;
		locFitRangeMax = 14.5;
		locFitFlag = 9;

		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.3;  locInputParams[3] = 0.003;
	}
	if(locParticleFlag == 3){ //k+
		locMass = 0.493677;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.0;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.2;
		locdEdxAxisMax_Mean = 2.6;
		locdEdxAxisMin_Sigma = 0.1;
		locdEdxAxisMax_Sigma = 0.65;
		locNumBinsToMergeX = 8;
		locNumBinsToMergeY = 1;
		locFitRangeMin = (locCDCFlag == 0) ? 0.9 : 0.65;
		locFitRangeMax = 4.0;
		locFitFlag = 9;

		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.3;  locInputParams[3] = 0.003;
	}
	locMaxBeta = locMaximumMomentum/sqrt(locMaximumMomentum*locMaximumMomentum + locMass*locMass);
	locMaxBetaGamma = locMaxBeta/sqrt(1.0 - locMaxBeta*locMaxBeta);

	if(locCDCFlag == 1)
		locHist = (TH2F*)locFile->Get("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC");
	else
		locHist = (TH2F*)locFile->Get("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC");
	locHist->Rebin2D(locNumBinsToMergeX, locNumBinsToMergeY);
	locHist->GetXaxis()->SetRangeUser(locMinBetaGamma, locMaxBetaGamma);
	locHist->GetYaxis()->SetRangeUser(0.0, locdEdxAxisMax_2D);
	Fit_ToGaussians(locHist, locGraphErrors_Mean, locGraphErrors_Sigma, locFitFlag);
	locGraphErrors_Mean->SetTitle("");
	locGraphErrors_Sigma->SetTitle("");


	if((locParticleFlag == 1) || (locParticleFlag == 3)){ //p, k+
		Set_dqfDrawFuncsFlag(3);
		Fit_dqfSignalPeak_Graph(locGraphErrors_Mean, 5, 1, locFitRangeMin, locFitRangeMax, locInputParams);
		locFuncName.str("");
		locFuncName << locGraphErrors_Mean->GetName() << "_fitfunc";
		locFunc = locGraphErrors_Mean->GetFunction(locFuncName.str().c_str());
		cout << locFunc->GetParameter(1) << ", " << locFunc->GetParameter(2) << ", " << locFunc->GetParameter(4) << ", " << locFunc->GetParameter(5) << endl;
		Set_dqfDrawFuncsFlag(0);
	}

/*
   locFunc = new TF1("gead", "[0]*log(x*x)/(x*x) + [1] + [2]*x", 3.0, 14.0);
	locGraphErrors_Mean->Fit(locFunc, "QRN+SEX0");
	locGraphErrors_Mean->GetListOfFunctions()->Add(locFunc);
*/

	locDrawInstructions = new DecayDrawInstructions();
	*locDrawInstructions = *(Get_ddhTemplateDrawInstructions(1));

	if((gDrawFlag == 1) || (gDrawFlag == 2))
		DrawAndSave_ddhHistOrGraph<TH2F>(locHist, 0, "COLZ");
	if((gDrawFlag == 1) || (gDrawFlag == 2) || (gDrawFlag == 3)){
		locDrawInstructions->Set_YRangeMin(locdEdxAxisMin_Mean);
		locDrawInstructions->Set_YRangeMax(locdEdxAxisMax_Mean);
		DrawAndSave_ddhHistOrGraph<TGraphErrors>(locGraphErrors_Mean, Get_ddhTemplateDrawInstructions(1), 1);
		locDrawInstructions->Set_YRangeMin(locdEdxAxisMin_Sigma);
		locDrawInstructions->Set_YRangeMax(locdEdxAxisMax_Sigma);
		DrawAndSave_ddhHistOrGraph<TGraphErrors>(locGraphErrors_Sigma, Get_ddhTemplateDrawInstructions(1), 1);
	}

}

void Fit_dEdxUncertainty_FDC_NumHitsBinCompare(const string& locFileName, int locParticleFlag, int locDrawFlag = 1);
void Fit_dEdxUncertainty_FDC_NumHitsBinCompare(const string& locFileName, int locParticleFlag, int locDrawFlag){
	Setup_Routines();
	gDrawFlag = locDrawFlag;

	unsigned int locNumNumHitsBins = 4;
	TFile *locFile = new TFile(locFileName.c_str(), "READ");
	TH2F *locHist;
	TF1 *locFunc;

	float locMass, locMaximumMomentum, locMinBetaGamma, locMaxBeta, locMaxBetaGamma;
	unsigned int locNumBinsToMergeX = 1, locNumBinsToMergeY = 1;
	int locFitFlag;
	TGraphErrors *locGraphErrors_Mean, *locGraphErrors_Sigma;
	TObjArray *locdEdxArray_2D, *locdEdxArray_Means, *locdEdxArray_Sigmas;
	vector<int> locHistColors(locNumNumHitsBins, 1);
	vector<int> locFitFlags(locNumNumHitsBins, 0);
	vector<string> locLegendNames(locNumNumHitsBins, "");
	DecayDrawInstructions* locDrawInstructions;
	ostringstream locHistName, locFuncName;
	double locFitRangeMin, locFitRangeMax;
	vector<double> locInputParams;
	double locdEdxAxisMin_Mean, locdEdxAxisMax_Mean, locdEdxAxisMin_Sigma, locdEdxAxisMax_Sigma, locdEdxAxisMax_2D;
	vector<int> locNumBinsToMergeXVector(locNumNumHitsBins, 1);
	vector<double> locFitRangeMinVector(locNumNumHitsBins, 0.0);

	if(locParticleFlag == 1){ //proton
		locMass = 0.93827203;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.3;
		locNumBinsToMergeX = 4;
		locNumBinsToMergeY = 2;
		locFitFlag = 1;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.0;
		locdEdxAxisMax_Mean = 15.0;
		locdEdxAxisMin_Sigma = 0.0;
		locdEdxAxisMax_Sigma = 4.0;
//		locFitRangeMin = 0.6;
		locFitRangeMax = 2.1;
		locFitRangeMinVector[0] = 0.5;  locFitRangeMinVector[1] = 0.5;  locFitRangeMinVector[2] = 0.35;  locFitRangeMinVector[3] = 0.48;

		locFitFlags[0] = 1;  locFitFlags[1] = 2;  locFitFlags[2] = 4;  locFitFlags[3] = 2;
		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.2;  locInputParams[3] = 0.0;
		locNumBinsToMergeXVector[0] = 8;  locNumBinsToMergeXVector[1] = 8;  locNumBinsToMergeXVector[2] = 8;  locNumBinsToMergeXVector[3] = 4;
	}
	if(locParticleFlag == 2){ //pi- or pi+
		locMass = 0.1395700;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.0;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.3;
		locdEdxAxisMax_Mean = 1.8;
		locdEdxAxisMin_Sigma = 0.1;
		locdEdxAxisMax_Sigma = 0.65;
		locNumBinsToMergeX = 32;
		locNumBinsToMergeY = 1;
//		locFitRangeMin = 1.0;
		locFitRangeMax = 14.5;
		locFitRangeMinVector[0] = 1.0;  locFitRangeMinVector[1] = 0.4;  locFitRangeMinVector[2] = 0.4;  locFitRangeMinVector[3] = 0.75;

		locFitFlags[0] = 8;  locFitFlags[1] = 3;  locFitFlags[2] = 9;  locFitFlags[3] = 9;
		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.3;  locInputParams[3] = 0.003;
		locNumBinsToMergeXVector[0] = 32;  locNumBinsToMergeXVector[1] = 32;  locNumBinsToMergeXVector[2] = 32;  locNumBinsToMergeXVector[3] = 32;
	}
	if(locParticleFlag == 3){ //K+
		locMass = 0.493677;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.0;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.3;
		locdEdxAxisMax_Mean = 1.8;
		locdEdxAxisMin_Sigma = 0.1;
		locdEdxAxisMax_Sigma = 0.65;
		locNumBinsToMergeX = 8;
		locNumBinsToMergeY = 1;
//		locFitRangeMin = 1.0;
		locFitRangeMax = 4.1;
		locFitRangeMinVector[0] = 1.0;  locFitRangeMinVector[1] = 1.0;  locFitRangeMinVector[2] = 1.0;  locFitRangeMinVector[3] = 1.0;

		locFitFlags[0] = 15;  locFitFlags[1] = 16;  locFitFlags[2] = 17;  locFitFlags[3] = 18;
		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.3;  locInputParams[3] = 0.003;
		locNumBinsToMergeXVector[0] = 16;  locNumBinsToMergeXVector[1] = 16;  locNumBinsToMergeXVector[2] = 16;  locNumBinsToMergeXVector[3] = 16;
	}
	locMaxBeta = locMaximumMomentum/sqrt(locMaximumMomentum*locMaximumMomentum + locMass*locMass);
	locMaxBetaGamma = locMaxBeta/sqrt(1.0 - locMaxBeta*locMaxBeta);

	TObjArray *locdEdxArray_2D = new TObjArray(locNumNumHitsBins);
	TObjArray *locdEdxArray_Means = new TObjArray(locNumNumHitsBins);
	TObjArray *locdEdxArray_Sigmas = new TObjArray(locNumNumHitsBins);

	for(unsigned int loc_i = 0; loc_i < locNumNumHitsBins; loc_i++){
		locHistName.str("");
		if(loc_i == 0)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_3Hits";
		if(loc_i == 1)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_6Hits";
		if(loc_i == 2)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_9Hits";
		if(loc_i == 3)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_12Hits";
		locHist = (TH2F*)locFile->Get(locHistName.str().c_str());
		locHist->Rebin2D(locNumBinsToMergeXVector[loc_i], locNumBinsToMergeY);
		locHist->GetXaxis()->SetRangeUser(locMinBetaGamma, locMaxBetaGamma);
		locHist->GetYaxis()->SetRangeUser(0.0, locdEdxAxisMax_2D);
		Fit_ToGaussians(locHist, locGraphErrors_Mean, locGraphErrors_Sigma, locFitFlags[loc_i]);
		locGraphErrors_Mean->SetTitle("");
		locGraphErrors_Sigma->SetTitle("");

		Set_dqfDrawFuncsFlag(3);
		if((locParticleFlag == 1) && ((loc_i == 1) || (loc_i == 0))){ //cheat!
			double *locErrors = locGraphErrors_Sigma->GetEY();
			locErrors[3] /= 40.0;
		}

		Fit_dqfSignalPeak_Graph(locGraphErrors_Sigma, 5, 1, locFitRangeMinVector[loc_i], locFitRangeMax, locInputParams);
//		Fit_dqfSignalPeak_Graph(locGraphErrors_Sigma, 5, -1, locFitRangeMin, locFitRangeMax, locInputParams); //for exponential
//		Fit_dqfSignalPeak_Graph(locGraphErrors_Sigma, 0, 1, locFitRangeMin, locFitRangeMax, locInputParams); //for linear
		Set_dqfDrawFuncsFlag(0);
		locFuncName.str("");
		locFuncName << locGraphErrors_Sigma->GetName() << "_fitfunc";
		locFunc = locGraphErrors_Sigma->GetFunction(locFuncName.str().c_str());
		cout << locFunc->GetParameter(1) << ", " << locFunc->GetParameter(2) << ", " << locFunc->GetParameter(4) << ", " << locFunc->GetParameter(5) << endl;

		locdEdxArray_2D->AddAt(locHist, loc_i);
		locdEdxArray_Means->AddAt(locGraphErrors_Mean, loc_i);
		locdEdxArray_Sigmas->AddAt(locGraphErrors_Sigma, loc_i);
		locHistColors[loc_i] = loc_i + 1;
		locLegendNames[loc_i] = ((TH2F*)(*locdEdxArray_2D)[loc_i])->GetTitle();
	}

	locDrawInstructions = new DecayDrawInstructions();
	*locDrawInstructions = *(Get_ddhTemplateDrawInstructions(1));
	locDrawInstructions->Set_HistOrGraphColorVector(locHistColors);
	locDrawInstructions->Set_LegendNameVector(locLegendNames);
	locDrawInstructions->Set_DrawOption("p");

//	if((gDrawFlag == 1) || (gDrawFlag == 2))
//		DrawAndSave_ddhHistOrGraphArray<TH2F>(locdEdxArray_2D, 0, "COLZ");
	if((gDrawFlag == 1) || (gDrawFlag == 2) || (gDrawFlag == 3)){
		locDrawInstructions->Set_YRangeMin(locdEdxAxisMin_Mean);
		locDrawInstructions->Set_YRangeMax(locdEdxAxisMax_Mean);
//		DrawAndSave_ddhHistOrGraphSuperimposition<TGraphErrors>(locdEdxArray_Means, locDrawInstructions, 1);
		locDrawInstructions->Set_YRangeMin(locdEdxAxisMin_Sigma);
		locDrawInstructions->Set_YRangeMax(locdEdxAxisMax_Sigma);
		DrawAndSave_ddhHistOrGraphSuperimposition<TGraphErrors>(locdEdxArray_Sigmas, locDrawInstructions, 1);
	}
}

void Fit_dEdxUncertainty_CDC_NumHitsBinCompare(const string& locFileName, int locParticleFlag, int locDrawFlag = 1);
void Fit_dEdxUncertainty_CDC_NumHitsBinCompare(const string& locFileName, int locParticleFlag, int locDrawFlag){
	Setup_Routines();
	gDrawFlag = locDrawFlag;

	unsigned int locNumNumHitsBins = 6;
	TFile *locFile = new TFile(locFileName.c_str(), "READ");
	TH2F *locHist;
	TF1 *locFunc;

	float locMass, locMaximumMomentum, locMinBetaGamma, locMaxBeta, locMaxBetaGamma;
	unsigned int locNumBinsToMergeX = 1, locNumBinsToMergeY = 1;
	int locFitFlag;
	TGraphErrors *locGraphErrors_Mean, *locGraphErrors_Sigma;
	TObjArray *locdEdxArray_2D, *locdEdxArray_Means, *locdEdxArray_Sigmas;
	vector<int> locHistColors(locNumNumHitsBins, 1);
	vector<int> locFitFlags(locNumNumHitsBins, 0);
	vector<string> locLegendNames(locNumNumHitsBins, "");
	DecayDrawInstructions* locDrawInstructions;
	ostringstream locHistName, locFuncName;
	double locFitRangeMin, locFitRangeMax;
	vector<double> locInputParams;
	double locdEdxAxisMin_Mean, locdEdxAxisMax_Mean, locdEdxAxisMin_Sigma, locdEdxAxisMax_Sigma, locdEdxAxisMax_2D;
	vector<int> locNumBinsToMergeXVector(locNumNumHitsBins, 1);
	vector<double> locFitRangeMinVector(locNumNumHitsBins, 0.0);

	if(locParticleFlag == 1){ //proton
		locMass = 0.93827203;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.3;
		locFitFlag = 1;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.0;
		locdEdxAxisMax_Mean = 15.0;
		locdEdxAxisMin_Sigma = 0.0;
		locdEdxAxisMax_Sigma = 3.5;
		locFitRangeMin = 0.45;
		locFitRangeMax = 2.15;
		locFitRangeMinVector[0] = 0.42;  locFitRangeMinVector[1] = 0.42;  locFitRangeMinVector[2] = 0.25;
		locFitRangeMinVector[3] = 0.25;  locFitRangeMinVector[4] = 0.25;  locFitRangeMinVector[5] = 0.25;

		locFitFlags[0] = 12;  locFitFlags[1] = 13;  locFitFlags[2] = 7;  locFitFlags[3] = 10;  locFitFlags[4] = 10;  locFitFlags[5] = 10;
		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.2;  locInputParams[3] = 0.0;
		locNumBinsToMergeXVector[0] = 8;  locNumBinsToMergeXVector[1] = 8;  locNumBinsToMergeXVector[2] = 4;
		locNumBinsToMergeXVector[3] = 4;  locNumBinsToMergeXVector[4] = 4;  locNumBinsToMergeXVector[5] = 4;
	}
	if(locParticleFlag == 2){ //pi- or pi+
		locMass = 0.1395700;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.0;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.3;
		locdEdxAxisMax_Mean = 1.8;
		locdEdxAxisMin_Sigma = 0.1;
		locdEdxAxisMax_Sigma = 0.65;
		locFitRangeMin = 1.0;
		locFitRangeMax = 14.5;

		locFitRangeMinVector[0] = 0.8;  locFitRangeMinVector[1] = 0.8;  locFitRangeMinVector[2] = 0.8;
		locFitRangeMinVector[3] = 0.8;  locFitRangeMinVector[4] = 0.8;  locFitRangeMinVector[5] = 0.8;
		locFitFlags[0] = 5;  locFitFlags[1] = 6;  locFitFlags[2] = 6;  locFitFlags[3] = 11;  locFitFlags[4] = 11;  locFitFlags[5] = 11;
		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.3;  locInputParams[3] = 0.003;
		locNumBinsToMergeXVector[0] = 64;  locNumBinsToMergeXVector[1] = 64;  locNumBinsToMergeXVector[2] = 64;
		locNumBinsToMergeXVector[3] = 64;  locNumBinsToMergeXVector[4] = 32;  locNumBinsToMergeXVector[5] = 32;
	}
	if(locParticleFlag == 3){ //k+
		locMass = 0.493677;
		locMaximumMomentum = 2.0;
		locMinBetaGamma = 0.0;
		locdEdxAxisMax_2D = 20.0;
		locdEdxAxisMin_Mean = 1.3;
		locdEdxAxisMax_Mean = 1.8;
		locdEdxAxisMin_Sigma = 0.1;
		locdEdxAxisMax_Sigma = 0.9;
		locFitRangeMin = 1.0;
		locFitRangeMax = 4.1;

		locFitRangeMinVector[0] = 1.0;  locFitRangeMinVector[1] = 1.0;  locFitRangeMinVector[2] = 1.0;
		locFitRangeMinVector[3] = 1.0;  locFitRangeMinVector[4] = 1.0;  locFitRangeMinVector[5] = 1.0;
		locFitFlags[0] = 19;  locFitFlags[1] = 20;  locFitFlags[2] = 21;  locFitFlags[3] = 22;  locFitFlags[4] = 22;  locFitFlags[5] = 22;
		locInputParams.resize(4);
		locInputParams[0] = 1.4;  locInputParams[1] = -1.5;  locInputParams[2] = 0.3;  locInputParams[3] = 0.003;
		locNumBinsToMergeXVector[0] = 16;  locNumBinsToMergeXVector[1] = 16;  locNumBinsToMergeXVector[2] = 16;
		locNumBinsToMergeXVector[3] = 16;  locNumBinsToMergeXVector[4] = 16;  locNumBinsToMergeXVector[5] = 16;
	}
	locMaxBeta = locMaximumMomentum/sqrt(locMaximumMomentum*locMaximumMomentum + locMass*locMass);
	locMaxBetaGamma = locMaxBeta/sqrt(1.0 - locMaxBeta*locMaxBeta);

	TObjArray *locdEdxArray_2D = new TObjArray(locNumNumHitsBins);
	TObjArray *locdEdxArray_Means = new TObjArray(locNumNumHitsBins);
	TObjArray *locdEdxArray_Sigmas = new TObjArray(locNumNumHitsBins);

	for(unsigned int loc_i = 0; loc_i < locNumNumHitsBins; loc_i++){
		locHistName.str("");

		if(loc_i == 0)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_4Hits";
		if(loc_i == 1)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_6Hits";
		if(loc_i == 2)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_8Hits";
		if(loc_i == 3)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_10Hits";
		if(loc_i == 4)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_12Hits";
		if(loc_i == 5)
			locHistName << "dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_14Hits";

		locHist = (TH2F*)locFile->Get(locHistName.str().c_str());
		locHist->Rebin2D(locNumBinsToMergeXVector[loc_i], locNumBinsToMergeY);
		locHist->GetXaxis()->SetRangeUser(locMinBetaGamma, locMaxBetaGamma);
		locHist->GetYaxis()->SetRangeUser(0.0, locdEdxAxisMax_2D);
		Fit_ToGaussians(locHist, locGraphErrors_Mean, locGraphErrors_Sigma, locFitFlags[loc_i]);
		locGraphErrors_Mean->SetTitle("");
		locGraphErrors_Sigma->SetTitle("");

		Set_dqfDrawFuncsFlag(3);
		Fit_dqfSignalPeak_Graph(locGraphErrors_Sigma, 5, 1, locFitRangeMinVector[loc_i], locFitRangeMax, locInputParams);
//		Fit_dqfSignalPeak_Graph(locGraphErrors_Sigma, 5, -1, locFitRangeMin, locFitRangeMax, locInputParams); //for exponential
//		Fit_dqfSignalPeak_Graph(locGraphErrors_Sigma, 0, 1, locFitRangeMin, locFitRangeMax, locInputParams); //for linear
		Set_dqfDrawFuncsFlag(0);
		locFuncName.str("");
		locFuncName << locGraphErrors_Sigma->GetName() << "_fitfunc";
		locFunc = locGraphErrors_Sigma->GetFunction(locFuncName.str().c_str());
		cout << locFunc->GetParameter(1) << ", " << locFunc->GetParameter(2) << ", " << locFunc->GetParameter(4) << ", " << locFunc->GetParameter(5) << endl;

		locdEdxArray_2D->AddAt(locHist, loc_i);
		locdEdxArray_Means->AddAt(locGraphErrors_Mean, loc_i);
		locdEdxArray_Sigmas->AddAt(locGraphErrors_Sigma, loc_i);
		locHistColors[loc_i] = loc_i + 1;
		locLegendNames[loc_i] = ((TH2F*)(*locdEdxArray_2D)[loc_i])->GetTitle();
	}

	locDrawInstructions = new DecayDrawInstructions();
	*locDrawInstructions = *(Get_ddhTemplateDrawInstructions(1));
	locDrawInstructions->Set_HistOrGraphColorVector(locHistColors);
	locDrawInstructions->Set_LegendNameVector(locLegendNames);
	locDrawInstructions->Set_DrawOption("p");

//	if((gDrawFlag == 1) || (gDrawFlag == 2))
//		DrawAndSave_ddhHistOrGraphArray<TH2F>(locdEdxArray_2D, 0, "COLZ");
	if((gDrawFlag == 1) || (gDrawFlag == 2) || (gDrawFlag == 3)){
		locDrawInstructions->Set_YRangeMin(locdEdxAxisMin_Mean);
		locDrawInstructions->Set_YRangeMax(locdEdxAxisMax_Mean);
//		DrawAndSave_ddhHistOrGraphSuperimposition<TGraphErrors>(locdEdxArray_Means, locDrawInstructions, 1);
		locDrawInstructions->Set_YRangeMin(locdEdxAxisMin_Sigma);
		locDrawInstructions->Set_YRangeMax(locdEdxAxisMax_Sigma);
		DrawAndSave_ddhHistOrGraphSuperimposition<TGraphErrors>(locdEdxArray_Sigmas, locDrawInstructions, 1);
	}
}

void Fit_ToGaussians(const TH2F *loc2FHist, TGraphErrors *&locGraphErrors_Mean, TGraphErrors *&locGraphErrors_Sigma, int locFitFlag){

	TH1D *locFitHist;
	TF1 *locFunction;
	TGraphErrors *locGraphErrors;
	TGraph *locGraph;
	string locOriginalTitle = loc2FHist->GetTitle();
	ostringstream locHistName, locHistTitle, locFuncName, locGraphName, locGraphTitle;
	unsigned int loc_i, locNDF;
	float locFitRangeMin, locFitRangeMax, locFitRange, locXOfPeak, locChiSq;
	float locFitRangeMin_DistanceFromPeak, locFitRangeMax_DistanceFromPeak;

	int locNumMinCountsToFit = 200;

	int locMinBin = loc2FHist->GetXaxis()->GetFirst();
	int locMaxBin = loc2FHist->GetXaxis()->GetLast();
	int locNumBins = locMaxBin - locMinBin + 1;

	double *locXArray = new double[locNumBins];
	double *locXUncertaintyArray = new double[locNumBins];
	double *locMeanArray = new double[locNumBins];
	double *locMeanUncertaintyArray = new double[locNumBins];
	double *locSigmaArray = new double[locNumBins];
	double *locSigmaUncertaintyArray = new double[locNumBins];
	double *locChiSqPerNDFArray = new double[locNumBins];

	unsigned int locNumFitHists = 0;
	TObjArray *locFitHistArray = new TObjArray(1);
	for(loc_i = locMinBin; loc_i <= locMaxBin; loc_i++){
		locHistName.str("");
		locHistName << loc2FHist->GetName() << "_Projection_" << loc_i;
		locFitHist = (TH1D*)loc2FHist->ProjectionY(locHistName.str().c_str(), loc_i, loc_i);
		locHistTitle.str("");
		if(locOriginalTitle != "")
			locHistTitle << loc2FHist->GetTitle() << ", ";
		locHistTitle << loc2FHist->GetXaxis()->GetBinLowEdge(loc_i) << " <= " << loc2FHist->GetXaxis()->GetTitle() << " < " << loc2FHist->GetXaxis()->GetBinUpEdge(loc_i);
		locFitHist->SetTitle(locHistTitle.str().c_str());

		//rebinning
		if(locFitHist->GetEntries() < locNumMinCountsToFit)
			continue;

		//fit range
		locXOfPeak = locFitHist->GetXaxis()->GetBinCenter(locFitHist->GetMaximumBin());

		if(locFitFlag == 1){
			if(locNumFitHists >= 11){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if((locNumFitHists >= 9) && (locNumFitHists < 11)){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
			if((locNumFitHists > 6) && (locNumFitHists < 9)){
				locFitRangeMin_DistanceFromPeak = 1.35;
				locFitRangeMax_DistanceFromPeak = 1.35;
			}
			if((locNumFitHists >= 2) && (locNumFitHists <= 6)){
				locFitRangeMin_DistanceFromPeak = 1.8;
				locFitRangeMax_DistanceFromPeak = 2.1;
			}
			if(locNumFitHists < 2){
				locFitRangeMin_DistanceFromPeak = 2.5;
				locFitRangeMax_DistanceFromPeak = 2.5;
			}
		}
		if(locFitFlag == 2){
			locFitRangeMin_DistanceFromPeak = 0.45;
			locFitRangeMax_DistanceFromPeak = 0.45;
			if((locNumFitHists >= 11) && (locNumFitHists <= 20)){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if((locNumFitHists >= 9) && (locNumFitHists < 11)){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
			if((locNumFitHists > 6) && (locNumFitHists < 9)){
				locFitRangeMin_DistanceFromPeak = 1.35;
				locFitRangeMax_DistanceFromPeak = 1.35;
			}
			if((locNumFitHists >= 2) && (locNumFitHists <= 6)){
				locFitRangeMin_DistanceFromPeak = 1.8;
				locFitRangeMax_DistanceFromPeak = 2.1;
			}
			if(locNumFitHists < 2){
				locFitRangeMin_DistanceFromPeak = 2.5;
				locFitRangeMax_DistanceFromPeak = 2.5;
			}
		}
		if(locFitFlag == 3){
			locFitRangeMin_DistanceFromPeak = 0.5;
			locFitRangeMax_DistanceFromPeak = 0.5;
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
		}
		if(locFitFlag == 4){
			locFitRangeMin_DistanceFromPeak = 0.45;
			locFitRangeMax_DistanceFromPeak = 0.45;
			if((locNumFitHists >= 11) && (locNumFitHists <= 20)){
				locFitRangeMin_DistanceFromPeak = 0.55;
				locFitRangeMax_DistanceFromPeak = 0.55;
			}
			if((locNumFitHists >= 9) && (locNumFitHists < 11)){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
			if((locNumFitHists > 5) && (locNumFitHists < 9)){
				locFitRangeMin_DistanceFromPeak = 0.65;
				locFitRangeMax_DistanceFromPeak = 1.0;
			}
			if(locNumFitHists == 3){
				locFitRangeMin_DistanceFromPeak = 1.8;
				locFitRangeMax_DistanceFromPeak = 2.1;
			}
			if(locNumFitHists < 2){
				locFitRangeMin_DistanceFromPeak = 2.5;
				locFitRangeMax_DistanceFromPeak = 2.9;
			}
			if(locNumFitHists == 2){
				locFitRangeMin_DistanceFromPeak = 2.8;
				locFitRangeMax_DistanceFromPeak = 2.1;
			}
			if(locNumFitHists == 4){
				locFitHist->Rebin(4);
				locFitRangeMin_DistanceFromPeak = 1.5;
				locFitRangeMax_DistanceFromPeak = 1.5;
			}
			if(locNumFitHists == 5){
				locFitHist->Rebin(4);
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 1.8;
			}
		}
		if(locFitFlag == 5){
			locFitRangeMin_DistanceFromPeak = 0.8;
			locFitRangeMax_DistanceFromPeak = 0.6;
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
		}
		if(locFitFlag == 6){
			locFitRangeMin_DistanceFromPeak = 0.6;
			locFitRangeMax_DistanceFromPeak = 0.6;
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
			if((locNumFitHists >= 4) && (locNumFitHists < 18)){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
		}
		if(locFitFlag == 7){
			locFitRangeMin_DistanceFromPeak = 0.65;
			locFitRangeMax_DistanceFromPeak = 0.6;
			if((locNumFitHists >= 14) && (locNumFitHists < 21)){
				locFitRangeMin_DistanceFromPeak = 1.0;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if((locNumFitHists >= 2) && (locNumFitHists < 14)){
				locFitRangeMin_DistanceFromPeak = 2.2;
				locFitRangeMax_DistanceFromPeak = 2.4;
			}
			if(locNumFitHists < 2){
				locFitRangeMin_DistanceFromPeak = 3.5;
				locFitRangeMax_DistanceFromPeak = 3.5;
			}
		}
		if(locFitFlag == 8){
			locFitRangeMin_DistanceFromPeak = 0.7;
			locFitRangeMax_DistanceFromPeak = 0.7;
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 1.2;
			}
		}
		if(locFitFlag == 9){
			locFitRangeMin_DistanceFromPeak = 0.4;
			locFitRangeMax_DistanceFromPeak = 0.37;
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
		}
		if(locFitFlag == 10){
			if(locNumFitHists >= 11){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if((locNumFitHists >= 9) && (locNumFitHists < 11)){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
			if((locNumFitHists > 6) && (locNumFitHists < 9)){
				locFitRangeMin_DistanceFromPeak = 1.35;
				locFitRangeMax_DistanceFromPeak = 1.35;
			}
			if((locNumFitHists >= 2) && (locNumFitHists <= 6)){
				locFitRangeMin_DistanceFromPeak = 1.8;
				locFitRangeMax_DistanceFromPeak = 2.1;
			}
			if(locNumFitHists < 2){
				locFitRangeMin_DistanceFromPeak = 2.5;
				locFitRangeMax_DistanceFromPeak = 2.5;
			}
		}
		if(locFitFlag == 11){
			locFitRangeMin_DistanceFromPeak = 0.4;
			locFitRangeMax_DistanceFromPeak = 0.37;
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
		}
		if(locFitFlag == 12){
			locFitHist->Rebin(2);
			if(locNumFitHists >= 14){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if(locNumFitHists == 13){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if((locNumFitHists >= 11) && (locNumFitHists <= 12)){
				locFitRangeMin_DistanceFromPeak = 1.1;
				locFitRangeMax_DistanceFromPeak = 1.0;
			}
			if(locNumFitHists == 10){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 1.0;
			}
			if((locNumFitHists > 6) && (locNumFitHists <= 9)){
				locFitRangeMin_DistanceFromPeak = 2.0;
				locFitRangeMax_DistanceFromPeak = 1.3;
			}
			if((locNumFitHists >= 5) && (locNumFitHists <= 6)){
				locFitRangeMin_DistanceFromPeak = 3.0;
				locFitRangeMax_DistanceFromPeak = 2.0;
			}
			if(locNumFitHists == 4){
				locFitRangeMin_DistanceFromPeak = 3.0;
				locFitRangeMax_DistanceFromPeak = 3.0;
			}
			if((locNumFitHists >= 2) && (locNumFitHists <= 3)){
				locFitRangeMin_DistanceFromPeak = 4.5;
				locFitRangeMax_DistanceFromPeak = 4.5;
			}
			if(locNumFitHists == 1){
				locFitRangeMin_DistanceFromPeak = 4.5;
				locFitRangeMax_DistanceFromPeak = 4.5;
			}
			if(locNumFitHists == 0){
				locFitRangeMin_DistanceFromPeak = 6.5;
				locFitRangeMax_DistanceFromPeak = 3.5;
			}
		}
		if(locFitFlag == 13){
			locFitRangeMin_DistanceFromPeak = 0.55;
			locFitRangeMax_DistanceFromPeak = 0.6;
			if((locNumFitHists >= 14) && (locNumFitHists < 21)){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if((locNumFitHists >= 9) && (locNumFitHists < 14)){
				locFitRangeMin_DistanceFromPeak = 1.3;
				locFitRangeMax_DistanceFromPeak = 1.0;
			}
			if((locNumFitHists >= 7) && (locNumFitHists <= 8)){
				locFitRangeMin_DistanceFromPeak = 1.8;
				locFitRangeMax_DistanceFromPeak = 1.2;
			}
			if((locNumFitHists >= 5) && (locNumFitHists <= 6)){
				locFitRangeMin_DistanceFromPeak = 2.2;
				locFitRangeMax_DistanceFromPeak = 2.1;
			}
			if(locNumFitHists == 4){
				locFitRangeMin_DistanceFromPeak = 2.2;
				locFitRangeMax_DistanceFromPeak = 2.4;
			}
			if(locNumFitHists == 3){
				locFitRangeMin_DistanceFromPeak = 2.2;
				locFitRangeMax_DistanceFromPeak = 3.4;
			}
			if(locNumFitHists == 2){
				locFitRangeMin_DistanceFromPeak = 3.4;
				locFitRangeMax_DistanceFromPeak = 3.4;
			}
			if(locNumFitHists == 1){
				locFitRangeMin_DistanceFromPeak = 3.5;
				locFitRangeMax_DistanceFromPeak = 3.5;
			}
			if(locNumFitHists == 0){
				locFitRangeMin_DistanceFromPeak = 5.5;
				locFitRangeMax_DistanceFromPeak = 3.5;
			}
		}
		if(locFitFlag == 14){
			locFitRangeMin_DistanceFromPeak = 0.45;
			locFitRangeMax_DistanceFromPeak = 0.45;
			if((locNumFitHists >= 11) && (locNumFitHists <= 20)){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if((locNumFitHists >= 9) && (locNumFitHists < 11)){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
			if((locNumFitHists > 6) && (locNumFitHists < 9)){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 1.2;
			}
			if((locNumFitHists >= 5) && (locNumFitHists <= 6)){
				locFitRangeMin_DistanceFromPeak = 1.6;
				locFitRangeMax_DistanceFromPeak = 1.9;
			}
			if((locNumFitHists >= 3) && (locNumFitHists <= 4)){
				locFitRangeMin_DistanceFromPeak = 2.1;
				locFitRangeMax_DistanceFromPeak = 2.1;
			}
			if(locNumFitHists <= 2){
				locFitRangeMin_DistanceFromPeak = 4.0;
				locFitRangeMax_DistanceFromPeak = 4.0;
			}
		}
		if(locFitFlag == 15){
			locFitRangeMin_DistanceFromPeak = 0.7;
			locFitRangeMax_DistanceFromPeak = 0.7;
			if(locNumFitHists == 8){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
			if(locNumFitHists == 7){
				locFitRangeMin_DistanceFromPeak = 1.1;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists == 6){
				locFitRangeMin_DistanceFromPeak = 1.1;
				locFitRangeMax_DistanceFromPeak = 1.1;
			}
			if(locNumFitHists == 4){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.5;
			}
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
		}
		if(locFitFlag == 16){
			locFitRangeMin_DistanceFromPeak = 0.5;
			locFitRangeMax_DistanceFromPeak = 0.5;
			if(locNumFitHists == 8){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.5;
			}
			if(locNumFitHists == 7){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists == 6){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 1.3;
			}
			if(locNumFitHists == 3){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if(locNumFitHists < 3){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
		}
		if(locFitFlag == 17){
			locFitRangeMin_DistanceFromPeak = 0.4;
			locFitRangeMax_DistanceFromPeak = 0.37;
			if(locNumFitHists == 8){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists == 7){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if(locNumFitHists == 6){
				locFitRangeMin_DistanceFromPeak = 0.5;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if(locNumFitHists == 5){
				locFitRangeMin_DistanceFromPeak = 0.5;
				locFitRangeMax_DistanceFromPeak = 1.1;
			}
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists == 2){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
		}
		if(locFitFlag == 18){
			locFitRangeMin_DistanceFromPeak = 0.4;
			locFitRangeMax_DistanceFromPeak = 0.37;
			if((locNumFitHists >= 5) && (locNumFitHists <= 6)){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists < 3){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists == 1){
				locFitRangeMin_DistanceFromPeak = 1.0;
				locFitRangeMax_DistanceFromPeak = 1.3;
			}
		}
		if(locFitFlag == 19){
			locFitRangeMin_DistanceFromPeak = 0.8;
			locFitRangeMax_DistanceFromPeak = 0.6;
			if((locNumFitHists >= 10) && (locNumFitHists <= 14)){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if((locNumFitHists >= 8) && (locNumFitHists <= 9)){
				locFitRangeMin_DistanceFromPeak = 1.0;
				locFitRangeMax_DistanceFromPeak = 1.0;
			}
			if(locNumFitHists == 7){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 1.1;
			}
			if(locNumFitHists == 6){
				locFitRangeMin_DistanceFromPeak = 1.1;
				locFitRangeMax_DistanceFromPeak = 1.3;
			}
			if(locNumFitHists == 5){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 2.0;
			}
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
		}
		if(locFitFlag == 20){
			locFitRangeMin_DistanceFromPeak = 0.6;
			locFitRangeMax_DistanceFromPeak = 0.6;
			if((locNumFitHists >= 10) && (locNumFitHists < 18)){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if((locNumFitHists >= 8) && (locNumFitHists < 9)){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.9;
			}
			if(locNumFitHists == 7){
				locFitRangeMin_DistanceFromPeak = 1.1;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
			if(locNumFitHists == 6){
				locFitRangeMin_DistanceFromPeak = 1.1;
				locFitRangeMax_DistanceFromPeak = 1.2;
			}
			if(locNumFitHists == 5){
				locFitRangeMin_DistanceFromPeak = 1.2;
				locFitRangeMax_DistanceFromPeak = 1.2;
			}
			if(locNumFitHists == 4){
				locFitRangeMin_DistanceFromPeak = 1.4;
				locFitRangeMax_DistanceFromPeak = 1.5;
			}
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
		}
		if(locFitFlag == 21){
			locFitRangeMin_DistanceFromPeak = 0.6;
			locFitRangeMax_DistanceFromPeak = 0.6;
			if((locNumFitHists >= 8) && (locNumFitHists < 18)){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists == 7){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
			if(locNumFitHists == 6){
				locFitRangeMin_DistanceFromPeak = 1.0;
				locFitRangeMax_DistanceFromPeak = 1.0;
			}
			if(locNumFitHists == 5){
				locFitRangeMin_DistanceFromPeak = 1.3;
				locFitRangeMax_DistanceFromPeak = 1.0;
			}
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.8;
				locFitRangeMax_DistanceFromPeak = 0.8;
			}
		}
		if(locFitFlag == 22){
			locFitRangeMin_DistanceFromPeak = 0.4;
			locFitRangeMax_DistanceFromPeak = 0.37;
			if(locNumFitHists == 9){
				locFitRangeMin_DistanceFromPeak = 0.4;
				locFitRangeMax_DistanceFromPeak = 0.55;
			}
			if((locNumFitHists >= 7) && (locNumFitHists <= 8)){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
			if(locNumFitHists == 6){
				locFitRangeMin_DistanceFromPeak = 0.7;
				locFitRangeMax_DistanceFromPeak = 0.7;
			}
			if(locNumFitHists == 5){
				locFitRangeMin_DistanceFromPeak = 0.9;
				locFitRangeMax_DistanceFromPeak = 0.5;
			}
			if(locNumFitHists < 4){
				locFitRangeMin_DistanceFromPeak = 0.6;
				locFitRangeMax_DistanceFromPeak = 0.6;
			}
		}

		locFitRangeMin = locXOfPeak - locFitRangeMin_DistanceFromPeak;
		locFitRangeMax = locXOfPeak + locFitRangeMax_DistanceFromPeak;
		locFitRange = locFitRangeMax - locFitRangeMin;
		Fit_dqfSignalPeak(locFitHist, 1, -1, locFitRangeMin, locFitRangeMax);
		locFitHist->GetXaxis()->SetRangeUser(locFitRangeMin - locFitRange, locFitRangeMax + locFitRange);

      locFuncName.str("");
      locFuncName << locFitHist->GetName() << "_fitfunc";
      locFunction = locFitHist->GetFunction(locFuncName.str().c_str());

		locXArray[locNumFitHists] = loc2FHist->GetXaxis()->GetBinCenter(loc_i);
		locXUncertaintyArray[locNumFitHists] = 0.0;
		locMeanArray[locNumFitHists] = locFunction->GetParameter(2);
		locMeanUncertaintyArray[locNumFitHists] = locFunction->GetParError(2);
		locSigmaArray[locNumFitHists] = locFunction->GetParameter(3);
		locSigmaUncertaintyArray[locNumFitHists] = locFunction->GetParError(3);
		locChiSq = locFunction->GetChisquare();
		locNDF = locFunction->GetNDF();
		locChiSqPerNDFArray[locNumFitHists] = locChiSq/(float(locNDF));

		locFitHistArray->AddAtAndExpand(locFitHist, locNumFitHists);
		++locNumFitHists;
	}

	if(gDrawFlag == 1)
		DrawAndSave_ddhHistOrGraphArray<TH1D>(locFitHistArray, Get_ddhTemplateDrawInstructions(0), 0);


	for(loc_i = 0; loc_i < locNumFitHists; loc_i++)
		cout << "  locBetaGamma[" << loc_i << "] = " << locXArray[loc_i] << ";";
	cout << endl;

	for(loc_i = 0; loc_i < locNumFitHists; loc_i++)
		cout << "  locdEdxMean[" << loc_i << "] = " << locMeanArray[loc_i] << ";";
	cout << endl;

	for(loc_i = 0; loc_i < locNumFitHists; loc_i++)
		cout << "  locdEdxSigma[" << loc_i << "] = " << locSigmaArray[loc_i] << ";";
	cout << endl;


	locGraphErrors_Mean = new TGraphErrors(locNumFitHists, locXArray, locMeanArray, locXUncertaintyArray, locMeanUncertaintyArray);
	locGraphName.str("");
	locGraphName << loc2FHist->GetName() << "_Means";
	locGraphErrors_Mean->SetName(locGraphName.str().c_str());
	locGraphErrors_Mean->SetTitle(loc2FHist->GetTitle());
	locGraphErrors_Mean->GetXaxis()->SetTitle(loc2FHist->GetXaxis()->GetTitle());
	locGraphTitle.str("");
	locGraphTitle << loc2FHist->GetYaxis()->GetTitle() << " #mu";
	locGraphErrors_Mean->GetYaxis()->SetTitle(locGraphTitle.str().c_str());

	locGraphErrors_Sigma = new TGraphErrors(locNumFitHists, locXArray, locSigmaArray, locXUncertaintyArray, locSigmaUncertaintyArray);
	locGraphName.str("");
	locGraphName << loc2FHist->GetName() << "_Sigmas";
	locGraphErrors_Sigma->SetName(locGraphName.str().c_str());
	locGraphErrors_Sigma->SetTitle(loc2FHist->GetTitle());
	locGraphErrors_Sigma->GetXaxis()->SetTitle(loc2FHist->GetXaxis()->GetTitle());
	locGraphTitle.str("");
	locGraphTitle << loc2FHist->GetYaxis()->GetTitle() << " #sigma";
	locGraphErrors_Sigma->GetYaxis()->SetTitle(locGraphTitle.str().c_str());

	locGraph = new TGraph(locNumFitHists, locXArray, locChiSqPerNDFArray);
	locGraphName.str("");
	locGraphName << loc2FHist->GetName() << "_ChiSqPerNDF";
	locGraph->SetName(locGraphName.str().c_str());
	locGraph->SetTitle(loc2FHist->GetTitle());
	locGraph->GetXaxis()->SetTitle(loc2FHist->GetXaxis()->GetTitle());
	locGraph->GetYaxis()->SetTitle("Fit #chi^{2}/NDF");
//	if((gDrawFlag == 1) || (gDrawFlag == 2))
//		DrawAndSave_ddhHistOrGraph<TGraph>(locGraph, Get_ddhTemplateDrawInstructions(1), 1);
}

void Draw_ProtonPiPlusOverlap(){
	Setup_Routines();

	double *locBetaGamma, *locPArray, *locdEdxMean, *locdEdxSigma, *locXUncertaintyArray;
	unsigned int locNumPoints;
	double locMass;
	TGraphErrors *locGraphErrors;
	TObjArray *locVsPArray = new TObjArray(6);
	TObjArray *locVsBetaGammaArray = new TObjArray(6);


	DecayDrawInstructions* locDrawInstructions;
	vector<int> locHistColors(6, 1);
	vector<string> locLegendNames(6, "");

	locHistColors[0] = 1;  locHistColors[1] = 2;  locHistColors[2] = 3;  locHistColors[3] = 4;  locHistColors[4] = 5;  locHistColors[5] = 6;
	locLegendNames[0] = "FDC Protons";  locLegendNames[1] = "CDC Protons";  locLegendNames[2] = "FDC #pi^{+}'s";  locLegendNames[3] = "CDC #pi^{+}'s";
	locLegendNames[4] = "FDC K^{+}'s";  locLegendNames[5] = "CDC K^{+}'s";
	locDrawInstructions = new DecayDrawInstructions();
	*locDrawInstructions = *(Get_ddhTemplateDrawInstructions(1));
	locDrawInstructions->Set_HistOrGraphColorVector(locHistColors);
	locDrawInstructions->Set_LegendNameVector(locLegendNames);
	locDrawInstructions->Set_DrawOption("p");

	//PROTON, FDC
	locNumPoints = 47;
	locMass = 0.93827203;
	locBetaGamma = new double[locNumPoints];
	locdEdxMean = new double[locNumPoints];
	locdEdxSigma = new double[locNumPoints];
	locPArray = new double[locNumPoints];
	locXUncertaintyArray = new double[locNumPoints];
	locBetaGamma[0] = 0.3;  locBetaGamma[1] = 0.34;  locBetaGamma[2] = 0.38;  locBetaGamma[3] = 0.42;  locBetaGamma[4] = 0.46;  locBetaGamma[5] = 0.5;  locBetaGamma[6] = 0.54;  locBetaGamma[7] = 0.58;  locBetaGamma[8] = 0.62;  locBetaGamma[9] = 0.66;  locBetaGamma[10] = 0.7;  locBetaGamma[11] = 0.74;  locBetaGamma[12] = 0.78;  locBetaGamma[13] = 0.82;  locBetaGamma[14] = 0.86;  locBetaGamma[15] = 0.9;  locBetaGamma[16] = 0.94;  locBetaGamma[17] = 0.98;  locBetaGamma[18] = 1.02;  locBetaGamma[19] = 1.06;  locBetaGamma[20] = 1.1;  locBetaGamma[21] = 1.14;  locBetaGamma[22] = 1.18;  locBetaGamma[23] = 1.22;  locBetaGamma[24] = 1.26;  locBetaGamma[25] = 1.3;  locBetaGamma[26] = 1.34;  locBetaGamma[27] = 1.38;  locBetaGamma[28] = 1.42;  locBetaGamma[29] = 1.46;  locBetaGamma[30] = 1.5;  locBetaGamma[31] = 1.54;  locBetaGamma[32] = 1.58;  locBetaGamma[33] = 1.62;  locBetaGamma[34] = 1.66;  locBetaGamma[35] = 1.7;  locBetaGamma[36] = 1.74;  locBetaGamma[37] = 1.78;  locBetaGamma[38] = 1.82;  locBetaGamma[39] = 1.86;  locBetaGamma[40] = 1.9;  locBetaGamma[41] = 1.94;  locBetaGamma[42] = 1.98;  locBetaGamma[43] = 2.02;  locBetaGamma[44] = 2.06;  locBetaGamma[45] = 2.1;  locBetaGamma[46] = 2.14;
	locdEdxMean[0] = 14.3659;  locdEdxMean[1] = 11.6386;  locdEdxMean[2] = 10.1905;  locdEdxMean[3] = 9.04111;  locdEdxMean[4] = 7.78774;  locdEdxMean[5] = 6.78761;  locdEdxMean[6] = 5.8374;  locdEdxMean[7] = 5.06118;  locdEdxMean[8] = 4.4717;  locdEdxMean[9] = 4.0293;  locdEdxMean[10] = 3.64625;  locdEdxMean[11] = 3.36446;  locdEdxMean[12] = 3.37726;  locdEdxMean[13] = 2.91828;  locdEdxMean[14] = 2.72075;  locdEdxMean[15] = 2.57965;  locdEdxMean[16] = 2.46508;  locdEdxMean[17] = 2.35995;  locdEdxMean[18] = 2.29682;  locdEdxMean[19] = 2.22206;  locdEdxMean[20] = 2.14525;  locdEdxMean[21] = 2.07374;  locdEdxMean[22] = 2.01049;  locdEdxMean[23] = 1.96039;  locdEdxMean[24] = 1.91539;  locdEdxMean[25] = 1.86384;  locdEdxMean[26] = 1.82845;  locdEdxMean[27] = 1.79351;  locdEdxMean[28] = 1.75973;  locdEdxMean[29] = 1.73746;  locdEdxMean[30] = 1.71082;  locdEdxMean[31] = 1.68394;  locdEdxMean[32] = 1.65828;  locdEdxMean[33] = 1.64309;  locdEdxMean[34] = 1.62569;  locdEdxMean[35] = 1.61136;  locdEdxMean[36] = 1.59284;  locdEdxMean[37] = 1.58204;  locdEdxMean[38] = 1.56831;  locdEdxMean[39] = 1.55495;  locdEdxMean[40] = 1.54746;  locdEdxMean[41] = 1.53769;  locdEdxMean[42] = 1.52752;  locdEdxMean[43] = 1.52087;  locdEdxMean[44] = 1.51577;  locdEdxMean[45] = 1.50627;  locdEdxMean[46] = 1.50409;
	locdEdxSigma[0] = 2.93573;  locdEdxSigma[1] = 2.33335;  locdEdxSigma[2] = 2.08781;  locdEdxSigma[3] = 1.71181;  locdEdxSigma[4] = 1.50988;  locdEdxSigma[5] = 1.33544;  locdEdxSigma[6] = 1.04385;  locdEdxSigma[7] = 0.850859;  locdEdxSigma[8] = 0.741105;  locdEdxSigma[9] = 0.658413;  locdEdxSigma[10] = 0.579022;  locdEdxSigma[11] = 0.518601;  locdEdxSigma[12] = 2.4537;  locdEdxSigma[13] = 0.437453;  locdEdxSigma[14] = 0.423342;  locdEdxSigma[15] = 0.399433;  locdEdxSigma[16] = 0.374628;  locdEdxSigma[17] = 0.360932;  locdEdxSigma[18] = 0.346619;  locdEdxSigma[19] = 0.3374;  locdEdxSigma[20] = 0.328811;  locdEdxSigma[21] = 0.297946;  locdEdxSigma[22] = 0.286214;  locdEdxSigma[23] = 0.277289;  locdEdxSigma[24] = 0.273103;  locdEdxSigma[25] = 0.263827;  locdEdxSigma[26] = 0.255441;  locdEdxSigma[27] = 0.25133;  locdEdxSigma[28] = 0.248591;  locdEdxSigma[29] = 0.243123;  locdEdxSigma[30] = 0.242105;  locdEdxSigma[31] = 0.239135;  locdEdxSigma[32] = 0.232915;  locdEdxSigma[33] = 0.233682;  locdEdxSigma[34] = 0.229795;  locdEdxSigma[35] = 0.227354;  locdEdxSigma[36] = 0.225178;  locdEdxSigma[37] = 0.224466;  locdEdxSigma[38] = 0.221981;  locdEdxSigma[39] = 0.217982;  locdEdxSigma[40] = 0.216642;  locdEdxSigma[41] = 0.217325;  locdEdxSigma[42] = 0.216148;  locdEdxSigma[43] = 0.215665;  locdEdxSigma[44] = 0.214442;  locdEdxSigma[45] = 0.213614;  locdEdxSigma[46] = 0.217714;
	for(unsigned int loc_i = 0; loc_i < locNumPoints; loc_i++){
		locXUncertaintyArray[loc_i] = 0.0;
		locPArray[loc_i] = locMass*locBetaGamma[loc_i];
	}
	locGraphErrors = new TGraphErrors(locNumPoints, locBetaGamma, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("#beta#gamma");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsBetaGammaArray->AddAt(locGraphErrors, 0);
	locGraphErrors = new TGraphErrors(locNumPoints, locPArray, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("p (GeV/c)");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsPArray->AddAt(locGraphErrors, 0);

	//PROTON, CDC
	locNumPoints = 47;
	locMass = 0.93827203;
	locBetaGamma = new double[locNumPoints];
	locdEdxMean = new double[locNumPoints];
	locdEdxSigma = new double[locNumPoints];
	locPArray = new double[locNumPoints];
	locXUncertaintyArray = new double[locNumPoints];
	locBetaGamma[0] = 0.3;  locBetaGamma[1] = 0.34;  locBetaGamma[2] = 0.38;  locBetaGamma[3] = 0.42;  locBetaGamma[4] = 0.46;  locBetaGamma[5] = 0.5;  locBetaGamma[6] = 0.54;  locBetaGamma[7] = 0.58;  locBetaGamma[8] = 0.62;  locBetaGamma[9] = 0.66;  locBetaGamma[10] = 0.7;  locBetaGamma[11] = 0.74;  locBetaGamma[12] = 0.78;  locBetaGamma[13] = 0.82;  locBetaGamma[14] = 0.86;  locBetaGamma[15] = 0.9;  locBetaGamma[16] = 0.94;  locBetaGamma[17] = 0.98;  locBetaGamma[18] = 1.02;  locBetaGamma[19] = 1.06;  locBetaGamma[20] = 1.1;  locBetaGamma[21] = 1.14;  locBetaGamma[22] = 1.18;  locBetaGamma[23] = 1.22;  locBetaGamma[24] = 1.26;  locBetaGamma[25] = 1.3;  locBetaGamma[26] = 1.34;  locBetaGamma[27] = 1.38;  locBetaGamma[28] = 1.42;  locBetaGamma[29] = 1.46;  locBetaGamma[30] = 1.5;  locBetaGamma[31] = 1.54;  locBetaGamma[32] = 1.58;  locBetaGamma[33] = 1.62;  locBetaGamma[34] = 1.66;  locBetaGamma[35] = 1.7;  locBetaGamma[36] = 1.74;  locBetaGamma[37] = 1.78;  locBetaGamma[38] = 1.82;  locBetaGamma[39] = 1.86;  locBetaGamma[40] = 1.9;  locBetaGamma[41] = 1.94;  locBetaGamma[42] = 1.98;  locBetaGamma[43] = 2.02;  locBetaGamma[44] = 2.06;  locBetaGamma[45] = 2.1;  locBetaGamma[46] = 2.14;
	locdEdxMean[0] = 13.2769;  locdEdxMean[1] = 11.1142;  locdEdxMean[2] = 10.0385;  locdEdxMean[3] = 9.13612;  locdEdxMean[4] = 8.09766;  locdEdxMean[5] = 6.75607;  locdEdxMean[6] = 5.93831;  locdEdxMean[7] = 5.2655;  locdEdxMean[8] = 4.75216;  locdEdxMean[9] = 4.30025;  locdEdxMean[10] = 3.96891;  locdEdxMean[11] = 3.65299;  locdEdxMean[12] = 3.41423;  locdEdxMean[13] = 3.19939;  locdEdxMean[14] = 3.00491;  locdEdxMean[15] = 2.8645;  locdEdxMean[16] = 2.71925;  locdEdxMean[17] = 2.58784;  locdEdxMean[18] = 2.49657;  locdEdxMean[19] = 2.40005;  locdEdxMean[20] = 2.30766;  locdEdxMean[21] = 2.23464;  locdEdxMean[22] = 2.17738;  locdEdxMean[23] = 2.11424;  locdEdxMean[24] = 2.05425;  locdEdxMean[25] = 2.00659;  locdEdxMean[26] = 1.97077;  locdEdxMean[27] = 1.93114;  locdEdxMean[28] = 1.89473;  locdEdxMean[29] = 1.85834;  locdEdxMean[30] = 1.82908;  locdEdxMean[31] = 1.80608;  locdEdxMean[32] = 1.78316;  locdEdxMean[33] = 1.75923;  locdEdxMean[34] = 1.73596;  locdEdxMean[35] = 1.71507;  locdEdxMean[36] = 1.69805;  locdEdxMean[37] = 1.68427;  locdEdxMean[38] = 1.67062;  locdEdxMean[39] = 1.65653;  locdEdxMean[40] = 1.64473;  locdEdxMean[41] = 1.63283;  locdEdxMean[42] = 1.62268;  locdEdxMean[43] = 1.61341;  locdEdxMean[44] = 1.60417;  locdEdxMean[45] = 1.59654;  locdEdxMean[46] = 1.59605;
	locdEdxSigma[0] = 3.19788;  locdEdxSigma[1] = 2.53706;  locdEdxSigma[2] = 2.02203;  locdEdxSigma[3] = 1.31974;  locdEdxSigma[4] = 1.1821;  locdEdxSigma[5] = 0.916922;  locdEdxSigma[6] = 0.775039;  locdEdxSigma[7] = 0.6374;  locdEdxSigma[8] = 0.547259;  locdEdxSigma[9] = 0.473546;  locdEdxSigma[10] = 0.426132;  locdEdxSigma[11] = 0.386837;  locdEdxSigma[12] = 0.351433;  locdEdxSigma[13] = 0.341537;  locdEdxSigma[14] = 0.318642;  locdEdxSigma[15] = 0.306191;  locdEdxSigma[16] = 0.29701;  locdEdxSigma[17] = 0.282137;  locdEdxSigma[18] = 0.274505;  locdEdxSigma[19] = 0.267861;  locdEdxSigma[20] = 0.260041;  locdEdxSigma[21] = 0.235912;  locdEdxSigma[22] = 0.23393;  locdEdxSigma[23] = 0.228528;  locdEdxSigma[24] = 0.225075;  locdEdxSigma[25] = 0.216987;  locdEdxSigma[26] = 0.214381;  locdEdxSigma[27] = 0.211059;  locdEdxSigma[28] = 0.209439;  locdEdxSigma[29] = 0.205865;  locdEdxSigma[30] = 0.202402;  locdEdxSigma[31] = 0.200022;  locdEdxSigma[32] = 0.198706;  locdEdxSigma[33] = 0.195332;  locdEdxSigma[34] = 0.191766;  locdEdxSigma[35] = 0.191046;  locdEdxSigma[36] = 0.187931;  locdEdxSigma[37] = 0.187142;  locdEdxSigma[38] = 0.186757;  locdEdxSigma[39] = 0.184169;  locdEdxSigma[40] = 0.182682;  locdEdxSigma[41] = 0.181118;  locdEdxSigma[42] = 0.179834;  locdEdxSigma[43] = 0.177071;  locdEdxSigma[44] = 0.175781;  locdEdxSigma[45] = 0.175693;  locdEdxSigma[46] = 0.18061;
	for(unsigned int loc_i = 0; loc_i < locNumPoints; loc_i++){
		locXUncertaintyArray[loc_i] = 0.0;
		locPArray[loc_i] = locMass*locBetaGamma[loc_i];
	}
	locGraphErrors = new TGraphErrors(locNumPoints, locBetaGamma, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("#beta#gamma");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsBetaGammaArray->AddAt(locGraphErrors, 1);
	locGraphErrors = new TGraphErrors(locNumPoints, locPArray, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("p (GeV/c)");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsPArray->AddAt(locGraphErrors, 1);

	//PIPLUS, FDC
	locNumPoints = 44;
	locMass = 0.1395700;
	locBetaGamma = new double[locNumPoints];
	locdEdxMean = new double[locNumPoints];
	locdEdxSigma = new double[locNumPoints];
	locPArray = new double[locNumPoints];
	locXUncertaintyArray = new double[locNumPoints];
	locBetaGamma[0] = 0.48;  locBetaGamma[1] = 0.8;  locBetaGamma[2] = 1.12;  locBetaGamma[3] = 1.44;  locBetaGamma[4] = 1.76;  locBetaGamma[5] = 2.08;  locBetaGamma[6] = 2.4;  locBetaGamma[7] = 2.72;  locBetaGamma[8] = 3.04;  locBetaGamma[9] = 3.36;  locBetaGamma[10] = 3.68;  locBetaGamma[11] = 4;  locBetaGamma[12] = 4.32;  locBetaGamma[13] = 4.64;  locBetaGamma[14] = 4.96;  locBetaGamma[15] = 5.28;  locBetaGamma[16] = 5.6;  locBetaGamma[17] = 5.92;  locBetaGamma[18] = 6.24;  locBetaGamma[19] = 6.56;  locBetaGamma[20] = 6.88;  locBetaGamma[21] = 7.2;  locBetaGamma[22] = 7.52;  locBetaGamma[23] = 7.84;  locBetaGamma[24] = 8.16;  locBetaGamma[25] = 8.48;  locBetaGamma[26] = 8.8;  locBetaGamma[27] = 9.12;  locBetaGamma[28] = 9.44;  locBetaGamma[29] = 9.76;  locBetaGamma[30] = 10.08;  locBetaGamma[31] = 10.4;  locBetaGamma[32] = 10.72;  locBetaGamma[33] = 11.04;  locBetaGamma[34] = 11.36;  locBetaGamma[35] = 11.68;  locBetaGamma[36] = 12;  locBetaGamma[37] = 12.32;  locBetaGamma[38] = 12.64;  locBetaGamma[39] = 12.96;  locBetaGamma[40] = 13.28;  locBetaGamma[41] = 13.6;  locBetaGamma[42] = 13.92;  locBetaGamma[43] = 14.24;
	locdEdxMean[0] = 1.32105;  locdEdxMean[1] = 1.50904;  locdEdxMean[2] = 1.50873;  locdEdxMean[3] = 1.57426;  locdEdxMean[4] = 1.52213;  locdEdxMean[5] = 1.45914;  locdEdxMean[6] = 1.4271;  locdEdxMean[7] = 1.42623;  locdEdxMean[8] = 1.42268;  locdEdxMean[9] = 1.42632;  locdEdxMean[10] = 1.43233;  locdEdxMean[11] = 1.43742;  locdEdxMean[12] = 1.4467;  locdEdxMean[13] = 1.45642;  locdEdxMean[14] = 1.46337;  locdEdxMean[15] = 1.4719;  locdEdxMean[16] = 1.48252;  locdEdxMean[17] = 1.49443;  locdEdxMean[18] = 1.50234;  locdEdxMean[19] = 1.513;  locdEdxMean[20] = 1.52394;  locdEdxMean[21] = 1.5289;  locdEdxMean[22] = 1.53878;  locdEdxMean[23] = 1.54396;  locdEdxMean[24] = 1.55427;  locdEdxMean[25] = 1.56105;  locdEdxMean[26] = 1.5674;  locdEdxMean[27] = 1.57221;  locdEdxMean[28] = 1.58027;  locdEdxMean[29] = 1.5865;  locdEdxMean[30] = 1.59519;  locdEdxMean[31] = 1.60211;  locdEdxMean[32] = 1.60751;  locdEdxMean[33] = 1.61522;  locdEdxMean[34] = 1.62066;  locdEdxMean[35] = 1.6277;  locdEdxMean[36] = 1.63227;  locdEdxMean[37] = 1.63661;  locdEdxMean[38] = 1.64457;  locdEdxMean[39] = 1.65043;  locdEdxMean[40] = 1.6557;  locdEdxMean[41] = 1.66085;  locdEdxMean[42] = 1.66486;  locdEdxMean[43] = 1.66916;
	locdEdxSigma[0] = 0.384194;  locdEdxSigma[1] = 0.481736;  locdEdxSigma[2] = 0.408882;  locdEdxSigma[3] = 0.308842;  locdEdxSigma[4] = 0.238201;  locdEdxSigma[5] = 0.205383;  locdEdxSigma[6] = 0.199131;  locdEdxSigma[7] = 0.183296;  locdEdxSigma[8] = 0.17939;  locdEdxSigma[9] = 0.180875;  locdEdxSigma[10] = 0.180979;  locdEdxSigma[11] = 0.181938;  locdEdxSigma[12] = 0.185807;  locdEdxSigma[13] = 0.186497;  locdEdxSigma[14] = 0.188448;  locdEdxSigma[15] = 0.190472;  locdEdxSigma[16] = 0.189878;  locdEdxSigma[17] = 0.189112;  locdEdxSigma[18] = 0.189443;  locdEdxSigma[19] = 0.192163;  locdEdxSigma[20] = 0.19386;  locdEdxSigma[21] = 0.194046;  locdEdxSigma[22] = 0.194451;  locdEdxSigma[23] = 0.194133;  locdEdxSigma[24] = 0.195598;  locdEdxSigma[25] = 0.198318;  locdEdxSigma[26] = 0.198263;  locdEdxSigma[27] = 0.198669;  locdEdxSigma[28] = 0.19729;  locdEdxSigma[29] = 0.197867;  locdEdxSigma[30] = 0.199663;  locdEdxSigma[31] = 0.201134;  locdEdxSigma[32] = 0.201602;  locdEdxSigma[33] = 0.202371;  locdEdxSigma[34] = 0.202888;  locdEdxSigma[35] = 0.201155;  locdEdxSigma[36] = 0.205072;  locdEdxSigma[37] = 0.201526;  locdEdxSigma[38] = 0.206744;  locdEdxSigma[39] = 0.209571;  locdEdxSigma[40] = 0.207218;  locdEdxSigma[41] = 0.204573;  locdEdxSigma[42] = 0.20254;  locdEdxSigma[43] = 0.208395;
	for(unsigned int loc_i = 0; loc_i < locNumPoints; loc_i++){
		locXUncertaintyArray[loc_i] = 0.0;
		locPArray[loc_i] = locMass*locBetaGamma[loc_i];
	}
	locGraphErrors = new TGraphErrors(locNumPoints, locBetaGamma, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("#beta#gamma");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsBetaGammaArray->AddAt(locGraphErrors, 2);
	locGraphErrors = new TGraphErrors(locNumPoints, locPArray, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("p (GeV/c)");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsPArray->AddAt(locGraphErrors, 2);

	//PIPLUS, CDC
	locNumPoints = 44;
	locMass = 0.1395700;
	locBetaGamma = new double[locNumPoints];
	locdEdxMean = new double[locNumPoints];
	locdEdxSigma = new double[locNumPoints];
	locPArray = new double[locNumPoints];
	locXUncertaintyArray = new double[locNumPoints];
	locBetaGamma[0] = 0.48;  locBetaGamma[1] = 0.8;  locBetaGamma[2] = 1.12;  locBetaGamma[3] = 1.44;  locBetaGamma[4] = 1.76;  locBetaGamma[5] = 2.08;  locBetaGamma[6] = 2.4;  locBetaGamma[7] = 2.72;  locBetaGamma[8] = 3.04;  locBetaGamma[9] = 3.36;  locBetaGamma[10] = 3.68;  locBetaGamma[11] = 4;  locBetaGamma[12] = 4.32;  locBetaGamma[13] = 4.64;  locBetaGamma[14] = 4.96;  locBetaGamma[15] = 5.28;  locBetaGamma[16] = 5.6;  locBetaGamma[17] = 5.92;  locBetaGamma[18] = 6.24;  locBetaGamma[19] = 6.56;  locBetaGamma[20] = 6.88;  locBetaGamma[21] = 7.2;  locBetaGamma[22] = 7.52;  locBetaGamma[23] = 7.84;  locBetaGamma[24] = 8.16;  locBetaGamma[25] = 8.48;  locBetaGamma[26] = 8.8;  locBetaGamma[27] = 9.12;  locBetaGamma[28] = 9.44;  locBetaGamma[29] = 9.76;  locBetaGamma[30] = 10.08;  locBetaGamma[31] = 10.4;  locBetaGamma[32] = 10.72;  locBetaGamma[33] = 11.04;  locBetaGamma[34] = 11.36;  locBetaGamma[35] = 11.68;  locBetaGamma[36] = 12;  locBetaGamma[37] = 12.32;  locBetaGamma[38] = 12.64;  locBetaGamma[39] = 12.96;  locBetaGamma[40] = 13.28;  locBetaGamma[41] = 13.6;  locBetaGamma[42] = 13.92;  locBetaGamma[43] = 14.24;
	locdEdxMean[0] = 2.40705;  locdEdxMean[1] = 2.64016;  locdEdxMean[2] = 1.97719;  locdEdxMean[3] = 1.74434;  locdEdxMean[4] = 1.60673;  locdEdxMean[5] = 1.52877;  locdEdxMean[6] = 1.47792;  locdEdxMean[7] = 1.44914;  locdEdxMean[8] = 1.44456;  locdEdxMean[9] = 1.44089;  locdEdxMean[10] = 1.44077;  locdEdxMean[11] = 1.44643;  locdEdxMean[12] = 1.46221;  locdEdxMean[13] = 1.47146;  locdEdxMean[14] = 1.48451;  locdEdxMean[15] = 1.49671;  locdEdxMean[16] = 1.50718;  locdEdxMean[17] = 1.51892;  locdEdxMean[18] = 1.53035;  locdEdxMean[19] = 1.54224;  locdEdxMean[20] = 1.5531;  locdEdxMean[21] = 1.56325;  locdEdxMean[22] = 1.57347;  locdEdxMean[23] = 1.58422;  locdEdxMean[24] = 1.59239;  locdEdxMean[25] = 1.60057;  locdEdxMean[26] = 1.61224;  locdEdxMean[27] = 1.62025;  locdEdxMean[28] = 1.6288;  locdEdxMean[29] = 1.63748;  locdEdxMean[30] = 1.64582;  locdEdxMean[31] = 1.65254;  locdEdxMean[32] = 1.65864;  locdEdxMean[33] = 1.66693;  locdEdxMean[34] = 1.67503;  locdEdxMean[35] = 1.68328;  locdEdxMean[36] = 1.68964;  locdEdxMean[37] = 1.69732;  locdEdxMean[38] = 1.70258;  locdEdxMean[39] = 1.70849;  locdEdxMean[40] = 1.71469;  locdEdxMean[41] = 1.72359;  locdEdxMean[42] = 1.73032;  locdEdxMean[43] = 1.73672;
	locdEdxSigma[0] = 0.00923813;  locdEdxSigma[1] = 0.753794;  locdEdxSigma[2] = 0.461342;  locdEdxSigma[3] = 0.273398;  locdEdxSigma[4] = 0.239905;  locdEdxSigma[5] = 0.23336;  locdEdxSigma[6] = 0.231064;  locdEdxSigma[7] = 0.216753;  locdEdxSigma[8] = 0.205846;  locdEdxSigma[9] = 0.202349;  locdEdxSigma[10] = 0.19949;  locdEdxSigma[11] = 0.19614;  locdEdxSigma[12] = 0.192053;  locdEdxSigma[13] = 0.189132;  locdEdxSigma[14] = 0.18807;  locdEdxSigma[15] = 0.186939;  locdEdxSigma[16] = 0.186402;  locdEdxSigma[17] = 0.18636;  locdEdxSigma[18] = 0.185533;  locdEdxSigma[19] = 0.185677;  locdEdxSigma[20] = 0.184546;  locdEdxSigma[21] = 0.183511;  locdEdxSigma[22] = 0.183899;  locdEdxSigma[23] = 0.183173;  locdEdxSigma[24] = 0.182368;  locdEdxSigma[25] = 0.182132;  locdEdxSigma[26] = 0.183335;  locdEdxSigma[27] = 0.182698;  locdEdxSigma[28] = 0.1839;  locdEdxSigma[29] = 0.184523;  locdEdxSigma[30] = 0.183395;  locdEdxSigma[31] = 0.182982;  locdEdxSigma[32] = 0.181914;  locdEdxSigma[33] = 0.183204;  locdEdxSigma[34] = 0.182819;  locdEdxSigma[35] = 0.184335;  locdEdxSigma[36] = 0.183594;  locdEdxSigma[37] = 0.183738;  locdEdxSigma[38] = 0.181587;  locdEdxSigma[39] = 0.183005;  locdEdxSigma[40] = 0.181546;  locdEdxSigma[41] = 0.183761;  locdEdxSigma[42] = 0.181959;  locdEdxSigma[43] = 0.185626;
	for(unsigned int loc_i = 0; loc_i < locNumPoints; loc_i++){
		locXUncertaintyArray[loc_i] = 0.0;
		locPArray[loc_i] = locMass*locBetaGamma[loc_i];
	}
	locGraphErrors = new TGraphErrors(locNumPoints, locBetaGamma, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("#beta#gamma");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsBetaGammaArray->AddAt(locGraphErrors, 3);
	locGraphErrors = new TGraphErrors(locNumPoints, locPArray, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("p (GeV/c)");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsPArray->AddAt(locGraphErrors, 3);


	//KPLUS, FDC
	locNumPoints = 50;
	locMass = 0.493677;
	locBetaGamma = new double[locNumPoints];
	locdEdxMean = new double[locNumPoints];
	locdEdxSigma = new double[locNumPoints];
	locPArray = new double[locNumPoints];
	locXUncertaintyArray = new double[locNumPoints];
	locBetaGamma[0] = 0.12;  locBetaGamma[1] = 0.2;  locBetaGamma[2] = 0.28;  locBetaGamma[3] = 0.36;  locBetaGamma[4] = 0.44;  locBetaGamma[5] = 0.52;  locBetaGamma[6] = 0.6;  locBetaGamma[7] = 0.68;  locBetaGamma[8] = 0.76;  locBetaGamma[9] = 0.84;  locBetaGamma[10] = 0.92;  locBetaGamma[11] = 1;  locBetaGamma[12] = 1.08;  locBetaGamma[13] = 1.16;  locBetaGamma[14] = 1.24;  locBetaGamma[15] = 1.32;  locBetaGamma[16] = 1.4;  locBetaGamma[17] = 1.48;  locBetaGamma[18] = 1.56;  locBetaGamma[19] = 1.64;  locBetaGamma[20] = 1.72;  locBetaGamma[21] = 1.8;  locBetaGamma[22] = 1.88;  locBetaGamma[23] = 1.96;  locBetaGamma[24] = 2.04;  locBetaGamma[25] = 2.12;  locBetaGamma[26] = 2.2;  locBetaGamma[27] = 2.28;  locBetaGamma[28] = 2.36;  locBetaGamma[29] = 2.44;  locBetaGamma[30] = 2.52;  locBetaGamma[31] = 2.6;  locBetaGamma[32] = 2.68;  locBetaGamma[33] = 2.76;  locBetaGamma[34] = 2.84;  locBetaGamma[35] = 2.92;  locBetaGamma[36] = 3;  locBetaGamma[37] = 3.08;  locBetaGamma[38] = 3.16;  locBetaGamma[39] = 3.24;  locBetaGamma[40] = 3.32;  locBetaGamma[41] = 3.4;  locBetaGamma[42] = 3.48;  locBetaGamma[43] = 3.56;  locBetaGamma[44] = 3.64;  locBetaGamma[45] = 3.72;  locBetaGamma[46] = 3.8;  locBetaGamma[47] = 3.88;  locBetaGamma[48] = 3.96;  locBetaGamma[49] = 4.04;
	locdEdxMean[0] = 1.35012;  locdEdxMean[1] = 1.38705;  locdEdxMean[2] = 1.4757;  locdEdxMean[3] = 1.39719;  locdEdxMean[4] = 1.35641;  locdEdxMean[5] = 1.38561;  locdEdxMean[6] = 1.39104;  locdEdxMean[7] = 1.37202;  locdEdxMean[8] = 1.37453;  locdEdxMean[9] = 1.39048;  locdEdxMean[10] = 2.47808;  locdEdxMean[11] = 2.27568;  locdEdxMean[12] = 2.1479;  locdEdxMean[13] = 2.00315;  locdEdxMean[14] = 1.88838;  locdEdxMean[15] = 1.78769;  locdEdxMean[16] = 1.72742;  locdEdxMean[17] = 1.68317;  locdEdxMean[18] = 1.64201;  locdEdxMean[19] = 1.61215;  locdEdxMean[20] = 1.59323;  locdEdxMean[21] = 1.56803;  locdEdxMean[22] = 1.54727;  locdEdxMean[23] = 1.5315;  locdEdxMean[24] = 1.5186;  locdEdxMean[25] = 1.50544;  locdEdxMean[26] = 1.49661;  locdEdxMean[27] = 1.48962;  locdEdxMean[28] = 1.48178;  locdEdxMean[29] = 1.47437;  locdEdxMean[30] = 1.46332;  locdEdxMean[31] = 1.46012;  locdEdxMean[32] = 1.45291;  locdEdxMean[33] = 1.45176;  locdEdxMean[34] = 1.44809;  locdEdxMean[35] = 1.44591;  locdEdxMean[36] = 1.44005;  locdEdxMean[37] = 1.43917;  locdEdxMean[38] = 1.43704;  locdEdxMean[39] = 1.43816;  locdEdxMean[40] = 1.43037;  locdEdxMean[41] = 1.43408;  locdEdxMean[42] = 1.43274;  locdEdxMean[43] = 1.43441;  locdEdxMean[44] = 1.43;  locdEdxMean[45] = 1.43307;  locdEdxMean[46] = 1.4322;  locdEdxMean[47] = 1.42994;  locdEdxMean[48] = 1.43297;  locdEdxMean[49] = 1.43424;
	locdEdxSigma[0] = 0.342043;  locdEdxSigma[1] = 0.415316;  locdEdxSigma[2] = 0.469585;  locdEdxSigma[3] = 0.426567;  locdEdxSigma[4] = 0.311762;  locdEdxSigma[5] = 0.276792;  locdEdxSigma[6] = 0.258224;  locdEdxSigma[7] = 0.240704;  locdEdxSigma[8] = 0.234693;  locdEdxSigma[9] = 0.234922;  locdEdxSigma[10] = 0.43263;  locdEdxSigma[11] = 0.371096;  locdEdxSigma[12] = 0.39021;  locdEdxSigma[13] = 0.323567;  locdEdxSigma[14] = 0.344277;  locdEdxSigma[15] = 0.332502;  locdEdxSigma[16] = 0.297985;  locdEdxSigma[17] = 0.266577;  locdEdxSigma[18] = 0.249613;  locdEdxSigma[19] = 0.235291;  locdEdxSigma[20] = 0.227965;  locdEdxSigma[21] = 0.221766;  locdEdxSigma[22] = 0.216916;  locdEdxSigma[23] = 0.210419;  locdEdxSigma[24] = 0.205736;  locdEdxSigma[25] = 0.208651;  locdEdxSigma[26] = 0.206319;  locdEdxSigma[27] = 0.204752;  locdEdxSigma[28] = 0.20364;  locdEdxSigma[29] = 0.200502;  locdEdxSigma[30] = 0.203598;  locdEdxSigma[31] = 0.200443;  locdEdxSigma[32] = 0.19678;  locdEdxSigma[33] = 0.199707;  locdEdxSigma[34] = 0.201759;  locdEdxSigma[35] = 0.203394;  locdEdxSigma[36] = 0.197385;  locdEdxSigma[37] = 0.196626;  locdEdxSigma[38] = 0.198206;  locdEdxSigma[39] = 0.198421;  locdEdxSigma[40] = 0.197648;  locdEdxSigma[41] = 0.198405;  locdEdxSigma[42] = 0.195213;  locdEdxSigma[43] = 0.197544;  locdEdxSigma[44] = 0.196865;  locdEdxSigma[45] = 0.198364;  locdEdxSigma[46] = 0.197348;  locdEdxSigma[47] = 0.197956;  locdEdxSigma[48] = 0.193743;  locdEdxSigma[49] = 0.195151;
	for(unsigned int loc_i = 0; loc_i < locNumPoints; loc_i++){
		locXUncertaintyArray[loc_i] = 0.0;
		locPArray[loc_i] = locMass*locBetaGamma[loc_i];
	}
	locGraphErrors = new TGraphErrors(locNumPoints, locBetaGamma, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("#beta#gamma");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsBetaGammaArray->AddAt(locGraphErrors, 4);
	locGraphErrors = new TGraphErrors(locNumPoints, locPArray, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("p (GeV/c)");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsPArray->AddAt(locGraphErrors, 4);

	//KPLUS, CDC
	locNumPoints = 49;
	locMass = 0.493677;
	locBetaGamma = new double[locNumPoints];
	locdEdxMean = new double[locNumPoints];
	locdEdxSigma = new double[locNumPoints];
	locPArray = new double[locNumPoints];
	locXUncertaintyArray = new double[locNumPoints];
	locBetaGamma[0] = 0.2;  locBetaGamma[1] = 0.28;  locBetaGamma[2] = 0.36;  locBetaGamma[3] = 0.44;  locBetaGamma[4] = 0.52;  locBetaGamma[5] = 0.6;  locBetaGamma[6] = 0.68;  locBetaGamma[7] = 0.76;  locBetaGamma[8] = 0.84;  locBetaGamma[9] = 0.92;  locBetaGamma[10] = 1;  locBetaGamma[11] = 1.08;  locBetaGamma[12] = 1.16;  locBetaGamma[13] = 1.24;  locBetaGamma[14] = 1.32;  locBetaGamma[15] = 1.4;  locBetaGamma[16] = 1.48;  locBetaGamma[17] = 1.56;  locBetaGamma[18] = 1.64;  locBetaGamma[19] = 1.72;  locBetaGamma[20] = 1.8;  locBetaGamma[21] = 1.88;  locBetaGamma[22] = 1.96;  locBetaGamma[23] = 2.04;  locBetaGamma[24] = 2.12;  locBetaGamma[25] = 2.2;  locBetaGamma[26] = 2.28;  locBetaGamma[27] = 2.36;  locBetaGamma[28] = 2.44;  locBetaGamma[29] = 2.52;  locBetaGamma[30] = 2.6;  locBetaGamma[31] = 2.68;  locBetaGamma[32] = 2.76;  locBetaGamma[33] = 2.84;  locBetaGamma[34] = 2.92;  locBetaGamma[35] = 3;  locBetaGamma[36] = 3.08;  locBetaGamma[37] = 3.16;  locBetaGamma[38] = 3.24;  locBetaGamma[39] = 3.32;  locBetaGamma[40] = 3.4;  locBetaGamma[41] = 3.48;  locBetaGamma[42] = 3.56;  locBetaGamma[43] = 3.64;  locBetaGamma[44] = 3.72;  locBetaGamma[45] = 3.8;  locBetaGamma[46] = 3.88;  locBetaGamma[47] = 3.96;  locBetaGamma[48] = 4.04;
	locdEdxMean[0] = 3.2117;  locdEdxMean[1] = 1.92621;  locdEdxMean[2] = 1.63171;  locdEdxMean[3] = 1.54946;  locdEdxMean[4] = 1.43884;  locdEdxMean[5] = 1.43401;  locdEdxMean[6] = 4.03109;  locdEdxMean[7] = 3.45989;  locdEdxMean[8] = 3.04638;  locdEdxMean[9] = 2.74302;  locdEdxMean[10] = 2.49407;  locdEdxMean[11] = 2.30716;  locdEdxMean[12] = 2.15957;  locdEdxMean[13] = 2.03802;  locdEdxMean[14] = 1.95242;  locdEdxMean[15] = 1.87438;  locdEdxMean[16] = 1.80751;  locdEdxMean[17] = 1.76246;  locdEdxMean[18] = 1.71855;  locdEdxMean[19] = 1.67958;  locdEdxMean[20] = 1.65377;  locdEdxMean[21] = 1.62775;  locdEdxMean[22] = 1.60503;  locdEdxMean[23] = 1.58512;  locdEdxMean[24] = 1.57129;  locdEdxMean[25] = 1.55657;  locdEdxMean[26] = 1.54396;  locdEdxMean[27] = 1.53245;  locdEdxMean[28] = 1.52701;  locdEdxMean[29] = 1.5202;  locdEdxMean[30] = 1.51276;  locdEdxMean[31] = 1.50879;  locdEdxMean[32] = 1.5018;  locdEdxMean[33] = 1.49804;  locdEdxMean[34] = 1.49589;  locdEdxMean[35] = 1.494;  locdEdxMean[36] = 1.49093;  locdEdxMean[37] = 1.4921;  locdEdxMean[38] = 1.49017;  locdEdxMean[39] = 1.48867;  locdEdxMean[40] = 1.48755;  locdEdxMean[41] = 1.48844;  locdEdxMean[42] = 1.48861;  locdEdxMean[43] = 1.48829;  locdEdxMean[44] = 1.49078;  locdEdxMean[45] = 1.49107;  locdEdxMean[46] = 1.49341;  locdEdxMean[47] = 1.49496;  locdEdxMean[48] = 1.49551;
	locdEdxSigma[0] = 0.00990457;  locdEdxSigma[1] = 0.713109;  locdEdxSigma[2] = 0.422844;  locdEdxSigma[3] = 0.319996;  locdEdxSigma[4] = 0.228746;  locdEdxSigma[5] = 0.299985;  locdEdxSigma[6] = 0.534387;  locdEdxSigma[7] = 0.446107;  locdEdxSigma[8] = 0.381853;  locdEdxSigma[9] = 0.334705;  locdEdxSigma[10] = 0.292848;  locdEdxSigma[11] = 0.279025;  locdEdxSigma[12] = 0.256161;  locdEdxSigma[13] = 0.246785;  locdEdxSigma[14] = 0.233125;  locdEdxSigma[15] = 0.227021;  locdEdxSigma[16] = 0.218084;  locdEdxSigma[17] = 0.210811;  locdEdxSigma[18] = 0.205133;  locdEdxSigma[19] = 0.198922;  locdEdxSigma[20] = 0.195733;  locdEdxSigma[21] = 0.191139;  locdEdxSigma[22] = 0.188482;  locdEdxSigma[23] = 0.184758;  locdEdxSigma[24] = 0.183807;  locdEdxSigma[25] = 0.18163;  locdEdxSigma[26] = 0.179246;  locdEdxSigma[27] = 0.178835;  locdEdxSigma[28] = 0.178264;  locdEdxSigma[29] = 0.176596;  locdEdxSigma[30] = 0.174917;  locdEdxSigma[31] = 0.174964;  locdEdxSigma[32] = 0.172598;  locdEdxSigma[33] = 0.171251;  locdEdxSigma[34] = 0.172524;  locdEdxSigma[35] = 0.171673;  locdEdxSigma[36] = 0.171231;  locdEdxSigma[37] = 0.170838;  locdEdxSigma[38] = 0.16876;  locdEdxSigma[39] = 0.167956;  locdEdxSigma[40] = 0.167727;  locdEdxSigma[41] = 0.168108;  locdEdxSigma[42] = 0.167474;  locdEdxSigma[43] = 0.16609;  locdEdxSigma[44] = 0.166928;  locdEdxSigma[45] = 0.164927;  locdEdxSigma[46] = 0.166005;  locdEdxSigma[47] = 0.166381;  locdEdxSigma[48] = 0.165866;
	for(unsigned int loc_i = 0; loc_i < locNumPoints; loc_i++){
		locXUncertaintyArray[loc_i] = 0.0;
		locPArray[loc_i] = locMass*locBetaGamma[loc_i];
	}
	locGraphErrors = new TGraphErrors(locNumPoints, locBetaGamma, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("#beta#gamma");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsBetaGammaArray->AddAt(locGraphErrors, 5);
	locGraphErrors = new TGraphErrors(locNumPoints, locPArray, locdEdxMean, locXUncertaintyArray, locdEdxSigma);
	locGraphErrors->GetXaxis()->SetTitle("p (GeV/c)");
	locGraphErrors->GetYaxis()->SetTitle("#frac{dE}{dx} #mu");
	locVsPArray->AddAt(locGraphErrors, 5);

	DrawAndSave_ddhHistOrGraphSuperimposition<TGraphErrors>(locVsBetaGammaArray, locDrawInstructions, 1);
	DrawAndSave_ddhHistOrGraphSuperimposition<TGraphErrors>(locVsPArray, locDrawInstructions, 1);
}

