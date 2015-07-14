#define DCdEdxSelector_cxx
// The class definition in DCdEdxSelector.h has been generated automatically
// by the ROOT utility TTree::MakeSelector(). This class is derived
// from the ROOT class TSelector. For more information on the TSelector
// framework see $ROOTSYS/README/README.SELECTOR or the ROOT User Manual.

// The following methods are defined in this file:
//    Begin():        called every time a loop on the tree starts,
//                    a convenient place to create your histograms.
//    SlaveBegin():   called after Begin(), when on PROOF called only on the
//                    slave servers.
//    Process():      called for each event, in this function you decide what
//                    to read and fill your histograms.
//    SlaveTerminate: called at the end of the loop on the tree, when on PROOF
//                    called only on the slave servers.
//    Terminate():    called at the end of the loop on the tree,
//                    a convenient place to draw/fit your histograms.
//
// To use this file, try the following session on your Tree T:
//
// Root > T->Process("DCdEdxSelector.C")
// Root > T->Process("DCdEdxSelector.C","some options")
// Root > T->Process("DCdEdxSelector.C+")
//

#include "DCdEdxSelector.h"
#include <TH2.h>
#include <TStyle.h>

void DCdEdxSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
	TF1 *locFunc;

	dRhoZoverA_CDC = 1.29409e-07/0.1535E-3;
	dRhoZoverA_FDC = 1.3277e-07/0.1535E-3;

	dCalcFOMManuallyFlag = true;
//	dParticleName = "Proton";
//	dParticleName = "PiPlus";
	dParticleName = "KPlus";

	dOutputFile = new TFile("dh_DCdEdxHists.root", "RECREATE");

   ddEdxMeanFunc_FDC_Proton = new TF1("dPID_dEdxMeanFunc_FDC_Proton", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_FDC_Proton->SetParameters(47.5481, -4.70207, 2.22134, -0.357701);
	
   ddEdxMeanFunc_CDC_Proton = new TF1("dPID_dEdxMeanFunc_CDC_Proton", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_CDC_Proton->SetParameters(40.3217, -4.33611, 2.36034, -0.381628);

   ddEdxMeanFunc_FDC_KPlus = new TF1("dPID_dEdxMeanFunc_FDC_KPlus", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_FDC_KPlus->SetParameters(11.1221, -2.65063, 1.50373, -0.0202496);
	
   ddEdxMeanFunc_CDC_KPlus = new TF1("dPID_dEdxMeanFunc_CDC_KPlus", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   ddEdxMeanFunc_CDC_KPlus->SetParameters(12.3025, -2.51467, 1.53035, -0.0120788);

	dSigmaFuncArray_FDC_Proton = new TObjArray(4);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(15.0407, -4.60489, 0.501318, -0.0913598);
	dSigmaFuncArray_FDC_Proton->AddAt(locFunc, 0);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(26.6243, -5.92033, 0.384732, -0.0820965);
	dSigmaFuncArray_FDC_Proton->AddAt(locFunc, 1);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.0415, -4.956, 0.29667, -0.0551989);
	dSigmaFuncArray_FDC_Proton->AddAt(locFunc, 2);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_Proton4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(14.5991, -5.9726, 0.274239, -0.0546458);
	dSigmaFuncArray_FDC_Proton->AddAt(locFunc, 3);
	ddEdxSigmaNumHitsVector_FDC_Proton.resize(4);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_Proton[loc_i] = 3 + 3*loc_i;

	dSigmaFuncArray_CDC_Proton = new TObjArray(6);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(14.263, -3.91119, 0.471623, -0.0867709);
	dSigmaFuncArray_CDC_Proton->AddAt(locFunc, 0);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.5091, -4.2701, 0.638759, -0.191789);
	dSigmaFuncArray_CDC_Proton->AddAt(locFunc, 1);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(5.74486, -3.46166, 0.353456, -0.0594251);
	dSigmaFuncArray_CDC_Proton->AddAt(locFunc, 2);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(6.93564, -4.55521, 0.372744, -0.0809587);
	dSigmaFuncArray_CDC_Proton->AddAt(locFunc, 3);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton5", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(11.3186, -6.17783, 0.269247, -0.0525885);
	dSigmaFuncArray_CDC_Proton->AddAt(locFunc, 4);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_Proton6", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(8.41072, -5.90604, 0.23806, -0.0406181);
	dSigmaFuncArray_CDC_Proton->AddAt(locFunc, 5);
	ddEdxSigmaNumHitsVector_CDC_Proton.resize(6);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_Proton[loc_i] = 4 + 2*loc_i;

	dSigmaFuncArray_CDC_PiPlus = new TObjArray(6);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.8644, -1.92051, 0.401207, -0.00510045);
	dSigmaFuncArray_CDC_PiPlus->AddAt(locFunc, 0);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(3.08596, -2.55591, 0.326127, -0.00439305);
	dSigmaFuncArray_CDC_PiPlus->AddAt(locFunc, 1);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.63183, -2.18999, 0.276187, -0.00277214);
	dSigmaFuncArray_CDC_PiPlus->AddAt(locFunc, 2);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.42353, -1.12885, 0.218266, 0.000475549);
	dSigmaFuncArray_CDC_PiPlus->AddAt(locFunc, 3);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus5", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.599551, -1.27373, 0.15483, 0.00125788);
	dSigmaFuncArray_CDC_PiPlus->AddAt(locFunc, 4);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_PiPlus6", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.02319, -2.15814, 0.146108, 0.00126295);
	dSigmaFuncArray_CDC_PiPlus->AddAt(locFunc, 5);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_PiPlus.push_back(4 + 2*loc_i);

	dSigmaFuncArray_FDC_PiPlus = new TObjArray(4);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.972821, -1.13778, 0.28177, 0.00414764);
	dSigmaFuncArray_FDC_PiPlus->AddAt(locFunc, 0);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.72011, -1.13622, 0.205834, 0.00248226);
	dSigmaFuncArray_FDC_PiPlus->AddAt(locFunc, 1);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.783602, -1.58954, 0.175932, 0.00139313);
	dSigmaFuncArray_FDC_PiPlus->AddAt(locFunc, 2);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_PiPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(0.666107, -1.50611, 0.144197, 0.0019479);
	dSigmaFuncArray_FDC_PiPlus->AddAt(locFunc, 3);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_PiPlus.push_back(3 + 3*loc_i);

	//sigmas, k+
	dSigmaFuncArray_CDC_KPlus = new TObjArray(6);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.3149, -3.48486, 0.590082, -0.0757512);
	dSigmaFuncArray_CDC_KPlus->AddAt(locFunc, 0);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(10.9306, -3.49101, 0.443805, -0.0527794);
	dSigmaFuncArray_CDC_KPlus->AddAt(locFunc, 1);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(7.71149, -3.19237, 0.319949, -0.0248138);
	dSigmaFuncArray_CDC_KPlus->AddAt(locFunc, 2);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(2.70214, -2.51235, 0.224439, -0.00564491);
	dSigmaFuncArray_CDC_KPlus->AddAt(locFunc, 3);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus5", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.12165, -2.16051, 0.155717, 0.000195505);
	dSigmaFuncArray_CDC_KPlus->AddAt(locFunc, 4);
   locFunc = new TF1("dPID_dEdxSigmaFunc_CDC_KPlus6", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.60267, -2.62908, 0.153127, -0.000775554);
	dSigmaFuncArray_CDC_KPlus->AddAt(locFunc, 5);
	for(unsigned int loc_i = 0; loc_i < 6; loc_i++)
		ddEdxSigmaNumHitsVector_CDC_KPlus.push_back(4 + 2*loc_i);

	dSigmaFuncArray_FDC_KPlus = new TObjArray(4);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus1", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(3.4624, -2.51387, 0.328651, -0.00906703);
	dSigmaFuncArray_FDC_KPlus->AddAt(locFunc, 0);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus2", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(8.77765, -3.16076, 0.20651, 0.00175279);
	dSigmaFuncArray_FDC_KPlus->AddAt(locFunc, 1);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus3", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(6.44965, -3.01797, 0.175224, -0.000185835);
	dSigmaFuncArray_FDC_KPlus->AddAt(locFunc, 2);
   locFunc = new TF1("dPID_dEdxSigmaFunc_FDC_KPlus4", "[0]*exp([1]*x) + [2] + [3]*x", 0.0, 100.0);
   locFunc->SetParameters(1.37552, -2.21704, 0.148445, 0.0014621);
	dSigmaFuncArray_FDC_KPlus->AddAt(locFunc, 3);
	for(unsigned int loc_i = 0; loc_i < 4; loc_i++)
		ddEdxSigmaNumHitsVector_FDC_KPlus.push_back(3 + 3*loc_i);

	unsigned int locNumBetaBins = 240;
	unsigned int locNumdEdxBins = 5000;
	unsigned int locNumMomentumBins = 500;
	unsigned int locNumBetaGammaBins = 1500;
	unsigned int locNumHitsBins = 20;
	unsigned int locThetaBins = 600;
	unsigned int locNumdxBins = 500;
	float locBetaMin = 0.0, locBetaMax = 1.0;
	float locdEdxMin = 0.0, locdEdxMax = 100.0;
	float locMomentumMin = 0.0, locMomentumMax = 3.0;
	float locBetaGammaMin = 0.0, locBetaGammaMax = 15.0;
	float locNumHitsMin = -0.5, locNumHitsMax = 19.5;
	float locThetaMin = 0.0, locThetaMax = 150.0;
	float locdxMin = 0.0, locdxMax = 50.0;

	//no dx bins, no #hits cutoff
	//also do: dx bins (probably very small dependence), different #hits cutoff
	dSelectorHist_dEdxVsBeta_HitsCutoff_FDC = new TH2F("dSelectorHist_dEdxVsBeta_HitsCutoff_FDC", ";#beta;FDC #frac{dE}{dx} (keV/cm)", locNumBetaBins, locBetaMin, locBetaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsP_HitsCutoff_FDC = new TH2F("dSelectorHist_dEdxVsP_HitsCutoff_FDC", ";p (GeV/c);FDC #frac{dE}{dx} (keV/cm)", locNumMomentumBins, locMomentumMin, locMomentumMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC", ";#beta#gamma;FDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBeta_HitsCutoff_CDC = new TH2F("dSelectorHist_dEdxVsBeta_HitsCutoff_CDC", ";#beta;CDC #frac{dE}{dx} (keV/cm)", locNumBetaBins, locBetaMin, locBetaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsP_HitsCutoff_CDC = new TH2F("dSelectorHist_dEdxVsP_HitsCutoff_CDC", ";p (GeV/c);CDC #frac{dE}{dx} (keV/cm)", locNumMomentumBins, locMomentumMin, locMomentumMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC", ";#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);

	dSelectorHist_NumHitsVsTheta_CDC = new TH2F("dSelectorHist_NumHitsVsTheta_CDC", ";#theta#circ;#CDC Hits / 2", locThetaBins, locThetaMin, locThetaMax, locNumHitsBins, locNumHitsMin, locNumHitsMax);
	dSelectorHist_dxVsNumHits_CDC = new TH2F("dSelectorHist_dxVsNumHits_CDC", ";#CDC Hits / 2;dx (cm)", locNumHitsBins, locNumHitsMin, locNumHitsMax, locNumdxBins, locdxMin, locdxMax);
	dSelectorHist_dxVsTheta_CDC = new TH2F("dSelectorHist_dxVsTheta_CDC", ";#theta#circ;CDC dx (cm)", locThetaBins, locThetaMin, locThetaMax, locNumdxBins, locdxMin, locdxMax);
	dSelectorHist_NumHitsVsTheta_FDC = new TH2F("dSelectorHist_NumHitsVsTheta_FDC", ";#theta#circ;#FDC Hits / 2", locThetaBins, locThetaMin, locThetaMax, locNumHitsBins, locNumHitsMin, locNumHitsMax);
	dSelectorHist_dxVsNumHits_FDC = new TH2F("dSelectorHist_dxVsNumHits_FDC", ";#FDC Hits / 2;dx (cm)", locNumHitsBins, locNumHitsMin, locNumHitsMax, locNumdxBins, locdxMin, locdxMax);
	dSelectorHist_dxVsTheta_FDC = new TH2F("dSelectorHist_dxVsTheta_FDC", ";#theta#circ;FDC dx (cm)", locThetaBins, locThetaMin, locThetaMax, locNumdxBins, locdxMin, locdxMax);

	dSelectorHist_NumHitsFDCVsNumHitsCDC = new TH2F("dSelectorHist_NumHitsFDCVsNumHitsCDC", ";#FDC Hits / 2;#CDC Hits / 2", locNumHitsBins, locNumHitsMin, locNumHitsMax, locNumHitsBins, locNumHitsMin, locNumHitsMax);


	dSelectorHist_NotEnoughHitsInTheta = new TH1F("dSelectorHist_NotEnoughHitsInTheta", "Not enough dE/dx Info;#theta#circ", locThetaBins, locThetaMin, locThetaMax);
	dSelectorHist_ThetaDistribution = new TH1F("dSelectorHist_ThetaDistribution", "All Tracks;#theta#circ", locThetaBins, locThetaMin, locThetaMax);
	dSelectorHist_PercentageNotEnoughHitsInTheta = new TH1F("dSelectorHist_PercentageNotEnoughHitsInTheta", "% Not enough dE/dx Info;#theta#circ", locThetaBins, locThetaMin, locThetaMax);

	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_3Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_3Hits", "# FDC Hits Used = 3;#beta#gamma;FDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_6Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_6Hits", "# FDC Hits Used = 6;#beta#gamma;FDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_9Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_9Hits", "# FDC Hits Used = 9;#beta#gamma;FDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_12Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_12Hits", "# FDC Hits Used = 12;#beta#gamma;FDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);

	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_2Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_2Hits", "# CDC Hits Used = 2;#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_4Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_4Hits", "# CDC Hits Used = 4;#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_6Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_6Hits", "# CDC Hits Used = 6;#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_8Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_8Hits", "# CDC Hits Used = 8;#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_10Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_10Hits", "# CDC Hits Used = 10;#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_12Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_12Hits", "# CDC Hits Used = 12;#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);
	dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_14Hits = new TH2F("dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_14Hits", "# CDC Hits Used = 14;#beta#gamma;CDC #frac{dE}{dx} (keV/cm)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEdxBins, locdEdxMin, locdEdxMax);

	//confidence level
	unsigned int locNumConfidenceLevelBins = 400;
	dSelectorHist_ConfidenceLevel_FDC = new TH1F("dSelectorHist_ConfidenceLevel_FDC", ";Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC = new TH1F("dSelectorHist_ConfidenceLevel_CDC", ";Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_Both = new TH1F("dSelectorHist_ConfidenceLevel_Both", ";Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);

	dSelectorHist_ConfidenceLevel_FDC_3Hits = new TH1F("dSelectorHist_ConfidenceLevel_FDC_3Hits", "# FDC Hits Used = 3;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_6Hits = new TH1F("dSelectorHist_ConfidenceLevel_FDC_6Hits", "# FDC Hits Used = 6;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_9Hits = new TH1F("dSelectorHist_ConfidenceLevel_FDC_9Hits", "# FDC Hits Used = 9;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_12Hits = new TH1F("dSelectorHist_ConfidenceLevel_FDC_12Hits", "# FDC Hits Used = 12;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);

	dSelectorHist_ConfidenceLevel_CDC_4Hits = new TH1F("dSelectorHist_ConfidenceLevel_CDC_4Hits", "# CDC Hits Used = 4;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_6Hits = new TH1F("dSelectorHist_ConfidenceLevel_CDC_6Hits", "# CDC Hits Used = 6;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_8Hits = new TH1F("dSelectorHist_ConfidenceLevel_CDC_8Hits", "# CDC Hits Used = 8;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_10Hits = new TH1F("dSelectorHist_ConfidenceLevel_CDC_10Hits", "# CDC Hits Used = 10;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_12Hits = new TH1F("dSelectorHist_ConfidenceLevel_CDC_12Hits", "# CDC Hits Used = 12;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_14Hits = new TH1F("dSelectorHist_ConfidenceLevel_CDC_14Hits", "# CDC Hits Used = 14;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);

	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin1 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin1", "0.3 < #beta#gamma < 0.6;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin2 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin2", "0.6 < #beta#gamma < 0.9;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin3 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin3", "0.9 < #beta#gamma < 1.2;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin4 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin4", "1.2 < #beta#gamma < 1.5;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin5 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin5", "1.5 < #beta#gamma < 1.8;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin6 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin6", "1.8 < #beta#gamma < 2.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin7 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin7", "2.1 < #beta#gamma < 4.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin8 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin8", "4.1 < #beta#gamma < 6.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin9 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin9", "6.1 < #beta#gamma < 8.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin10 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin10", "8.1 < #beta#gamma < 10.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin11 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin11", "10.1 < #beta#gamma < 12.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin12 = new TH1F("dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin12", "12.1 < #beta#gamma < 14.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);

	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin1 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin1", "0.3 < #beta#gamma < 0.6;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin2 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin2", "0.6 < #beta#gamma < 0.9;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin3 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin3", "0.9 < #beta#gamma < 1.2;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin4 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin4", "1.2 < #beta#gamma < 1.5;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin5 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin5", "1.5 < #beta#gamma < 1.8;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin6 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin6", "1.8 < #beta#gamma < 2.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin7 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin7", "2.1 < #beta#gamma < 4.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin8 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin8", "4.1 < #beta#gamma < 6.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin9 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin9", "6.1 < #beta#gamma < 8.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin10 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin10", "8.1 < #beta#gamma < 10.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin11 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin11", "10.1 < #beta#gamma < 12.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
	dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin12 = new TH1F("dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin12", "12.1 < #beta#gamma < 14.1;Confidence Level", locNumConfidenceLevelBins, 0.0, 1.0);
}

void DCdEdxSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t DCdEdxSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either DCdEdxSelector::GetEntry() or TBranch::GetEntry()
   // to read either all or the required parts of the data. When processing
   // keyed objects with PROOF, the object is already loaded and is available
   // via the fObject pointer.
   //
   // This function should contain the "body" of the analysis. It can contain
   // simple or elaborate selection criteria, run algorithms on the data
   // of the event and typically fill histograms.
   //
   // The processing can be stopped by calling Abort().
   //
   // Use fStatus to set the return value of TTree::Process().
   //
   // The return value is currently not used.

	GetEntry(entry);
	if((entry % 1000) == 0)
		cout << "entry = " << entry << endl;

	double ddEdxkeV_FDC = 1000000.0*ddEdx_FDC; //convert GeV to keV (gas!)
	double ddEdxkeV_CDC = 1000000.0*ddEdx_CDC; //convert GeV to keV (gas!)

	float locBetaGamma = dBeta/sqrt(1.0 - dBeta*dBeta);
	float locTheta_Degrees = dTheta*180.0/3.141592654;

	unsigned int locMinimumNumHitsUsed = 3; //minimum num hits = 2x this

	dSelectorHist_NumHitsVsTheta_CDC->Fill(locTheta_Degrees, dNumHitsUsedFordEdx_CDC);
	dSelectorHist_dxVsNumHits_CDC->Fill(dNumHitsUsedFordEdx_CDC, ddx_CDC);
	dSelectorHist_dxVsTheta_CDC->Fill(locTheta_Degrees, ddx_CDC);
	dSelectorHist_NumHitsVsTheta_FDC->Fill(locTheta_Degrees, dNumHitsUsedFordEdx_FDC);
	dSelectorHist_dxVsNumHits_FDC->Fill(dNumHitsUsedFordEdx_FDC, ddx_FDC);
	dSelectorHist_dxVsTheta_FDC->Fill(locTheta_Degrees, ddx_FDC);

	dSelectorHist_NumHitsFDCVsNumHitsCDC->Fill(dNumHitsUsedFordEdx_CDC, dNumHitsUsedFordEdx_FDC);

	if((dNumHitsUsedFordEdx_FDC < locMinimumNumHitsUsed) && (dNumHitsUsedFordEdx_CDC < locMinimumNumHitsUsed)) //not enough hits (6+) in either!
		dSelectorHist_NotEnoughHitsInTheta->Fill(locTheta_Degrees);
	dSelectorHist_ThetaDistribution->Fill(locTheta_Degrees);

	if(dNumHitsUsedFordEdx_FDC >= locMinimumNumHitsUsed){ //6+ original hits
		dSelectorHist_dEdxVsBeta_HitsCutoff_FDC->Fill(dBeta, ddEdxkeV_FDC);
		dSelectorHist_dEdxVsP_HitsCutoff_FDC->Fill(dMomentum, ddEdxkeV_FDC);
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC->Fill(locBetaGamma, ddEdxkeV_FDC);
	}

	if(dNumHitsUsedFordEdx_CDC >= locMinimumNumHitsUsed){ //6+ original hits
		dSelectorHist_dEdxVsBeta_HitsCutoff_CDC->Fill(dBeta, ddEdxkeV_CDC);
		dSelectorHist_dEdxVsP_HitsCutoff_CDC->Fill(dMomentum, ddEdxkeV_CDC);
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC->Fill(locBetaGamma, ddEdxkeV_CDC);
	}

	if(dNumHitsUsedFordEdx_FDC == 3)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_3Hits->Fill(locBetaGamma, ddEdxkeV_FDC);
	if(dNumHitsUsedFordEdx_FDC == 6)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_6Hits->Fill(locBetaGamma, ddEdxkeV_FDC);
	if(dNumHitsUsedFordEdx_FDC == 9)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_9Hits->Fill(locBetaGamma, ddEdxkeV_FDC);
	if(dNumHitsUsedFordEdx_FDC == 12)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_FDC_12Hits->Fill(locBetaGamma, ddEdxkeV_FDC);

	if(dNumHitsUsedFordEdx_CDC == 2)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_2Hits->Fill(locBetaGamma, ddEdxkeV_CDC);
	if(dNumHitsUsedFordEdx_CDC == 4)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_4Hits->Fill(locBetaGamma, ddEdxkeV_CDC);
	if(dNumHitsUsedFordEdx_CDC == 6)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_6Hits->Fill(locBetaGamma, ddEdxkeV_CDC);
	if(dNumHitsUsedFordEdx_CDC == 8)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_8Hits->Fill(locBetaGamma, ddEdxkeV_CDC);
	if(dNumHitsUsedFordEdx_CDC == 10)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_10Hits->Fill(locBetaGamma, ddEdxkeV_CDC);
	if(dNumHitsUsedFordEdx_CDC == 12)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_12Hits->Fill(locBetaGamma, ddEdxkeV_CDC);
	if(dNumHitsUsedFordEdx_CDC == 14)
		dSelectorHist_dEdxVsBetaGamma_HitsCutoff_CDC_14Hits->Fill(locBetaGamma, ddEdxkeV_CDC);

	double locFOM;
	bool locCalcStatus = Calc_FOM(locFOM);
/*
	if(fabs(dFOM - locFOM) > 0.001){
		cout << "orig fom, new fom, diff = " << dFOM << ", " << locFOM << ", " << (dFOM - locFOM) << endl;
		cout << "#hits CDC, #hits FDC, dE/dx CDC, dE/dx FDC = " << dNumHitsUsedFordEdx_CDC << ", " << dNumHitsUsedFordEdx_FDC << ", " << ddEdx_CDC << ", " << ddEdx_FDC << endl;
	}
*/
	if(dCalcFOMManuallyFlag == true){
		if(locCalcStatus == false)
			return kTRUE;
		dFOM = locFOM;
	}


	if((dNumHitsUsedFordEdx_FDC >= 3) && (dNumHitsUsedFordEdx_CDC < 3)){ //FDC only
		dSelectorHist_ConfidenceLevel_FDC->Fill(dFOM);
		if(dNumHitsUsedFordEdx_FDC == 3)
			dSelectorHist_ConfidenceLevel_FDC_3Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_FDC == 6)
			dSelectorHist_ConfidenceLevel_FDC_6Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_FDC == 9)
			dSelectorHist_ConfidenceLevel_FDC_9Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_FDC == 12)
			dSelectorHist_ConfidenceLevel_FDC_12Hits->Fill(dFOM);

		if((locBetaGamma >= 0.3) && (locBetaGamma < 0.6))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin1->Fill(dFOM);
		if((locBetaGamma >= 0.6) && (locBetaGamma < 0.9))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin2->Fill(dFOM);
		if((locBetaGamma >= 0.9) && (locBetaGamma < 1.2))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin3->Fill(dFOM);
		if((locBetaGamma >= 1.2) && (locBetaGamma < 1.5))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin4->Fill(dFOM);
		if((locBetaGamma >= 1.5) && (locBetaGamma < 1.8))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin5->Fill(dFOM);
		if((locBetaGamma >= 1.8) && (locBetaGamma < 2.1))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin6->Fill(dFOM);
		if((locBetaGamma >= 2.1) && (locBetaGamma < 4.1))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin7->Fill(dFOM);
		if((locBetaGamma >= 4.1) && (locBetaGamma < 6.1))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin8->Fill(dFOM);
		if((locBetaGamma >= 6.1) && (locBetaGamma < 8.1))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin9->Fill(dFOM);
		if((locBetaGamma >= 8.1) && (locBetaGamma < 10.1))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin10->Fill(dFOM);
		if((locBetaGamma >= 10.1) && (locBetaGamma < 12.1))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin11->Fill(dFOM);
		if((locBetaGamma >= 12.1) && (locBetaGamma < 14.1))
			dSelectorHist_ConfidenceLevel_FDC_BetaGammaBin12->Fill(dFOM);
	}
	if((dNumHitsUsedFordEdx_CDC >= 3) && (dNumHitsUsedFordEdx_FDC < 3)){ //CDC only
		dSelectorHist_ConfidenceLevel_CDC->Fill(dFOM);
		if(dNumHitsUsedFordEdx_CDC == 4)
			dSelectorHist_ConfidenceLevel_CDC_4Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_CDC == 6)
			dSelectorHist_ConfidenceLevel_CDC_6Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_CDC == 8)
			dSelectorHist_ConfidenceLevel_CDC_8Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_CDC == 10)
			dSelectorHist_ConfidenceLevel_CDC_10Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_CDC == 12)
			dSelectorHist_ConfidenceLevel_CDC_12Hits->Fill(dFOM);
		if(dNumHitsUsedFordEdx_CDC == 14)
			dSelectorHist_ConfidenceLevel_CDC_14Hits->Fill(dFOM);

		if((locBetaGamma >= 0.3) && (locBetaGamma < 0.6))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin1->Fill(dFOM);
		if((locBetaGamma >= 0.6) && (locBetaGamma < 0.9))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin2->Fill(dFOM);
		if((locBetaGamma >= 0.9) && (locBetaGamma < 1.2))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin3->Fill(dFOM);
		if((locBetaGamma >= 1.2) && (locBetaGamma < 1.5))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin4->Fill(dFOM);
		if((locBetaGamma >= 1.5) && (locBetaGamma < 1.8))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin5->Fill(dFOM);
		if((locBetaGamma >= 1.8) && (locBetaGamma < 2.1))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin6->Fill(dFOM);
		if((locBetaGamma >= 2.1) && (locBetaGamma < 4.1))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin7->Fill(dFOM);
		if((locBetaGamma >= 4.1) && (locBetaGamma < 6.1))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin8->Fill(dFOM);
		if((locBetaGamma >= 6.1) && (locBetaGamma < 8.1))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin9->Fill(dFOM);
		if((locBetaGamma >= 8.1) && (locBetaGamma < 10.1))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin10->Fill(dFOM);
		if((locBetaGamma >= 10.1) && (locBetaGamma < 12.1))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin11->Fill(dFOM);
		if((locBetaGamma >= 12.1) && (locBetaGamma < 14.1))
			dSelectorHist_ConfidenceLevel_CDC_BetaGammaBin12->Fill(dFOM);


	}

	if((dNumHitsUsedFordEdx_CDC >= 3) && (dNumHitsUsedFordEdx_FDC >= 3)){ //Both
		dSelectorHist_ConfidenceLevel_Both->Fill(dFOM);

		double locMeandEdx_CDC, locMeandEdx_FDC;
		if(!GetdEdxMean_CDC(dBeta, dNumHitsUsedFordEdx_CDC, locMeandEdx_CDC))
			return kTRUE;
		if(!GetdEdxMean_FDC(dBeta, dNumHitsUsedFordEdx_FDC, locMeandEdx_FDC))
			return kTRUE;

		//Correlation coefficient:
		//http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Geometric_interpretation


	}

   return kTRUE;
}

void DCdEdxSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void DCdEdxSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	float locNumNotEnoughHits, locNumGenerated, locNotEnoughHitsRatio, locNotEnoughHitsRatioUncertainty;
	for(unsigned int loc_i = 1; loc_i <= dSelectorHist_ThetaDistribution->GetNbinsX(); loc_i++){
		locNumNotEnoughHits = float(dSelectorHist_NotEnoughHitsInTheta->GetBinContent(loc_i));
		locNumGenerated = float(dSelectorHist_ThetaDistribution->GetBinContent(loc_i));
		locNotEnoughHitsRatio = locNumNotEnoughHits/locNumGenerated;
		locNotEnoughHitsRatioUncertainty = sqrt(locNumNotEnoughHits*(1.0 - locNotEnoughHitsRatio))/locNumGenerated;
		dSelectorHist_PercentageNotEnoughHitsInTheta->SetBinContent(loc_i, locNotEnoughHitsRatio);
		dSelectorHist_PercentageNotEnoughHitsInTheta->SetBinError(loc_i, locNotEnoughHitsRatioUncertainty);
	}

	dOutputFile->Write();
}

bool DCdEdxSelector::GetdEdxMean_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx){
	double locBetaGammaValue = locBeta/sqrt(1.0 - locBeta*locBeta);
	if(dParticleName == "Proton"){
		if(locBetaGammaValue < 0.3)
			return false;
		if(locBetaGammaValue > 2.15)
			return false;

		locMeandEdx = (ddEdxMeanFunc_CDC_Proton->Eval(locBetaGammaValue))/1000000.0;
		return true;
	}
	if(dParticleName == "KPlus"){
		if(locBetaGammaValue < 0.65)
			return false; //K+ very likely to decay, don't compute chisq!!
		if(locBetaGammaValue > 4.1)
			return false;

		locMeandEdx = (ddEdxMeanFunc_CDC_KPlus->Eval(locBetaGammaValue))/1000000.0;
		return true;
	}

	double locBetaGamma[44];
	double locdEdxMean[44];
	locBetaGamma[0] = 0.48;  locBetaGamma[1] = 0.8;  locBetaGamma[2] = 1.12;  locBetaGamma[3] = 1.44;  locBetaGamma[4] = 1.76;  locBetaGamma[5] = 2.08;  locBetaGamma[6] = 2.4;  locBetaGamma[7] = 2.72;  locBetaGamma[8] = 3.04;  locBetaGamma[9] = 3.36;  locBetaGamma[10] = 3.68;  locBetaGamma[11] = 4;  locBetaGamma[12] = 4.32;  locBetaGamma[13] = 4.64;  locBetaGamma[14] = 4.96;  locBetaGamma[15] = 5.28;  locBetaGamma[16] = 5.6;  locBetaGamma[17] = 5.92;  locBetaGamma[18] = 6.24;  locBetaGamma[19] = 6.56;  locBetaGamma[20] = 6.88;  locBetaGamma[21] = 7.2;  locBetaGamma[22] = 7.52;  locBetaGamma[23] = 7.84;  locBetaGamma[24] = 8.16;  locBetaGamma[25] = 8.48;  locBetaGamma[26] = 8.8;  locBetaGamma[27] = 9.12;  locBetaGamma[28] = 9.44;  locBetaGamma[29] = 9.76;  locBetaGamma[30] = 10.08;  locBetaGamma[31] = 10.4;  locBetaGamma[32] = 10.72;  locBetaGamma[33] = 11.04;  locBetaGamma[34] = 11.36;  locBetaGamma[35] = 11.68;  locBetaGamma[36] = 12;  locBetaGamma[37] = 12.32;  locBetaGamma[38] = 12.64;  locBetaGamma[39] = 12.96;  locBetaGamma[40] = 13.28;  locBetaGamma[41] = 13.6;  locBetaGamma[42] = 13.92;  locBetaGamma[43] = 14.24;
	locdEdxMean[0] = 2.40705;  locdEdxMean[1] = 2.64016;  locdEdxMean[2] = 1.97719;  locdEdxMean[3] = 1.74434;  locdEdxMean[4] = 1.60673;  locdEdxMean[5] = 1.52877;  locdEdxMean[6] = 1.47792;  locdEdxMean[7] = 1.44914;  locdEdxMean[8] = 1.44456;  locdEdxMean[9] = 1.44089;  locdEdxMean[10] = 1.44077;  locdEdxMean[11] = 1.44643;  locdEdxMean[12] = 1.46221;  locdEdxMean[13] = 1.47146;  locdEdxMean[14] = 1.48451;  locdEdxMean[15] = 1.49671;  locdEdxMean[16] = 1.50718;  locdEdxMean[17] = 1.51892;  locdEdxMean[18] = 1.53035;  locdEdxMean[19] = 1.54224;  locdEdxMean[20] = 1.5531;  locdEdxMean[21] = 1.56325;  locdEdxMean[22] = 1.57347;  locdEdxMean[23] = 1.58422;  locdEdxMean[24] = 1.59239;  locdEdxMean[25] = 1.60057;  locdEdxMean[26] = 1.61224;  locdEdxMean[27] = 1.62025;  locdEdxMean[28] = 1.6288;  locdEdxMean[29] = 1.63748;  locdEdxMean[30] = 1.64582;  locdEdxMean[31] = 1.65254;  locdEdxMean[32] = 1.65864;  locdEdxMean[33] = 1.66693;  locdEdxMean[34] = 1.67503;  locdEdxMean[35] = 1.68328;  locdEdxMean[36] = 1.68964;  locdEdxMean[37] = 1.69732;  locdEdxMean[38] = 1.70258;  locdEdxMean[39] = 1.70849;  locdEdxMean[40] = 1.71469;  locdEdxMean[41] = 1.72359;  locdEdxMean[42] = 1.73032;  locdEdxMean[43] = 1.73672;

	if(locBetaGammaValue < 0.81)
		return false;
	if(locBetaGammaValue > 14.2)
		return false;

	double locSlope, locIntercept;
	for(unsigned int loc_i = 0; loc_i < 43; loc_i++){
		if(locBetaGammaValue > locBetaGamma[loc_i + 1])
			continue;
		locSlope = (locdEdxMean[loc_i + 1] - locdEdxMean[loc_i])/(locBetaGamma[loc_i + 1] - locBetaGamma[loc_i]);
		locIntercept = locdEdxMean[loc_i] - locSlope*locBetaGamma[loc_i]; //y = mx + b, b = y - mx
		locMeandEdx = (locSlope*locBetaGammaValue + locIntercept)/1000000.0;
		return true;
	}

	return false;
}

bool DCdEdxSelector::GetdEdxMean_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locMeandEdx){
	double locBetaGammaValue = locBeta/sqrt(1.0 - locBeta*locBeta);

	if(dParticleName == "Proton"){
		if(locBetaGammaValue < 0.3)
			return false;
		if(locBetaGammaValue > 2.15)
			return false;

		locMeandEdx = (ddEdxMeanFunc_FDC_Proton->Eval(locBetaGammaValue))/1000000.0;
		return true;
	}
	if(dParticleName == "KPlus"){
		if(locBetaGammaValue < 0.9)
			return false; //K+ very likely to decay, don't compute chisq!!
		if(locBetaGammaValue > 4.1)
			return false;

		locMeandEdx = (ddEdxMeanFunc_FDC_KPlus->Eval(locBetaGammaValue))/1000000.0;
		return true;
	}

	double locBetaGamma[44];
	double locdEdxMean[44];
	locBetaGamma[0] = 0.48;  locBetaGamma[1] = 0.8;  locBetaGamma[2] = 1.12;  locBetaGamma[3] = 1.44;  locBetaGamma[4] = 1.76;  locBetaGamma[5] = 2.08;  locBetaGamma[6] = 2.4;  locBetaGamma[7] = 2.72;  locBetaGamma[8] = 3.04;  locBetaGamma[9] = 3.36;  locBetaGamma[10] = 3.68;  locBetaGamma[11] = 4;  locBetaGamma[12] = 4.32;  locBetaGamma[13] = 4.64;  locBetaGamma[14] = 4.96;  locBetaGamma[15] = 5.28;  locBetaGamma[16] = 5.6;  locBetaGamma[17] = 5.92;  locBetaGamma[18] = 6.24;  locBetaGamma[19] = 6.56;  locBetaGamma[20] = 6.88;  locBetaGamma[21] = 7.2;  locBetaGamma[22] = 7.52;  locBetaGamma[23] = 7.84;  locBetaGamma[24] = 8.16;  locBetaGamma[25] = 8.48;  locBetaGamma[26] = 8.8;  locBetaGamma[27] = 9.12;  locBetaGamma[28] = 9.44;  locBetaGamma[29] = 9.76;  locBetaGamma[30] = 10.08;  locBetaGamma[31] = 10.4;  locBetaGamma[32] = 10.72;  locBetaGamma[33] = 11.04;  locBetaGamma[34] = 11.36;  locBetaGamma[35] = 11.68;  locBetaGamma[36] = 12;  locBetaGamma[37] = 12.32;  locBetaGamma[38] = 12.64;  locBetaGamma[39] = 12.96;  locBetaGamma[40] = 13.28;  locBetaGamma[41] = 13.6;  locBetaGamma[42] = 13.92;  locBetaGamma[43] = 14.24;
	locdEdxMean[0] = 1.32105;  locdEdxMean[1] = 1.50904;  locdEdxMean[2] = 1.50873;  locdEdxMean[3] = 1.57426;  locdEdxMean[4] = 1.52213;  locdEdxMean[5] = 1.45914;  locdEdxMean[6] = 1.4271;  locdEdxMean[7] = 1.42623;  locdEdxMean[8] = 1.42268;  locdEdxMean[9] = 1.42632;  locdEdxMean[10] = 1.43233;  locdEdxMean[11] = 1.43742;  locdEdxMean[12] = 1.4467;  locdEdxMean[13] = 1.45642;  locdEdxMean[14] = 1.46337;  locdEdxMean[15] = 1.4719;  locdEdxMean[16] = 1.48252;  locdEdxMean[17] = 1.49443;  locdEdxMean[18] = 1.50234;  locdEdxMean[19] = 1.513;  locdEdxMean[20] = 1.52394;  locdEdxMean[21] = 1.5289;  locdEdxMean[22] = 1.53878;  locdEdxMean[23] = 1.54396;  locdEdxMean[24] = 1.55427;  locdEdxMean[25] = 1.56105;  locdEdxMean[26] = 1.5674;  locdEdxMean[27] = 1.57221;  locdEdxMean[28] = 1.58027;  locdEdxMean[29] = 1.5865;  locdEdxMean[30] = 1.59519;  locdEdxMean[31] = 1.60211;  locdEdxMean[32] = 1.60751;  locdEdxMean[33] = 1.61522;  locdEdxMean[34] = 1.62066;  locdEdxMean[35] = 1.6277;  locdEdxMean[36] = 1.63227;  locdEdxMean[37] = 1.63661;  locdEdxMean[38] = 1.64457;  locdEdxMean[39] = 1.65043;  locdEdxMean[40] = 1.6557;  locdEdxMean[41] = 1.66085;  locdEdxMean[42] = 1.66486;  locdEdxMean[43] = 1.66916;

	if(locBetaGammaValue < 1.77)
		return false;
	if(locBetaGammaValue > 14.2)
		return false;

	double locSlope, locIntercept;
	for(unsigned int loc_i = 0; loc_i < 43; loc_i++){
		if(locBetaGammaValue > locBetaGamma[loc_i + 1])
			continue;
		locSlope = (locdEdxMean[loc_i + 1] - locdEdxMean[loc_i])/(locBetaGamma[loc_i + 1] - locBetaGamma[loc_i]);
		locIntercept = locdEdxMean[loc_i] - locSlope*locBetaGamma[loc_i]; //y = mx + b, b = y - mx
		locMeandEdx = (locSlope*locBetaGammaValue + locIntercept)/1000000.0;
		return true;
	}
	return false;
}


bool DCdEdxSelector::GetdEdxSigma_FDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx){
	double locBetaGamma = locBeta/sqrt(1.0 - locBeta*locBeta);
	double locSigmadEdx_LowSide, locSigmadEdx_HighSide; //for linear interpolation/extrapolation

	vector<int> locNumHitsVector;
	TObjArray *locFuncArray;
	if(dParticleName == "Proton"){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_Proton;
		locFuncArray = dSigmaFuncArray_FDC_Proton;
	}
	if(dParticleName == "KPlus"){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_KPlus;
		locFuncArray = dSigmaFuncArray_FDC_KPlus;
	}
	if(dParticleName == "PiPlus"){
		locNumHitsVector = ddEdxSigmaNumHitsVector_FDC_PiPlus;
		locFuncArray = dSigmaFuncArray_FDC_PiPlus;
	}

	double locSlope, locIntercept;

	for(unsigned int loc_i = 0; loc_i < locNumHitsVector.size(); loc_i++){
		if(locNumHitsUsedFordEdx == locNumHitsVector[loc_i]){
			locSigmadEdx = (((TF1*)(*locFuncArray)[loc_i])->Eval(locBetaGamma))/1000000.0;
			return true;
		}
		if(loc_i == (locNumHitsVector.size() - 1)){
//extrapolate!!
			locSigmadEdx_LowSide = ((TF1*)(*locFuncArray)[loc_i - 1])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = ((TF1*)(*locFuncArray)[loc_i])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i] - locNumHitsVector[loc_i - 1]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return true;
		}
		if((locNumHitsUsedFordEdx > locNumHitsVector[loc_i]) && (locNumHitsUsedFordEdx < locNumHitsVector[loc_i + 1])){
			locSigmadEdx_LowSide = ((TF1*)(*locFuncArray)[loc_i])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = ((TF1*)(*locFuncArray)[loc_i + 1])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i + 1] - locNumHitsVector[loc_i]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i + 1]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return true;
		}
	}
	return true;
}

bool DCdEdxSelector::GetdEdxSigma_CDC(double locBeta, unsigned int locNumHitsUsedFordEdx, double& locSigmadEdx){
	double locBetaGamma = locBeta/sqrt(1.0 - locBeta*locBeta);
	double locSigmadEdx_LowSide, locSigmadEdx_HighSide; //for linear interpolation/extrapolation

	vector<int> locNumHitsVector;
	TObjArray *locFuncArray;
	if(dParticleName == "Proton"){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_Proton;
		locFuncArray = dSigmaFuncArray_CDC_Proton;
	}
	if(dParticleName == "KPlus"){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_KPlus;
		locFuncArray = dSigmaFuncArray_CDC_KPlus;
	}
	if(dParticleName == "PiPlus"){
		locNumHitsVector = ddEdxSigmaNumHitsVector_CDC_PiPlus;
		locFuncArray = dSigmaFuncArray_CDC_PiPlus;
	}

	double locSlope, locIntercept;

	for(unsigned int loc_i = 0; loc_i < locNumHitsVector.size(); loc_i++){
		if(locNumHitsUsedFordEdx == locNumHitsVector[loc_i]){
			locSigmadEdx = (((TF1*)(*locFuncArray)[loc_i])->Eval(locBetaGamma))/1000000.0;
			return true;
		}
		if(loc_i == (locNumHitsVector.size() - 1)){
//extrapolate!!
			locSigmadEdx_LowSide = ((TF1*)(*locFuncArray)[loc_i - 1])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = ((TF1*)(*locFuncArray)[loc_i])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i] - locNumHitsVector[loc_i - 1]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return true;
		}
		if((locNumHitsUsedFordEdx > locNumHitsVector[loc_i]) && (locNumHitsUsedFordEdx < locNumHitsVector[loc_i + 1])){
			locSigmadEdx_LowSide = ((TF1*)(*locFuncArray)[loc_i])->Eval(locBetaGamma);
			locSigmadEdx_HighSide = ((TF1*)(*locFuncArray)[loc_i + 1])->Eval(locBetaGamma);

			locSlope = (locSigmadEdx_HighSide - locSigmadEdx_LowSide)/(locNumHitsVector[loc_i + 1] - locNumHitsVector[loc_i]);
			locIntercept = locSigmadEdx_HighSide - locSlope*locNumHitsVector[loc_i + 1]; //y = mx + b, b = y - mx
			locSigmadEdx = (locSlope*locNumHitsUsedFordEdx + locIntercept)/1000000.0;

			return true;
		}
	}

	return true;
}


bool DCdEdxSelector::Calc_FOM(double& locFOM){
	double locDCdEdx, locChiSq;
	unsigned int locNDF;
	unsigned int locMinimumNumberUsedHitsForConfidence = 3; //dE/dx is landau-distributed, so to approximate Gaussian must remove hits with largest dE/dx //3 means 6 or more hits originally

	bool locUseCDCHitsFlag = (dNumHitsUsedFordEdx_CDC > locMinimumNumberUsedHitsForConfidence) ? true : false;
	bool locUseFDCHitsFlag = (dNumHitsUsedFordEdx_FDC > locMinimumNumberUsedHitsForConfidence) ? true : false;

//cout << "#hits cdc, fdc = " << locNumHitsUsedFordEdx_CDC << ", " << locNumHitsUsedFordEdx_FDC << endl;
	if((locUseCDCHitsFlag == false) && (locUseFDCHitsFlag == false))
		return false; //not enough hits, use other sources of information for PID

	double locMeandEdx_FDC, locMeandEdx_CDC, locSigmadEdx_FDC, locSigmadEdx_CDC;
	double locDeltadEdx_CDC = 0.0, locDeltadEdx_FDC = 0.0;

	if(locUseCDCHitsFlag == true){
		if(GetdEdxMean_CDC(dBeta, dNumHitsUsedFordEdx_CDC, locMeandEdx_CDC) != true)
			return false;
		if(GetdEdxSigma_CDC(dBeta, dNumHitsUsedFordEdx_CDC, locSigmadEdx_CDC) != true)
			return false;
		locDeltadEdx_CDC = ddEdx_CDC - locMeandEdx_CDC;
//cout << "CDC: #hits, beta, actual, mean, sigma = " << dNumHitsUsedFordEdx_CDC << ", " << dBeta << ", " << ddEdx_CDC << ", " << locMeandEdx_CDC << ", " << locSigmadEdx_CDC << endl;
	}

	if(locUseFDCHitsFlag == true){
		if(GetdEdxMean_FDC(dBeta, dNumHitsUsedFordEdx_FDC, locMeandEdx_FDC) != true)
			return false;
		if(GetdEdxSigma_FDC(dBeta, dNumHitsUsedFordEdx_FDC, locSigmadEdx_FDC) != true)
			return false;
		locDeltadEdx_FDC = ddEdx_FDC - locMeandEdx_FDC;
//cout << "FDC: #hits, beta, actual, mean, sigma = " << dNumHitsUsedFordEdx_FDC << ", " << dBeta << ", " << ddEdx_FDC << ", " << locMeandEdx_FDC << ", " << locSigmadEdx_FDC << endl;
	}

	if((locUseCDCHitsFlag == true) && ((locUseFDCHitsFlag == false) || (dNumHitsUsedFordEdx_CDC >= dNumHitsUsedFordEdx_FDC))){
		locDCdEdx = ddEdx_CDC;
		locChiSq = locDeltadEdx_CDC/locSigmadEdx_CDC;
		locChiSq *= locChiSq;
		locNDF = 1;
		locFOM = TMath::Prob(locChiSq, locNDF);
		return true;
	}

	if((locUseFDCHitsFlag == true) && ((locUseCDCHitsFlag == false) || (dNumHitsUsedFordEdx_FDC > dNumHitsUsedFordEdx_CDC)) ){
		locDCdEdx = ddEdx_FDC;
		locChiSq = locDeltadEdx_FDC/locSigmadEdx_FDC;
		locChiSq *= locChiSq;
		locNDF = 1;
		locFOM = TMath::Prob(locChiSq, locNDF);
		return true;
	}
	return true;
}

