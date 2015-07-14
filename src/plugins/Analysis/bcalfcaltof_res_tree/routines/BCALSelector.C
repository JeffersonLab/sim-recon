#define BCALSelector_cxx
// The class definition in BCALSelector.h has been generated automatically
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
// Root > T->Process("BCALSelector.C")
// Root > T->Process("BCALSelector.C","some options")
// Root > T->Process("BCALSelector.C+")
//

#include "BCALSelector.h"
#include <TH2.h>
#include <TStyle.h>


void BCALSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

	dOutputFile = new TFile("dh_BCALMCComparisonHists.root", "RECREATE");

//	float locRMin = 64.0, locRMax = 67.0; //pre depth correction
	float locRMin = 65.0, locRMax = 75.0;
	float locPhiMin = -3.2, locPhiMax = 3.2;
	float locZMin = 0.0, locZMax = 430.0;
	float locEMin = 0.0, locEMax = 9.1;
	float locTMin = 2.0, locTMax = 15.0;

//	float locDeltaRMin = -2.0, locDeltaRMax = 12.0; //pre depth correction
	float locDeltaRMin = -4.0, locDeltaRMax = 4.0;
	float locDeltaPhiMin = -0.05, locDeltaPhiMax = 0.05;
	float locDeltaZMin = -12.0, locDeltaZMax = 12.0;
	float locDeltaEMin = -1.0, locDeltaEMax = 1.0;
	float locDeltaTMin = -2.0, locDeltaTMax = 2.0;

	dPluginHist_BCAL_PathLengthCorrection = new TH2F("dPluginHist_BCAL_PathLengthCorrection", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);Path Length Correction (cm)", 100, locEMin, locEMax, 100, -1.0, 25.0);
	dPluginHist_BCAL_PathLengthCorrectionPostEVsZ = new TH2F("dPluginHist_BCAL_PathLengthCorrectionPostEVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);Path Length Correction (cm)", 400, locZMin, locZMax, 400, -30.0, 10.0);

	//DeltaR Dependence
	dPluginHist_BCAL_DeltaRVsR = new TH2F("dPluginHist_BCAL_DeltaRVsR", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured R (cm);BCAL Shower #DeltaR (cm)", 200, locRMin, locRMax, 200, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaRVsPhi = new TH2F("dPluginHist_BCAL_DeltaRVsPhi", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured #phi (rad);BCAL Shower #DeltaR (cm)", 100, locPhiMin, locPhiMax, 100, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaRVsZ = new TH2F("dPluginHist_BCAL_DeltaRVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #DeltaR (cm)", 100, locZMin, locZMax, 100, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaRVsE = new TH2F("dPluginHist_BCAL_DeltaRVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #DeltaR (cm)", 100, locEMin, locEMax, 100, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaRVsT = new TH2F("dPluginHist_BCAL_DeltaRVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #DeltaR (cm)", 100, locTMin, locTMax, 100, locDeltaRMin, locDeltaRMax);

	//DeltaPhi Dependence
	dPluginHist_BCAL_DeltaPhiVsR = new TH2F("dPluginHist_BCAL_DeltaPhiVsR", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured R (cm);BCAL Shower #Delta#Phi (rad)", 100, locRMin, locRMax, 100, locDeltaPhiMin, locDeltaPhiMax);
	dPluginHist_BCAL_DeltaPhiVsPhi = new TH2F("dPluginHist_BCAL_DeltaPhiVsPhi", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured #phi (rad);BCAL Shower #Delta#Phi (rad)", 100, locPhiMin, locPhiMax, 400, locDeltaPhiMin, locDeltaPhiMax);
	dPluginHist_BCAL_DeltaPhiVsZ = new TH2F("dPluginHist_BCAL_DeltaPhiVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #Delta#Phi (rad)", 100, locZMin, locZMax, 100, locDeltaPhiMin, locDeltaPhiMax);
	dPluginHist_BCAL_DeltaPhiVsE = new TH2F("dPluginHist_BCAL_DeltaPhiVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #Delta#Phi (rad)", 100, locEMin, locEMax, 100, locDeltaPhiMin, locDeltaPhiMax);
	dPluginHist_BCAL_DeltaPhiVsT = new TH2F("dPluginHist_BCAL_DeltaPhiVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #Delta#Phi (rad)", 100, locTMin, locTMax, 100, locDeltaPhiMin, locDeltaPhiMax);

	//DeltaZ Dependence
	dPluginHist_BCAL_DeltaZVsR = new TH2F("dPluginHist_BCAL_DeltaZVsR", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured R (cm);BCAL Shower #DeltaZ (cm)", 100, locRMin, locRMax, 100, locDeltaZMin, locDeltaZMax);
	dPluginHist_BCAL_DeltaZVsPhi = new TH2F("dPluginHist_BCAL_DeltaZVsPhi", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured #phi (rad);BCAL Shower #DeltaZ (cm)", 100, locPhiMin, locPhiMax, 100, locDeltaZMin, locDeltaZMax);
	dPluginHist_BCAL_DeltaZVsZ = new TH2F("dPluginHist_BCAL_DeltaZVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #DeltaZ (cm)", 400, locZMin, locZMax, 400, locDeltaZMin, locDeltaZMax);
	dPluginHist_BCAL_DeltaZVsE = new TH2F("dPluginHist_BCAL_DeltaZVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #DeltaZ (cm)", 100, locEMin, locEMax, 100, locDeltaZMin, locDeltaZMax);
	dPluginHist_BCAL_DeltaZVsT = new TH2F("dPluginHist_BCAL_DeltaZVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #DeltaZ (cm)", 100, locTMin, locTMax, 100, locDeltaZMin, locDeltaZMax);

	//DeltaE Dependence
	dPluginHist_BCAL_DeltaEVsR = new TH2F("dPluginHist_BCAL_DeltaEVsR", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured R (cm);BCAL Shower #DeltaE (GeV)", 100, locRMin, locRMax, 100, locDeltaEMin, locDeltaEMax);
	dPluginHist_BCAL_DeltaEVsPhi = new TH2F("dPluginHist_BCAL_DeltaEVsPhi", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured #phi (rad);BCAL Shower #DeltaE (GeV)", 100, locPhiMin, locPhiMax, 100, locDeltaEMin, locDeltaEMax);
	dPluginHist_BCAL_DeltaEVsZ = new TH2F("dPluginHist_BCAL_DeltaEVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #DeltaE (GeV)", 100, locZMin, locZMax, 100, locDeltaEMin, locDeltaEMax);
	dPluginHist_BCAL_DeltaEVsE = new TH2F("dPluginHist_BCAL_DeltaEVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #DeltaE (GeV)", 400, locEMin, locEMax, 400, locDeltaEMin, locDeltaEMax);
	dPluginHist_BCAL_DeltaEVsT = new TH2F("dPluginHist_BCAL_DeltaEVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #DeltaE (GeV)", 100, locTMin, locTMax, 100, locDeltaEMin, locDeltaEMax);

	//DeltaT Dependence
	dPluginHist_BCAL_DeltaTVsR = new TH2F("dPluginHist_BCAL_DeltaTVsR", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured R (cm);BCAL Shower #DeltaT (ns)", 100, locRMin, locRMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_BCAL_DeltaTVsPhi = new TH2F("dPluginHist_BCAL_DeltaTVsPhi", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured #phi (rad);BCAL Shower #DeltaT (ns)", 100, locPhiMin, locPhiMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_BCAL_DeltaTVsZ = new TH2F("dPluginHist_BCAL_DeltaTVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #DeltaT (ns)", 100, locZMin, locZMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_BCAL_DeltaTVsE = new TH2F("dPluginHist_BCAL_DeltaTVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #DeltaT (ns)", 100, locEMin, locEMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_BCAL_DeltaTVsT = new TH2F("dPluginHist_BCAL_DeltaTVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #DeltaT (ns)", 100, locTMin, locTMax, 400, locDeltaTMin, locDeltaTMax);

	//Common Dependence
	dPluginHist_BCAL_DeltaRVsDeltaPhi = new TH2F("dPluginHist_BCAL_DeltaRVsDeltaPhi", "Simulation Reconstructed Uncertainty Study;BCAL Shower #Delta#phi (rad);BCAL Shower #DeltaR (cm)", 100, locDeltaPhiMin, locDeltaPhiMax, 100, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaRVsDeltaZ = new TH2F("dPluginHist_BCAL_DeltaRVsDeltaZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaZ (cm);BCAL Shower #DeltaR (cm)", 100, locDeltaZMin, locDeltaZMax, 100, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaRVsDeltaE = new TH2F("dPluginHist_BCAL_DeltaRVsDeltaE", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaE (GeV);BCAL Shower #DeltaR (cm)", 100, locDeltaEMin, locDeltaEMax, 100, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaRVsDeltaT = new TH2F("dPluginHist_BCAL_DeltaRVsDeltaT", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaT (ns);BCAL Shower #DeltaR (cm)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaRMin, locDeltaRMax);
	dPluginHist_BCAL_DeltaPhiVsDeltaZ = new TH2F("dPluginHist_BCAL_DeltaPhiVsDeltaZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaZ (cm);BCAL Shower #Delta#phi (rad)", 100, locDeltaZMin, locDeltaZMax, 100, locDeltaPhiMin, locDeltaPhiMax);
	dPluginHist_BCAL_DeltaPhiVsDeltaE = new TH2F("dPluginHist_BCAL_DeltaPhiVsDeltaE", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaE (GeV);BCAL Shower #Delta#phi (rad)", 100, locDeltaEMin, locDeltaEMax, 100, locDeltaPhiMin, locDeltaPhiMax);
	dPluginHist_BCAL_DeltaPhiVsDeltaT = new TH2F("dPluginHist_BCAL_DeltaPhiVsDeltaT", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaT (ns);BCAL Shower #Delta#phi (rad)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaPhiMin, locDeltaPhiMax);
	dPluginHist_BCAL_DeltaZVsDeltaE = new TH2F("dPluginHist_BCAL_DeltaZVsDeltaE", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaE (GeV);BCAL Shower #DeltaZ (cm)", 100, locDeltaEMin, locDeltaEMax, 100, locDeltaZMin, locDeltaZMax);
	dPluginHist_BCAL_DeltaZVsDeltaT = new TH2F("dPluginHist_BCAL_DeltaZVsDeltaT", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaT (ns);BCAL Shower #DeltaZ (cm)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaZMin, locDeltaZMax);
	dPluginHist_BCAL_DeltaEVsDeltaT = new TH2F("dPluginHist_BCAL_DeltaEVsDeltaT", "Simulation Reconstructed Uncertainty Study;BCAL Shower #DeltaT (ns);BCAL Shower #DeltaE (GeV)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaEMin, locDeltaEMax);

	float locXMin = -75.0, locXMax = 75.0;
	float locYMin = -75.0, locYMax = 75.0;

	float locShowerSigmaXMin = 0.2, locShowerSigmaXMax = 1.0;
	float locShowerSigmaYMin = 0.2, locShowerSigmaYMax = 1.0;
	float locShowerSigmaZMin = 0.1, locShowerSigmaZMax = 1.1;
	float locShowerSigmaEMin = 0.0, locShowerSigmaEMax = 0.0;
	float locShowerSigmaTMin = 0.0, locShowerSigmaTMax = 2.0;

	//ShowerSigmaX Dependence
	dPluginHist_BCAL_ShowerSigmaXVsX = new TH2F("dPluginHist_BCAL_ShowerSigmaXVsX", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured X (cm);BCAL Shower #sigma X (cm)", 100, locXMin, locXMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_BCAL_ShowerSigmaXVsY = new TH2F("dPluginHist_BCAL_ShowerSigmaXVsY", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Y (cm);BCAL Shower #sigma X (cm)", 100, locYMin, locYMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_BCAL_ShowerSigmaXVsZ = new TH2F("dPluginHist_BCAL_ShowerSigmaXVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #sigma X (cm)", 100, locZMin, locZMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_BCAL_ShowerSigmaXVsE = new TH2F("dPluginHist_BCAL_ShowerSigmaXVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #sigma X (cm)", 100, locEMin, locEMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_BCAL_ShowerSigmaXVsT = new TH2F("dPluginHist_BCAL_ShowerSigmaXVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #sigma X (cm)", 100, locTMin, locTMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);

	//ShowerSigmaY Dependence
	dPluginHist_BCAL_ShowerSigmaYVsX = new TH2F("dPluginHist_BCAL_ShowerSigmaYVsX", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured X (cm);BCAL Shower #sigma Y (cm)", 100, locXMin, locXMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_BCAL_ShowerSigmaYVsY = new TH2F("dPluginHist_BCAL_ShowerSigmaYVsY", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Y (cm);BCAL Shower #sigma Y (cm)", 100, locYMin, locYMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_BCAL_ShowerSigmaYVsZ = new TH2F("dPluginHist_BCAL_ShowerSigmaYVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #sigma Y (cm)", 100, locZMin, locZMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_BCAL_ShowerSigmaYVsE = new TH2F("dPluginHist_BCAL_ShowerSigmaYVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #sigma Y (cm)", 100, locEMin, locEMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_BCAL_ShowerSigmaYVsT = new TH2F("dPluginHist_BCAL_ShowerSigmaYVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #sigma Y (cm)", 100, locTMin, locTMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);

	//ShowerSigmaZ Dependence
	dPluginHist_BCAL_ShowerSigmaZVsX = new TH2F("dPluginHist_BCAL_ShowerSigmaZVsX", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured X (cm);BCAL Shower #sigma Z (cm)", 100, locXMin, locXMax, 400, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_BCAL_ShowerSigmaZVsY = new TH2F("dPluginHist_BCAL_ShowerSigmaZVsY", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Y (cm);BCAL Shower #sigma Z (cm)", 100, locYMin, locYMax, 400, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_BCAL_ShowerSigmaZVsZ = new TH2F("dPluginHist_BCAL_ShowerSigmaZVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #sigma Z (cm)", 100, locZMin, locZMax, 400, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_BCAL_ShowerSigmaZVsE = new TH2F("dPluginHist_BCAL_ShowerSigmaZVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #sigma Z (cm)", 100, locEMin, locEMax, 400, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_BCAL_ShowerSigmaZVsT = new TH2F("dPluginHist_BCAL_ShowerSigmaZVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #sigma Z (cm)", 100, locTMin, locTMax, 400, locShowerSigmaZMin, locShowerSigmaZMax);

	//ShowerSigmaE Dependence
	dPluginHist_BCAL_ShowerSigmaEVsX = new TH2F("dPluginHist_BCAL_ShowerSigmaEVsX", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured X (cm);BCAL Shower #sigma E (GeV)", 100, locXMin, locXMax, 100, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_BCAL_ShowerSigmaEVsY = new TH2F("dPluginHist_BCAL_ShowerSigmaEVsY", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Y (cm);BCAL Shower #sigma E (GeV)", 100, locYMin, locYMax, 100, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_BCAL_ShowerSigmaEVsZ = new TH2F("dPluginHist_BCAL_ShowerSigmaEVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #sigma E (GeV)", 100, locZMin, locZMax, 100, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_BCAL_ShowerSigmaEVsE = new TH2F("dPluginHist_BCAL_ShowerSigmaEVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #sigma E (GeV)", 400, locEMin, locEMax, 400, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_BCAL_ShowerSigmaEVsT = new TH2F("dPluginHist_BCAL_ShowerSigmaEVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #sigma E (GeV)", 100, locTMin, locTMax, 100, locShowerSigmaEMin, locShowerSigmaEMax);

	//ShowerSigmaT Dependence
	dPluginHist_BCAL_ShowerSigmaTVsX = new TH2F("dPluginHist_BCAL_ShowerSigmaTVsX", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured X (cm);BCAL Shower #sigma T (ns)", 100, locXMin, locXMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_BCAL_ShowerSigmaTVsY = new TH2F("dPluginHist_BCAL_ShowerSigmaTVsY", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Y (cm);BCAL Shower #sigma T (ns)", 100, locYMin, locYMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_BCAL_ShowerSigmaTVsZ = new TH2F("dPluginHist_BCAL_ShowerSigmaTVsZ", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured Z (cm);BCAL Shower #sigma T (ns)", 100, locZMin, locZMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_BCAL_ShowerSigmaTVsE = new TH2F("dPluginHist_BCAL_ShowerSigmaTVsE", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured E (GeV);BCAL Shower #sigma T (ns)", 100, locEMin, locEMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_BCAL_ShowerSigmaTVsT = new TH2F("dPluginHist_BCAL_ShowerSigmaTVsT", "Simulation Reconstructed Uncertainty Study;BCAL Shower Measured T (ns);BCAL Shower #sigma T (ns)", 400, locTMin, locTMax, 400, locShowerSigmaTMin, locShowerSigmaTMax);

}

void BCALSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t BCALSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either BCALSelector::GetEntry() or TBranch::GetEntry()
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

	//Path Length Correction
	dPluginHist_BCAL_PathLengthCorrection->Fill(dTrueE, dPathLengthCorrection);
	float locPathLengthCorrectionPostE = dPathLengthCorrection - Calc_BCALPathLengthCorrection(dTrueE);
	dPluginHist_BCAL_PathLengthCorrectionPostEVsZ->Fill(dTrueZ, locPathLengthCorrectionPostE); //add (NOT subtract) the result of this fit to the path length

	bool locVsMeasuredFlag = true; //if true, swap true/meas
	if(locVsMeasuredFlag == true){
		//delta = meas - true
		dTrueR += dDeltaR;
		dTruePhi += dDeltaPhi;
		dTrueZ += dDeltaZ;
		dTrueE += dDeltaE;
		dTrueT += dDeltaT;
	}

	//DeltaR Dependence
	dPluginHist_BCAL_DeltaRVsR->Fill(dTrueR, dDeltaR);
	dPluginHist_BCAL_DeltaRVsPhi->Fill(dTruePhi, dDeltaR);
	dPluginHist_BCAL_DeltaRVsZ->Fill(dTrueZ, dDeltaR);
	dPluginHist_BCAL_DeltaRVsE->Fill(dTrueE, dDeltaR);
	dPluginHist_BCAL_DeltaRVsT->Fill(dTrueT, dDeltaR);

	//DeltaPhi Dependence
	dPluginHist_BCAL_DeltaPhiVsR->Fill(dTrueR, dDeltaPhi);
	dPluginHist_BCAL_DeltaPhiVsPhi->Fill(dTruePhi, dDeltaPhi);
	dPluginHist_BCAL_DeltaPhiVsZ->Fill(dTrueZ, dDeltaPhi);
	dPluginHist_BCAL_DeltaPhiVsE->Fill(dTrueE, dDeltaPhi);
	dPluginHist_BCAL_DeltaPhiVsT->Fill(dTrueT, dDeltaPhi);

	//DeltaZ Dependence
	dPluginHist_BCAL_DeltaZVsR->Fill(dTrueR, dDeltaZ);
	dPluginHist_BCAL_DeltaZVsPhi->Fill(dTruePhi, dDeltaZ);
	dPluginHist_BCAL_DeltaZVsZ->Fill(dTrueZ, dDeltaZ);
	dPluginHist_BCAL_DeltaZVsE->Fill(dTrueE, dDeltaZ);
	dPluginHist_BCAL_DeltaZVsT->Fill(dTrueT, dDeltaZ);

	//DeltaE Dependence
	dPluginHist_BCAL_DeltaEVsR->Fill(dTrueR, dDeltaE);
	dPluginHist_BCAL_DeltaEVsPhi->Fill(dTruePhi, dDeltaE);
	dPluginHist_BCAL_DeltaEVsZ->Fill(dTrueZ, dDeltaE);
	dPluginHist_BCAL_DeltaEVsE->Fill(dTrueE, dDeltaE);
	dPluginHist_BCAL_DeltaEVsT->Fill(dTrueT, dDeltaE);

	//DeltaT Dependence
	dPluginHist_BCAL_DeltaTVsR->Fill(dTrueR, dDeltaT);
	dPluginHist_BCAL_DeltaTVsPhi->Fill(dTruePhi, dDeltaT);
	dPluginHist_BCAL_DeltaTVsZ->Fill(dTrueZ, dDeltaT);
	dPluginHist_BCAL_DeltaTVsE->Fill(dTrueE, dDeltaT);
	dPluginHist_BCAL_DeltaTVsT->Fill(dTrueT, dDeltaT);

	//Common Dependence
	dPluginHist_BCAL_DeltaRVsDeltaPhi->Fill(dDeltaPhi, dDeltaR);
	dPluginHist_BCAL_DeltaRVsDeltaZ->Fill(dDeltaZ, dDeltaR);
	dPluginHist_BCAL_DeltaRVsDeltaE->Fill(dDeltaE, dDeltaR);
	dPluginHist_BCAL_DeltaRVsDeltaT->Fill(dDeltaT, dDeltaR);
	dPluginHist_BCAL_DeltaPhiVsDeltaZ->Fill(dDeltaZ, dDeltaPhi);
	dPluginHist_BCAL_DeltaPhiVsDeltaE->Fill(dDeltaE, dDeltaPhi);
	dPluginHist_BCAL_DeltaPhiVsDeltaT->Fill(dDeltaT, dDeltaPhi);
	dPluginHist_BCAL_DeltaZVsDeltaE->Fill(dDeltaE, dDeltaZ);
	dPluginHist_BCAL_DeltaZVsDeltaT->Fill(dDeltaT, dDeltaZ);
	dPluginHist_BCAL_DeltaEVsDeltaT->Fill(dDeltaT, dDeltaE);

	float dTrueX = dTrueR*cos(dTruePhi);
	float dTrueY = dTrueR*sin(dTruePhi);;

	//ShowerSigmaX Dependence
	dPluginHist_BCAL_ShowerSigmaXVsX->Fill(dTrueX, dShowerUncertaintyX);
	dPluginHist_BCAL_ShowerSigmaXVsY->Fill(dTrueY, dShowerUncertaintyX);
	dPluginHist_BCAL_ShowerSigmaXVsZ->Fill(dTrueZ, dShowerUncertaintyX);
	dPluginHist_BCAL_ShowerSigmaXVsE->Fill(dTrueE, dShowerUncertaintyX);
	dPluginHist_BCAL_ShowerSigmaXVsT->Fill(dTrueT, dShowerUncertaintyX);

	//ShowerSigmaY Dependence
	dPluginHist_BCAL_ShowerSigmaYVsX->Fill(dTrueX, dShowerUncertaintyY);
	dPluginHist_BCAL_ShowerSigmaYVsY->Fill(dTrueY, dShowerUncertaintyY);
	dPluginHist_BCAL_ShowerSigmaYVsZ->Fill(dTrueZ, dShowerUncertaintyY);
	dPluginHist_BCAL_ShowerSigmaYVsE->Fill(dTrueE, dShowerUncertaintyY);
	dPluginHist_BCAL_ShowerSigmaYVsT->Fill(dTrueT, dShowerUncertaintyY);

	//ShowerSigmaZ Dependence
	dPluginHist_BCAL_ShowerSigmaZVsX->Fill(dTrueX, dShowerUncertaintyZ);
	dPluginHist_BCAL_ShowerSigmaZVsY->Fill(dTrueY, dShowerUncertaintyZ);
	dPluginHist_BCAL_ShowerSigmaZVsZ->Fill(dTrueZ, dShowerUncertaintyZ);
	dPluginHist_BCAL_ShowerSigmaZVsE->Fill(dTrueE, dShowerUncertaintyZ);
	dPluginHist_BCAL_ShowerSigmaZVsT->Fill(dTrueT, dShowerUncertaintyZ);

	//ShowerSigmaE Dependence
	dPluginHist_BCAL_ShowerSigmaEVsX->Fill(dTrueX, dShowerUncertaintyE);
	dPluginHist_BCAL_ShowerSigmaEVsY->Fill(dTrueY, dShowerUncertaintyE);
	dPluginHist_BCAL_ShowerSigmaEVsZ->Fill(dTrueZ, dShowerUncertaintyE);
	dPluginHist_BCAL_ShowerSigmaEVsE->Fill(dTrueE, dShowerUncertaintyE);
	dPluginHist_BCAL_ShowerSigmaEVsT->Fill(dTrueT, dShowerUncertaintyE);

	//ShowerSigmaT Dependence
	dPluginHist_BCAL_ShowerSigmaTVsX->Fill(dTrueX, dShowerUncertaintyT);
	dPluginHist_BCAL_ShowerSigmaTVsY->Fill(dTrueY, dShowerUncertaintyT);
	dPluginHist_BCAL_ShowerSigmaTVsZ->Fill(dTrueZ, dShowerUncertaintyT);
	dPluginHist_BCAL_ShowerSigmaTVsE->Fill(dTrueE, dShowerUncertaintyT);
	dPluginHist_BCAL_ShowerSigmaTVsT->Fill(dTrueT, dShowerUncertaintyT);


   return kTRUE;
}

void BCALSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void BCALSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	dOutputFile->Write();
}

float BCALSelector::Calc_BCALPathLengthCorrection(float locEnergy){
	if(locEnergy >= 1.0){
		TF1 locFunction("df_BCAL_DepthCorrection", "[0] + [1]*x + [2]*exp([3]*x)", 0.0, 9.0);
		locFunction.SetParameters(9.95659, 0.142382, -2.9869, -0.56881);
		return locFunction.Eval(locEnergy);
	}
	TF1 locFunction("df_BCAL_DepthCorrection_LowE", "[0] + [1]*x + [2]*x*x", 0.0, 1.0);
	locFunction.SetParameters(4.89141, 5.98362, -2.56308);
	return locFunction.Eval(locEnergy);
}

