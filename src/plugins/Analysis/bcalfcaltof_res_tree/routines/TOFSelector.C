#define TOFSelector_cxx
// The class definition in TOFSelector.h has been generated automatically
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
// Root > T->Process("TOFSelector.C")
// Root > T->Process("TOFSelector.C","some options")
// Root > T->Process("TOFSelector.C+")
//

#include "TOFSelector.h"
#include <TH2.h>
#include <TStyle.h>


void TOFSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

	dOutputFile = new TFile("dh_TOFMCComparisonHists.root", "RECREATE");

	float locXMin = -100.0, locXMax = 100.0;
	float locYMin = -100.0, locYMax = 100.0;
	float locZMin = 616.0, locZMax = 622.0;
//	float locdEMin = 0.00392, locdEMax = 0.00405; //true
	float locdEMin = 0.0035, locdEMax = 0.020; //measured
	float locTMin = 18.0, locTMax = 23.0;
	float locBetaGammaMin = 0.0, locBetaGammaMax = 65.0; //pion
//	float locBetaGammaMin = 0.0, locBetaGammaMax = 10.0; //proton

	float locDeltaXMin = -4.0, locDeltaXMax = 4.0;
	float locDeltaYMin = -4.0, locDeltaYMax = 4.0;
	float locDeltaZMin = -0.4, locDeltaZMax = 0.4;
	float locDeltadEMin = 0.0, locDeltadEMax = 0.020;
	float locDeltaTMin = -0.2, locDeltaTMax = 0.2;

	int locNumXBins = 400, locNumDeltaXBins = 400;
	int locNumYBins = 400, locNumDeltaYBins = 400;
	int locNumZBins = 100, locNumDeltaZBins = 100;
	int locNumdEBins = 400, locNumDeltadEBins = 400;
	int locNumTBins = 400, locNumDeltaTBins = 400;
	int locNumBetaGammaBins = 400;

	//dEVsBetaGamma
	dPluginHist_TOF_dEVsBetaGamma = new TH2F("dPluginHist_TOF_dEVsBetaGamma", "Simulation Reconstructed Uncertainty Study;TOF Hit True #beta#gamma;TOF Hit dE (GeV)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEBins, locdEMin, locdEMax);

	//DeltaX Dependence
	dPluginHist_TOF_DeltaXVsX = new TH2F("dPluginHist_TOF_DeltaXVsX", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaX (cm)", locNumXBins, locXMin, locXMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsY = new TH2F("dPluginHist_TOF_DeltaXVsY", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaX (cm)", locNumYBins, locYMin, locYMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsZ = new TH2F("dPluginHist_TOF_DeltaXVsZ", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaX (cm)", locNumZBins, locZMin, locZMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsdE = new TH2F("dPluginHist_TOF_DeltaXVsdE", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaX (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsT = new TH2F("dPluginHist_TOF_DeltaXVsT", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaX (cm)", locNumTBins, locTMin, locTMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);

	//DeltaY Dependence
	dPluginHist_TOF_DeltaYVsX = new TH2F("dPluginHist_TOF_DeltaYVsX", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaY (cm)", locNumXBins, locXMin, locXMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsY = new TH2F("dPluginHist_TOF_DeltaYVsY", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaY (cm)", locNumYBins, locYMin, locYMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsZ = new TH2F("dPluginHist_TOF_DeltaYVsZ", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaY (cm)", locNumZBins, locZMin, locZMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsdE = new TH2F("dPluginHist_TOF_DeltaYVsdE", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaY (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsT = new TH2F("dPluginHist_TOF_DeltaYVsT", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaY (cm)", locNumTBins, locTMin, locTMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);

	//DeltaZ Dependence
	dPluginHist_TOF_DeltaZVsX = new TH2F("dPluginHist_TOF_DeltaZVsX", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaZ (cm)", locNumXBins, locXMin, locXMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsY = new TH2F("dPluginHist_TOF_DeltaZVsY", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaZ (cm)", locNumYBins, locYMin, locYMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsZ = new TH2F("dPluginHist_TOF_DeltaZVsZ", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaZ (cm)", locNumZBins, locZMin, locZMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsdE = new TH2F("dPluginHist_TOF_DeltaZVsdE", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaZ (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsT = new TH2F("dPluginHist_TOF_DeltaZVsT", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaZ (cm)", locNumTBins, locTMin, locTMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);

	//DeltadE Dependence
	dPluginHist_TOF_DeltadEVsX = new TH2F("dPluginHist_TOF_DeltadEVsX", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltadE (GeV)", locNumXBins, locXMin, locXMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsY = new TH2F("dPluginHist_TOF_DeltadEVsY", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltadE (GeV)", locNumYBins, locYMin, locYMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsZ = new TH2F("dPluginHist_TOF_DeltadEVsZ", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltadE (GeV)", locNumZBins, locZMin, locZMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsdE = new TH2F("dPluginHist_TOF_DeltadEVsdE", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltadE (GeV)", locNumdEBins, locdEMin, locdEMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsT = new TH2F("dPluginHist_TOF_DeltadEVsT", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltadE (GeV)", locNumTBins, locTMin, locTMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);

	//DeltaT Dependence
	dPluginHist_TOF_DeltaTVsX = new TH2F("dPluginHist_TOF_DeltaTVsX", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaT (ns)", locNumXBins, locXMin, locXMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsY = new TH2F("dPluginHist_TOF_DeltaTVsY", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaT (ns)", locNumYBins, locYMin, locYMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsZ = new TH2F("dPluginHist_TOF_DeltaTVsZ", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaT (ns)", locNumZBins, locZMin, locZMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsdE = new TH2F("dPluginHist_TOF_DeltaTVsdE", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaT (ns)", locNumdEBins, locdEMin, locdEMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsT = new TH2F("dPluginHist_TOF_DeltaTVsT", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaT (ns)", locNumTBins, locTMin, locTMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);

	//Common Dependence
	dPluginHist_TOF_DeltaXVsDeltaY = new TH2F("dPluginHist_TOF_DeltaXVsDeltaY", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaY (cm);TOF Hit #DeltaX (cm)", locNumDeltaYBins, locDeltaYMin, locDeltaYMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltaZ = new TH2F("dPluginHist_TOF_DeltaXVsDeltaZ", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaZ (cm);TOF Hit #DeltaX (cm)", locNumDeltaZBins, locDeltaZMin, locDeltaZMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltadE = new TH2F("dPluginHist_TOF_DeltaXVsDeltadE", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaX (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltaT = new TH2F("dPluginHist_TOF_DeltaXVsDeltaT", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaX (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaYVsDeltaZ = new TH2F("dPluginHist_TOF_DeltaYVsDeltaZ", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaZ (cm);TOF Hit #DeltaY (cm)", locNumDeltaZBins, locDeltaZMin, locDeltaZMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsDeltadE = new TH2F("dPluginHist_TOF_DeltaYVsDeltadE", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaY (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsDeltaT = new TH2F("dPluginHist_TOF_DeltaYVsDeltaT", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaY (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaZVsDeltadE = new TH2F("dPluginHist_TOF_DeltaZVsDeltadE", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaZ (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsDeltaT = new TH2F("dPluginHist_TOF_DeltaZVsDeltaT", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaZ (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltadEVsDeltaT = new TH2F("dPluginHist_TOF_DeltadEVsDeltaT", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltadE (GeV)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);

	//dEVsBetaGamma
	dPluginHist_TOF_dEVsBetaGamma_HorizontalOnly = new TH2F("dPluginHist_TOF_dEVsBetaGamma_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit True #beta#gamma;TOF Hit dE (GeV)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEBins, locdEMin, locdEMax);

	//DeltaX Dependence
	dPluginHist_TOF_DeltaXVsX_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsX_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaX (cm)", locNumXBins, locXMin, locXMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsY_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsY_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaX (cm)", locNumYBins, locYMin, locYMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsZ_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsZ_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaX (cm)", locNumZBins, locZMin, locZMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsdE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsdE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaX (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaX (cm)", locNumTBins, locTMin, locTMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);

	//DeltaY Dependence
	dPluginHist_TOF_DeltaYVsX_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsX_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaY (cm)", locNumXBins, locXMin, locXMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsY_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsY_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaY (cm)", locNumYBins, locYMin, locYMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsZ_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsZ_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaY (cm)", locNumZBins, locZMin, locZMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsdE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsdE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaY (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaY (cm)", locNumTBins, locTMin, locTMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);

	//DeltaZ Dependence
	dPluginHist_TOF_DeltaZVsX_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaZVsX_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaZ (cm)", locNumXBins, locXMin, locXMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsY_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaZVsY_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaZ (cm)", locNumYBins, locYMin, locYMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsZ_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaZVsZ_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaZ (cm)", locNumZBins, locZMin, locZMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsdE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaZVsdE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaZ (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaZVsT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaZ (cm)", locNumTBins, locTMin, locTMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);

	//DeltadE Dependence
	dPluginHist_TOF_DeltadEVsX_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltadEVsX_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltadE (GeV)", locNumXBins, locXMin, locXMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsY_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltadEVsY_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltadE (GeV)", locNumYBins, locYMin, locYMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsZ_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltadEVsZ_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltadE (GeV)", locNumZBins, locZMin, locZMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsdE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltadEVsdE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltadE (GeV)", locNumdEBins, locdEMin, locdEMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltadEVsT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltadE (GeV)", locNumTBins, locTMin, locTMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);

	//DeltaT Dependence
	dPluginHist_TOF_DeltaTVsX_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaTVsX_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaT (ns)", locNumXBins, locXMin, locXMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsY_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaTVsY_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaT (ns)", locNumYBins, locYMin, locYMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsZ_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaTVsZ_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaT (ns)", locNumZBins, locZMin, locZMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsdE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaTVsdE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaT (ns)", locNumdEBins, locdEMin, locdEMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaTVsT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaT (ns)", locNumTBins, locTMin, locTMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);

	//Common Dependence
	dPluginHist_TOF_DeltaXVsDeltaY_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltaY_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaY (cm);TOF Hit #DeltaX (cm)", locNumDeltaYBins, locDeltaYMin, locDeltaYMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltaZ_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltaZ_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaZ (cm);TOF Hit #DeltaX (cm)", locNumDeltaZBins, locDeltaZMin, locDeltaZMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltadE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltadE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaX (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltaT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltaT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaX (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaYVsDeltaZ_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsDeltaZ_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaZ (cm);TOF Hit #DeltaY (cm)", locNumDeltaZBins, locDeltaZMin, locDeltaZMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsDeltadE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsDeltadE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaY (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsDeltaT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaYVsDeltaT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaY (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaZVsDeltadE_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaZVsDeltadE_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaZ (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsDeltaT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltaZVsDeltaT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaZ (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltadEVsDeltaT_HorizontalOnly = new TH2F("dPluginHist_TOF_DeltadEVsDeltaT_HorizontalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltadE (GeV)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);

	//dEVsBetaGamma
	dPluginHist_TOF_dEVsBetaGamma_VerticalOnly = new TH2F("dPluginHist_TOF_dEVsBetaGamma_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit True #beta#gamma;TOF Hit dE (GeV)", locNumBetaGammaBins, locBetaGammaMin, locBetaGammaMax, locNumdEBins, locdEMin, locdEMax);

	//DeltaX Dependence
	dPluginHist_TOF_DeltaXVsX_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsX_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaX (cm)", locNumXBins, locXMin, locXMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsY_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsY_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaX (cm)", locNumYBins, locYMin, locYMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsZ_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsZ_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaX (cm)", locNumZBins, locZMin, locZMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsdE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsdE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaX (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaX (cm)", locNumTBins, locTMin, locTMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);

	//DeltaY Dependence
	dPluginHist_TOF_DeltaYVsX_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsX_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaY (cm)", locNumXBins, locXMin, locXMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsY_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsY_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaY (cm)", locNumYBins, locYMin, locYMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsZ_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsZ_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaY (cm)", locNumZBins, locZMin, locZMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsdE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsdE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaY (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaY (cm)", locNumTBins, locTMin, locTMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);

	//DeltaZ Dependence
	dPluginHist_TOF_DeltaZVsX_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaZVsX_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaZ (cm)", locNumXBins, locXMin, locXMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsY_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaZVsY_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaZ (cm)", locNumYBins, locYMin, locYMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsZ_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaZVsZ_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaZ (cm)", locNumZBins, locZMin, locZMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsdE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaZVsdE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaZ (cm)", locNumdEBins, locdEMin, locdEMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaZVsT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaZ (cm)", locNumTBins, locTMin, locTMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);

	//DeltadE Dependence
	dPluginHist_TOF_DeltadEVsX_VerticalOnly = new TH2F("dPluginHist_TOF_DeltadEVsX_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltadE (GeV)", locNumXBins, locXMin, locXMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsY_VerticalOnly = new TH2F("dPluginHist_TOF_DeltadEVsY_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltadE (GeV)", locNumYBins, locYMin, locYMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsZ_VerticalOnly = new TH2F("dPluginHist_TOF_DeltadEVsZ_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltadE (GeV)", locNumZBins, locZMin, locZMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsdE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltadEVsdE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltadE (GeV)", locNumdEBins, locdEMin, locdEMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
	dPluginHist_TOF_DeltadEVsT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltadEVsT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltadE (GeV)", locNumTBins, locTMin, locTMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);

	//DeltaT Dependence
	dPluginHist_TOF_DeltaTVsX_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaTVsX_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured X (cm);TOF Hit #DeltaT (ns)", locNumXBins, locXMin, locXMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsY_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaTVsY_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Y (cm);TOF Hit #DeltaT (ns)", locNumYBins, locYMin, locYMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsZ_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaTVsZ_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured Z (cm);TOF Hit #DeltaT (ns)", locNumZBins, locZMin, locZMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsdE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaTVsdE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured dE (GeV);TOF Hit #DeltaT (ns)", locNumdEBins, locdEMin, locdEMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);
	dPluginHist_TOF_DeltaTVsT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaTVsT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit Measured T (ns);TOF Hit #DeltaT (ns)", locNumTBins, locTMin, locTMax, locNumDeltaTBins, locDeltaTMin, locDeltaTMax);

	//Common Dependence
	dPluginHist_TOF_DeltaXVsDeltaY_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltaY_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaY (cm);TOF Hit #DeltaX (cm)", locNumDeltaYBins, locDeltaYMin, locDeltaYMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltaZ_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltaZ_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaZ (cm);TOF Hit #DeltaX (cm)", locNumDeltaZBins, locDeltaZMin, locDeltaZMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltadE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltadE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaX (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaXVsDeltaT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaXVsDeltaT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaX (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaXBins, locDeltaXMin, locDeltaXMax);
	dPluginHist_TOF_DeltaYVsDeltaZ_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsDeltaZ_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaZ (cm);TOF Hit #DeltaY (cm)", locNumDeltaZBins, locDeltaZMin, locDeltaZMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsDeltadE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsDeltadE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaY (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaYVsDeltaT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaYVsDeltaT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaY (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaYBins, locDeltaYMin, locDeltaYMax);
	dPluginHist_TOF_DeltaZVsDeltadE_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaZVsDeltadE_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltadE (GeV);TOF Hit #DeltaZ (cm)", locNumDeltadEBins, locDeltadEMin, locDeltadEMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltaZVsDeltaT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltaZVsDeltaT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltaZ (cm)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltaZBins, locDeltaZMin, locDeltaZMax);
	dPluginHist_TOF_DeltadEVsDeltaT_VerticalOnly = new TH2F("dPluginHist_TOF_DeltadEVsDeltaT_VerticalOnly", "Simulation Reconstructed Uncertainty Study;TOF Hit #DeltaT (ns);TOF Hit #DeltadE (GeV)", locNumDeltaTBins, locDeltaTMin, locDeltaTMax, locNumDeltadEBins, locDeltadEMin, locDeltadEMax);
}

void TOFSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t TOFSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either TOFSelector::GetEntry() or TBranch::GetEntry()
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

	bool locVsMeasuredFlag = true; //if true, swap true/meas
	if(locVsMeasuredFlag == true){
		//delta = meas - true
		dTrueX += dDeltaX;
		dTrueY += dDeltaY;
		dTrueZ += dDeltaZ;
		dTruedE += dDeltadE;
		dTrueT += dDeltaT;
	}

	if((dHorizontalPlaneFlag == true) && (dVerticalPlaneFlag == true)){

	//dEVsBetaGamma Dependence
	dPluginHist_TOF_dEVsBetaGamma->Fill(dTrueBetaGamma, dTruedE); //dTruedE is measured value if switched above

	//DeltaX Dependence
	dPluginHist_TOF_DeltaXVsX->Fill(dTrueX, dDeltaX);
	dPluginHist_TOF_DeltaXVsY->Fill(dTrueY, dDeltaX);
	dPluginHist_TOF_DeltaXVsZ->Fill(dTrueZ, dDeltaX);
	dPluginHist_TOF_DeltaXVsdE->Fill(dTruedE, dDeltaX);
	dPluginHist_TOF_DeltaXVsT->Fill(dTrueT, dDeltaX);

	//DeltaY Dependence
	dPluginHist_TOF_DeltaYVsX->Fill(dTrueX, dDeltaY);
	dPluginHist_TOF_DeltaYVsY->Fill(dTrueY, dDeltaY);
	dPluginHist_TOF_DeltaYVsZ->Fill(dTrueZ, dDeltaY);
	dPluginHist_TOF_DeltaYVsdE->Fill(dTruedE, dDeltaY);
	dPluginHist_TOF_DeltaYVsT->Fill(dTrueT, dDeltaY);

	//DeltaZ Dependence
	dPluginHist_TOF_DeltaZVsX->Fill(dTrueX, dDeltaZ);
	dPluginHist_TOF_DeltaZVsY->Fill(dTrueY, dDeltaZ);
	dPluginHist_TOF_DeltaZVsZ->Fill(dTrueZ, dDeltaZ);
	dPluginHist_TOF_DeltaZVsdE->Fill(dTruedE, dDeltaZ);
	dPluginHist_TOF_DeltaZVsT->Fill(dTrueT, dDeltaZ);

	//DeltadE Dependence
	dPluginHist_TOF_DeltadEVsX->Fill(dTrueX, dDeltadE);
	dPluginHist_TOF_DeltadEVsY->Fill(dTrueY, dDeltadE);
	dPluginHist_TOF_DeltadEVsZ->Fill(dTrueZ, dDeltadE);
	dPluginHist_TOF_DeltadEVsdE->Fill(dTruedE, dDeltadE);
	dPluginHist_TOF_DeltadEVsT->Fill(dTrueT, dDeltadE);

	//DeltaT Dependence
	dPluginHist_TOF_DeltaTVsX->Fill(dTrueX, dDeltaT);
	dPluginHist_TOF_DeltaTVsY->Fill(dTrueY, dDeltaT);
	dPluginHist_TOF_DeltaTVsZ->Fill(dTrueZ, dDeltaT);
	dPluginHist_TOF_DeltaTVsdE->Fill(dTruedE, dDeltaT);
	dPluginHist_TOF_DeltaTVsT->Fill(dTrueT, dDeltaT);

	//Common Dependence
	dPluginHist_TOF_DeltaXVsDeltaY->Fill(dDeltaY, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltaZ->Fill(dDeltaZ, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltadE->Fill(dDeltadE, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltaT->Fill(dDeltaT, dDeltaX);
	dPluginHist_TOF_DeltaYVsDeltaZ->Fill(dDeltaZ, dDeltaY);
	dPluginHist_TOF_DeltaYVsDeltadE->Fill(dDeltadE, dDeltaY);
	dPluginHist_TOF_DeltaYVsDeltaT->Fill(dDeltaT, dDeltaY);
	dPluginHist_TOF_DeltaZVsDeltadE->Fill(dDeltadE, dDeltaZ);
	dPluginHist_TOF_DeltaZVsDeltaT->Fill(dDeltaT, dDeltaZ);
	dPluginHist_TOF_DeltadEVsDeltaT->Fill(dDeltaT, dDeltadE);
	}

	if((dHorizontalPlaneFlag == true) && (dVerticalPlaneFlag == false)){
	//dEVsBetaGamma Dependence
	dPluginHist_TOF_dEVsBetaGamma_HorizontalOnly->Fill(dTrueBetaGamma, dTruedE); //dTruedE is measured value if switched above

	//DeltaX Dependence
	dPluginHist_TOF_DeltaXVsX_HorizontalOnly->Fill(dTrueX, dDeltaX);
	dPluginHist_TOF_DeltaXVsY_HorizontalOnly->Fill(dTrueY, dDeltaX);
	dPluginHist_TOF_DeltaXVsZ_HorizontalOnly->Fill(dTrueZ, dDeltaX);
	dPluginHist_TOF_DeltaXVsdE_HorizontalOnly->Fill(dTruedE, dDeltaX);
	dPluginHist_TOF_DeltaXVsT_HorizontalOnly->Fill(dTrueT, dDeltaX);

	//DeltaY Dependence
	dPluginHist_TOF_DeltaYVsX_HorizontalOnly->Fill(dTrueX, dDeltaY);
	dPluginHist_TOF_DeltaYVsY_HorizontalOnly->Fill(dTrueY, dDeltaY);
	dPluginHist_TOF_DeltaYVsZ_HorizontalOnly->Fill(dTrueZ, dDeltaY);
	dPluginHist_TOF_DeltaYVsdE_HorizontalOnly->Fill(dTruedE, dDeltaY);
	dPluginHist_TOF_DeltaYVsT_HorizontalOnly->Fill(dTrueT, dDeltaY);

	//DeltaZ Dependence
	dPluginHist_TOF_DeltaZVsX_HorizontalOnly->Fill(dTrueX, dDeltaZ);
	dPluginHist_TOF_DeltaZVsY_HorizontalOnly->Fill(dTrueY, dDeltaZ);
	dPluginHist_TOF_DeltaZVsZ_HorizontalOnly->Fill(dTrueZ, dDeltaZ);
	dPluginHist_TOF_DeltaZVsdE_HorizontalOnly->Fill(dTruedE, dDeltaZ);
	dPluginHist_TOF_DeltaZVsT_HorizontalOnly->Fill(dTrueT, dDeltaZ);

	//DeltadE Dependence
	dPluginHist_TOF_DeltadEVsX_HorizontalOnly->Fill(dTrueX, dDeltadE);
	dPluginHist_TOF_DeltadEVsY_HorizontalOnly->Fill(dTrueY, dDeltadE);
	dPluginHist_TOF_DeltadEVsZ_HorizontalOnly->Fill(dTrueZ, dDeltadE);
	dPluginHist_TOF_DeltadEVsdE_HorizontalOnly->Fill(dTruedE, dDeltadE);
	dPluginHist_TOF_DeltadEVsT_HorizontalOnly->Fill(dTrueT, dDeltadE);

	//DeltaT Dependence
	dPluginHist_TOF_DeltaTVsX_HorizontalOnly->Fill(dTrueX, dDeltaT);
	dPluginHist_TOF_DeltaTVsY_HorizontalOnly->Fill(dTrueY, dDeltaT);
	dPluginHist_TOF_DeltaTVsZ_HorizontalOnly->Fill(dTrueZ, dDeltaT);
	dPluginHist_TOF_DeltaTVsdE_HorizontalOnly->Fill(dTruedE, dDeltaT);
	dPluginHist_TOF_DeltaTVsT_HorizontalOnly->Fill(dTrueT, dDeltaT);

	//Common Dependence
	dPluginHist_TOF_DeltaXVsDeltaY_HorizontalOnly->Fill(dDeltaY, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltaZ_HorizontalOnly->Fill(dDeltaZ, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltadE_HorizontalOnly->Fill(dDeltadE, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltaT_HorizontalOnly->Fill(dDeltaT, dDeltaX);
	dPluginHist_TOF_DeltaYVsDeltaZ_HorizontalOnly->Fill(dDeltaZ, dDeltaY);
	dPluginHist_TOF_DeltaYVsDeltadE_HorizontalOnly->Fill(dDeltadE, dDeltaY);
	dPluginHist_TOF_DeltaYVsDeltaT_HorizontalOnly->Fill(dDeltaT, dDeltaY);
	dPluginHist_TOF_DeltaZVsDeltadE_HorizontalOnly->Fill(dDeltadE, dDeltaZ);
	dPluginHist_TOF_DeltaZVsDeltaT_HorizontalOnly->Fill(dDeltaT, dDeltaZ);
	dPluginHist_TOF_DeltadEVsDeltaT_HorizontalOnly->Fill(dDeltaT, dDeltadE);
	}

	if((dHorizontalPlaneFlag == false) && (dVerticalPlaneFlag == true)){
	//dEVsBetaGamma Dependence
	dPluginHist_TOF_dEVsBetaGamma_VerticalOnly->Fill(dTrueBetaGamma, dTruedE); //dTruedE is measured value if switched above

	//DeltaX Dependence
	dPluginHist_TOF_DeltaXVsX_VerticalOnly->Fill(dTrueX, dDeltaX);
	dPluginHist_TOF_DeltaXVsY_VerticalOnly->Fill(dTrueY, dDeltaX);
	dPluginHist_TOF_DeltaXVsZ_VerticalOnly->Fill(dTrueZ, dDeltaX);
	dPluginHist_TOF_DeltaXVsdE_VerticalOnly->Fill(dTruedE, dDeltaX);
	dPluginHist_TOF_DeltaXVsT_VerticalOnly->Fill(dTrueT, dDeltaX);

	//DeltaY Dependence
	dPluginHist_TOF_DeltaYVsX_VerticalOnly->Fill(dTrueX, dDeltaY);
	dPluginHist_TOF_DeltaYVsY_VerticalOnly->Fill(dTrueY, dDeltaY);
	dPluginHist_TOF_DeltaYVsZ_VerticalOnly->Fill(dTrueZ, dDeltaY);
	dPluginHist_TOF_DeltaYVsdE_VerticalOnly->Fill(dTruedE, dDeltaY);
	dPluginHist_TOF_DeltaYVsT_VerticalOnly->Fill(dTrueT, dDeltaY);

	//DeltaZ Dependence
	dPluginHist_TOF_DeltaZVsX_VerticalOnly->Fill(dTrueX, dDeltaZ);
	dPluginHist_TOF_DeltaZVsY_VerticalOnly->Fill(dTrueY, dDeltaZ);
	dPluginHist_TOF_DeltaZVsZ_VerticalOnly->Fill(dTrueZ, dDeltaZ);
	dPluginHist_TOF_DeltaZVsdE_VerticalOnly->Fill(dTruedE, dDeltaZ);
	dPluginHist_TOF_DeltaZVsT_VerticalOnly->Fill(dTrueT, dDeltaZ);

	//DeltadE Dependence
	dPluginHist_TOF_DeltadEVsX_VerticalOnly->Fill(dTrueX, dDeltadE);
	dPluginHist_TOF_DeltadEVsY_VerticalOnly->Fill(dTrueY, dDeltadE);
	dPluginHist_TOF_DeltadEVsZ_VerticalOnly->Fill(dTrueZ, dDeltadE);
	dPluginHist_TOF_DeltadEVsdE_VerticalOnly->Fill(dTruedE, dDeltadE);
	dPluginHist_TOF_DeltadEVsT_VerticalOnly->Fill(dTrueT, dDeltadE);

	//DeltaT Dependence
	dPluginHist_TOF_DeltaTVsX_VerticalOnly->Fill(dTrueX, dDeltaT);
	dPluginHist_TOF_DeltaTVsY_VerticalOnly->Fill(dTrueY, dDeltaT);
	dPluginHist_TOF_DeltaTVsZ_VerticalOnly->Fill(dTrueZ, dDeltaT);
	dPluginHist_TOF_DeltaTVsdE_VerticalOnly->Fill(dTruedE, dDeltaT);
	dPluginHist_TOF_DeltaTVsT_VerticalOnly->Fill(dTrueT, dDeltaT);

	//Common Dependence
	dPluginHist_TOF_DeltaXVsDeltaY_VerticalOnly->Fill(dDeltaY, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltaZ_VerticalOnly->Fill(dDeltaZ, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltadE_VerticalOnly->Fill(dDeltadE, dDeltaX);
	dPluginHist_TOF_DeltaXVsDeltaT_VerticalOnly->Fill(dDeltaT, dDeltaX);
	dPluginHist_TOF_DeltaYVsDeltaZ_VerticalOnly->Fill(dDeltaZ, dDeltaY);
	dPluginHist_TOF_DeltaYVsDeltadE_VerticalOnly->Fill(dDeltadE, dDeltaY);
	dPluginHist_TOF_DeltaYVsDeltaT_VerticalOnly->Fill(dDeltaT, dDeltaY);
	dPluginHist_TOF_DeltaZVsDeltadE_VerticalOnly->Fill(dDeltadE, dDeltaZ);
	dPluginHist_TOF_DeltaZVsDeltaT_VerticalOnly->Fill(dDeltaT, dDeltaZ);
	dPluginHist_TOF_DeltadEVsDeltaT_VerticalOnly->Fill(dDeltaT, dDeltadE);
	}

   return kTRUE;
}

void TOFSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void TOFSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	dOutputFile->Write();
}
