#define FCALSelector_cxx
// The class definition in FCALSelector.h has been generated automatically
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
// Root > T->Process("FCALSelector.C")
// Root > T->Process("FCALSelector.C","some options")
// Root > T->Process("FCALSelector.C+")
//

#include "FCALSelector.h"
#include <TH2.h>
#include <TStyle.h>


void FCALSelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

	dOutputFile = new TFile("dh_FCALMCComparisonHists.root", "RECREATE");

	float locXMin = -100.0, locXMax = 100.0;
	float locYMin = -100.0, locYMax = 100.0;
	float locZMin = 632.0, locZMax = 646.0;
	float locEMin = 0.0, locEMax = 9.1;
	float locTMin = 18.8, locTMax = 19.8;

	float locDeltaXMin = -8.0, locDeltaXMax = 8.0;
	float locDeltaYMin = -8.0, locDeltaYMax = 8.0;
//	float locDeltaZMin = -2.0, locDeltaZMax = 12.0; //pre depth correction
	float locDeltaZMin = -3.0, locDeltaZMax = 1.5;
	float locDeltaEMin = -1.5, locDeltaEMax = 1.5;
	float locDeltaTMin = -5.0, locDeltaTMax = 5.0;

	dPluginHist_FCAL_PathLengthCorrection = new TH2F("dPluginHist_FCAL_PathLengthCorrection", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);Path Length Correction (cm)", 100, locEMin, locEMax, 1000, -1.0, 25.0);

	//DeltaX Dependence
	dPluginHist_FCAL_DeltaXVsX = new TH2F("dPluginHist_FCAL_DeltaXVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #DeltaX (cm)", 500, locXMin, locXMax, 400, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaXVsY = new TH2F("dPluginHist_FCAL_DeltaXVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #DeltaX (cm)", 500, locYMin, locYMax, 400, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaXVsZ = new TH2F("dPluginHist_FCAL_DeltaXVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #DeltaX (cm)", 100, locZMin, locZMax, 400, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaXVsE = new TH2F("dPluginHist_FCAL_DeltaXVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #DeltaX (cm)", 100, locEMin, locEMax, 400, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaXVsT = new TH2F("dPluginHist_FCAL_DeltaXVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #DeltaX (cm)", 100, locTMin, locTMax, 400, locDeltaXMin, locDeltaXMax);

	//DeltaY Dependence
	dPluginHist_FCAL_DeltaYVsX = new TH2F("dPluginHist_FCAL_DeltaYVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #DeltaY (cm)", 500, locXMin, locXMax, 400, locDeltaYMin, locDeltaYMax);
	dPluginHist_FCAL_DeltaYVsY = new TH2F("dPluginHist_FCAL_DeltaYVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #DeltaY (cm)", 500, locYMin, locYMax, 400, locDeltaYMin, locDeltaYMax);
	dPluginHist_FCAL_DeltaYVsZ = new TH2F("dPluginHist_FCAL_DeltaYVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #DeltaY (cm)", 100, locZMin, locZMax, 400, locDeltaYMin, locDeltaYMax);
	dPluginHist_FCAL_DeltaYVsE = new TH2F("dPluginHist_FCAL_DeltaYVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #DeltaY (cm)", 100, locEMin, locEMax, 400, locDeltaYMin, locDeltaYMax);
	dPluginHist_FCAL_DeltaYVsT = new TH2F("dPluginHist_FCAL_DeltaYVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #DeltaY (cm)", 100, locTMin, locTMax, 400, locDeltaYMin, locDeltaYMax);

	//DeltaZ Dependence
	dPluginHist_FCAL_DeltaZVsX = new TH2F("dPluginHist_FCAL_DeltaZVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #DeltaZ (cm)", 100, locXMin, locXMax, 400, locDeltaZMin, locDeltaZMax);
	dPluginHist_FCAL_DeltaZVsY = new TH2F("dPluginHist_FCAL_DeltaZVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #DeltaZ (cm)", 100, locYMin, locYMax, 400, locDeltaZMin, locDeltaZMax);
	dPluginHist_FCAL_DeltaZVsZ = new TH2F("dPluginHist_FCAL_DeltaZVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #DeltaZ (cm)", 100, locZMin, locZMax, 400, locDeltaZMin, locDeltaZMax);
	dPluginHist_FCAL_DeltaZVsE = new TH2F("dPluginHist_FCAL_DeltaZVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #DeltaZ (cm)", 100, locEMin, locEMax, 400, locDeltaZMin, locDeltaZMax);
	dPluginHist_FCAL_DeltaZVsT = new TH2F("dPluginHist_FCAL_DeltaZVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #DeltaZ (cm)", 100, locTMin, locTMax, 400, locDeltaZMin, locDeltaZMax);

	//DeltaE Dependence
	dPluginHist_FCAL_DeltaEVsX = new TH2F("dPluginHist_FCAL_DeltaEVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #DeltaE (GeV)", 100, locXMin, locXMax, 400, locDeltaEMin, locDeltaEMax);
	dPluginHist_FCAL_DeltaEVsY = new TH2F("dPluginHist_FCAL_DeltaEVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #DeltaE (GeV)", 100, locYMin, locYMax, 400, locDeltaEMin, locDeltaEMax);
	dPluginHist_FCAL_DeltaEVsZ = new TH2F("dPluginHist_FCAL_DeltaEVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #DeltaE (GeV)", 100, locZMin, locZMax, 400, locDeltaEMin, locDeltaEMax);
	dPluginHist_FCAL_DeltaEVsE = new TH2F("dPluginHist_FCAL_DeltaEVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #DeltaE (GeV)", 100, locEMin, locEMax, 400, locDeltaEMin, locDeltaEMax);
	dPluginHist_FCAL_DeltaEVsT = new TH2F("dPluginHist_FCAL_DeltaEVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #DeltaE (GeV)", 100, locTMin, locTMax, 400, locDeltaEMin, locDeltaEMax);

	//DeltaT Dependence
	dPluginHist_FCAL_DeltaTVsX = new TH2F("dPluginHist_FCAL_DeltaTVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #DeltaT (ns)", 100, locXMin, locXMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_FCAL_DeltaTVsY = new TH2F("dPluginHist_FCAL_DeltaTVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #DeltaT (ns)", 100, locYMin, locYMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_FCAL_DeltaTVsZ = new TH2F("dPluginHist_FCAL_DeltaTVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #DeltaT (ns)", 100, locZMin, locZMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_FCAL_DeltaTVsE = new TH2F("dPluginHist_FCAL_DeltaTVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #DeltaT (ns)", 200, locEMin, locEMax, 400, locDeltaTMin, locDeltaTMax);
	dPluginHist_FCAL_DeltaTVsT = new TH2F("dPluginHist_FCAL_DeltaTVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #DeltaT (ns)", 100, locTMin, locTMax, 400, locDeltaTMin, locDeltaTMax);

	//Common Dependence
	dPluginHist_FCAL_DeltaXVsDeltaY = new TH2F("dPluginHist_FCAL_DeltaXVsDeltaY", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaY (cm);FCAL Shower #DeltaX (cm)", 100, locDeltaYMin, locDeltaYMax, 100, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaXVsDeltaZ = new TH2F("dPluginHist_FCAL_DeltaXVsDeltaZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaZ (cm);FCAL Shower #DeltaX (cm)", 100, locDeltaZMin, locDeltaZMax, 100, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaXVsDeltaE = new TH2F("dPluginHist_FCAL_DeltaXVsDeltaE", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaE (GeV);FCAL Shower #DeltaX (cm)", 100, locDeltaEMin, locDeltaEMax, 100, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaXVsDeltaT = new TH2F("dPluginHist_FCAL_DeltaXVsDeltaT", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaT (ns);FCAL Shower #DeltaX (cm)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaXMin, locDeltaXMax);
	dPluginHist_FCAL_DeltaYVsDeltaZ = new TH2F("dPluginHist_FCAL_DeltaYVsDeltaZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaZ (cm);FCAL Shower #DeltaY (cm)", 100, locDeltaZMin, locDeltaZMax, 100, locDeltaYMin, locDeltaYMax);
	dPluginHist_FCAL_DeltaYVsDeltaE = new TH2F("dPluginHist_FCAL_DeltaYVsDeltaE", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaE (GeV);FCAL Shower #DeltaY (cm)", 100, locDeltaEMin, locDeltaEMax, 100, locDeltaYMin, locDeltaYMax);
	dPluginHist_FCAL_DeltaYVsDeltaT = new TH2F("dPluginHist_FCAL_DeltaYVsDeltaT", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaT (ns);FCAL Shower #DeltaY (cm)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaYMin, locDeltaYMax);
	dPluginHist_FCAL_DeltaZVsDeltaE = new TH2F("dPluginHist_FCAL_DeltaZVsDeltaE", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaE (GeV);FCAL Shower #DeltaZ (cm)", 100, locDeltaEMin, locDeltaEMax, 100, locDeltaZMin, locDeltaZMax);
	dPluginHist_FCAL_DeltaZVsDeltaT = new TH2F("dPluginHist_FCAL_DeltaZVsDeltaT", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaT (ns);FCAL Shower #DeltaZ (cm)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaZMin, locDeltaZMax);
	dPluginHist_FCAL_DeltaEVsDeltaT = new TH2F("dPluginHist_FCAL_DeltaEVsDeltaT", "Simulation Reconstructed Uncertainty Study;FCAL Shower #DeltaT (ns);FCAL Shower #DeltaE (GeV)", 100, locDeltaTMin, locDeltaTMax, 100, locDeltaEMin, locDeltaEMax);

	float locShowerSigmaXMin = 0.0, locShowerSigmaXMax = 4.0;
	float locShowerSigmaYMin = 0.0, locShowerSigmaYMax = 4.0;
	float locShowerSigmaZMin = -0.3, locShowerSigmaZMax = 0.3;
	float locShowerSigmaEMin = 0.0, locShowerSigmaEMax = 0.15;
	float locShowerSigmaTMin = -0.3, locShowerSigmaTMax = 0.3;

	//ShowerSigmaX Dependence
	dPluginHist_FCAL_ShowerSigmaXVsX = new TH2F("dPluginHist_FCAL_ShowerSigmaXVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #sigma X (cm)", 500, locXMin, locXMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_FCAL_ShowerSigmaXVsY = new TH2F("dPluginHist_FCAL_ShowerSigmaXVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #sigma X (cm)", 100, locYMin, locYMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_FCAL_ShowerSigmaXVsZ = new TH2F("dPluginHist_FCAL_ShowerSigmaXVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #sigma X (cm)", 100, locZMin, locZMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_FCAL_ShowerSigmaXVsE = new TH2F("dPluginHist_FCAL_ShowerSigmaXVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #sigma X (cm)", 100, locEMin, locEMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);
	dPluginHist_FCAL_ShowerSigmaXVsT = new TH2F("dPluginHist_FCAL_ShowerSigmaXVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #sigma X (cm)", 100, locTMin, locTMax, 400, locShowerSigmaXMin, locShowerSigmaXMax);

	//ShowerSigmaY Dependence
	dPluginHist_FCAL_ShowerSigmaYVsX = new TH2F("dPluginHist_FCAL_ShowerSigmaYVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #sigma Y (cm)", 100, locXMin, locXMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_FCAL_ShowerSigmaYVsY = new TH2F("dPluginHist_FCAL_ShowerSigmaYVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #sigma Y (cm)", 500, locYMin, locYMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_FCAL_ShowerSigmaYVsZ = new TH2F("dPluginHist_FCAL_ShowerSigmaYVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #sigma Y (cm)", 100, locZMin, locZMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_FCAL_ShowerSigmaYVsE = new TH2F("dPluginHist_FCAL_ShowerSigmaYVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #sigma Y (cm)", 100, locEMin, locEMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);
	dPluginHist_FCAL_ShowerSigmaYVsT = new TH2F("dPluginHist_FCAL_ShowerSigmaYVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #sigma Y (cm)", 100, locTMin, locTMax, 400, locShowerSigmaYMin, locShowerSigmaYMax);

	//ShowerSigmaZ Dependence
	dPluginHist_FCAL_ShowerSigmaZVsX = new TH2F("dPluginHist_FCAL_ShowerSigmaZVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #sigma Z (cm)", 100, locXMin, locXMax, 100, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_FCAL_ShowerSigmaZVsY = new TH2F("dPluginHist_FCAL_ShowerSigmaZVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #sigma Z (cm)", 100, locYMin, locYMax, 100, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_FCAL_ShowerSigmaZVsZ = new TH2F("dPluginHist_FCAL_ShowerSigmaZVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #sigma Z (cm)", 100, locZMin, locZMax, 100, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_FCAL_ShowerSigmaZVsE = new TH2F("dPluginHist_FCAL_ShowerSigmaZVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #sigma Z (cm)", 100, locEMin, locEMax, 100, locShowerSigmaZMin, locShowerSigmaZMax);
	dPluginHist_FCAL_ShowerSigmaZVsT = new TH2F("dPluginHist_FCAL_ShowerSigmaZVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #sigma Z (cm)", 100, locTMin, locTMax, 100, locShowerSigmaZMin, locShowerSigmaZMax);

	//ShowerSigmaE Dependence
	dPluginHist_FCAL_ShowerSigmaEVsX = new TH2F("dPluginHist_FCAL_ShowerSigmaEVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #sigma E (GeV)", 100, locXMin, locXMax, 400, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_FCAL_ShowerSigmaEVsY = new TH2F("dPluginHist_FCAL_ShowerSigmaEVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #sigma E (GeV)", 100, locYMin, locYMax, 400, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_FCAL_ShowerSigmaEVsZ = new TH2F("dPluginHist_FCAL_ShowerSigmaEVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #sigma E (GeV)", 100, locZMin, locZMax, 400, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_FCAL_ShowerSigmaEVsE = new TH2F("dPluginHist_FCAL_ShowerSigmaEVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #sigma E (GeV)", 200, locEMin, locEMax, 400, locShowerSigmaEMin, locShowerSigmaEMax);
	dPluginHist_FCAL_ShowerSigmaEVsT = new TH2F("dPluginHist_FCAL_ShowerSigmaEVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #sigma E (GeV)", 100, locTMin, locTMax, 400, locShowerSigmaEMin, locShowerSigmaEMax);

	//ShowerSigmaT Dependence
	dPluginHist_FCAL_ShowerSigmaTVsX = new TH2F("dPluginHist_FCAL_ShowerSigmaTVsX", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured X (cm);FCAL Shower #sigma T (ns)", 100, locXMin, locXMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_FCAL_ShowerSigmaTVsY = new TH2F("dPluginHist_FCAL_ShowerSigmaTVsY", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Y (cm);FCAL Shower #sigma T (ns)", 100, locYMin, locYMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_FCAL_ShowerSigmaTVsZ = new TH2F("dPluginHist_FCAL_ShowerSigmaTVsZ", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured Z (cm);FCAL Shower #sigma T (ns)", 100, locZMin, locZMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_FCAL_ShowerSigmaTVsE = new TH2F("dPluginHist_FCAL_ShowerSigmaTVsE", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured E (GeV);FCAL Shower #sigma T (ns)", 100, locEMin, locEMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);
	dPluginHist_FCAL_ShowerSigmaTVsT = new TH2F("dPluginHist_FCAL_ShowerSigmaTVsT", "Simulation Reconstructed Uncertainty Study;FCAL Shower Measured T (ns);FCAL Shower #sigma T (ns)", 100, locTMin, locTMax, 100, locShowerSigmaTMin, locShowerSigmaTMax);

}

void FCALSelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t FCALSelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either FCALSelector::GetEntry() or TBranch::GetEntry()
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
	dPluginHist_FCAL_PathLengthCorrection->Fill(dTrueE, dPathLengthCorrection);

	bool locVsMeasuredFlag = true; //if true, swap true/meas
	if(locVsMeasuredFlag == true){
		//delta = meas - true
		dTrueX += dDeltaX;
		dTrueY += dDeltaY;
		dTrueZ += dDeltaZ;
		dTrueE += dDeltaE;
		dTrueT += dDeltaT;
	}

	//DeltaX Dependence
	dPluginHist_FCAL_DeltaXVsX->Fill(dTrueX, dDeltaX);
	dPluginHist_FCAL_DeltaXVsY->Fill(dTrueY, dDeltaX);
	dPluginHist_FCAL_DeltaXVsZ->Fill(dTrueZ, dDeltaX);
	dPluginHist_FCAL_DeltaXVsE->Fill(dTrueE, dDeltaX);
	dPluginHist_FCAL_DeltaXVsT->Fill(dTrueT, dDeltaX);

	//DeltaY Dependence
	dPluginHist_FCAL_DeltaYVsX->Fill(dTrueX, dDeltaY);
	dPluginHist_FCAL_DeltaYVsY->Fill(dTrueY, dDeltaY);
	dPluginHist_FCAL_DeltaYVsZ->Fill(dTrueZ, dDeltaY);
	dPluginHist_FCAL_DeltaYVsE->Fill(dTrueE, dDeltaY);
	dPluginHist_FCAL_DeltaYVsT->Fill(dTrueT, dDeltaY);

	//DeltaZ Dependence
	dPluginHist_FCAL_DeltaZVsX->Fill(dTrueX, dDeltaZ);
	dPluginHist_FCAL_DeltaZVsY->Fill(dTrueY, dDeltaZ);
	dPluginHist_FCAL_DeltaZVsZ->Fill(dTrueZ, dDeltaZ);
	dPluginHist_FCAL_DeltaZVsE->Fill(dTrueE, dDeltaZ);
	dPluginHist_FCAL_DeltaZVsT->Fill(dTrueT, dDeltaZ);

	//DeltaE Dependence
	dPluginHist_FCAL_DeltaEVsX->Fill(dTrueX, dDeltaE);
	dPluginHist_FCAL_DeltaEVsY->Fill(dTrueY, dDeltaE);
	dPluginHist_FCAL_DeltaEVsZ->Fill(dTrueZ, dDeltaE);
	dPluginHist_FCAL_DeltaEVsE->Fill(dTrueE, dDeltaE);
	dPluginHist_FCAL_DeltaEVsT->Fill(dTrueT, dDeltaE);

	//DeltaT Dependence
	dPluginHist_FCAL_DeltaTVsX->Fill(dTrueX, dDeltaT);
	dPluginHist_FCAL_DeltaTVsY->Fill(dTrueY, dDeltaT);
	dPluginHist_FCAL_DeltaTVsZ->Fill(dTrueZ, dDeltaT);
	dPluginHist_FCAL_DeltaTVsE->Fill(dTrueE, dDeltaT);
	dPluginHist_FCAL_DeltaTVsT->Fill(dTrueT, dDeltaT);

	//Common Dependence
	dPluginHist_FCAL_DeltaXVsDeltaY->Fill(dDeltaY, dDeltaX);
	dPluginHist_FCAL_DeltaXVsDeltaZ->Fill(dDeltaZ, dDeltaX);
	dPluginHist_FCAL_DeltaXVsDeltaE->Fill(dDeltaE, dDeltaX);
	dPluginHist_FCAL_DeltaXVsDeltaT->Fill(dDeltaT, dDeltaX);
	dPluginHist_FCAL_DeltaYVsDeltaZ->Fill(dDeltaZ, dDeltaY);
	dPluginHist_FCAL_DeltaYVsDeltaE->Fill(dDeltaE, dDeltaY);
	dPluginHist_FCAL_DeltaYVsDeltaT->Fill(dDeltaT, dDeltaY);
	dPluginHist_FCAL_DeltaZVsDeltaE->Fill(dDeltaE, dDeltaZ);
	dPluginHist_FCAL_DeltaZVsDeltaT->Fill(dDeltaT, dDeltaZ);
	dPluginHist_FCAL_DeltaEVsDeltaT->Fill(dDeltaT, dDeltaE);

	//ShowerSigmaX Dependence
	dPluginHist_FCAL_ShowerSigmaXVsX->Fill(dTrueX, dShowerUncertaintyX);
	dPluginHist_FCAL_ShowerSigmaXVsY->Fill(dTrueY, dShowerUncertaintyX);
	dPluginHist_FCAL_ShowerSigmaXVsZ->Fill(dTrueZ, dShowerUncertaintyX);
	dPluginHist_FCAL_ShowerSigmaXVsE->Fill(dTrueE, dShowerUncertaintyX);
	dPluginHist_FCAL_ShowerSigmaXVsT->Fill(dTrueT, dShowerUncertaintyX);

	//ShowerSigmaY Dependence
	dPluginHist_FCAL_ShowerSigmaYVsX->Fill(dTrueX, dShowerUncertaintyY);
	dPluginHist_FCAL_ShowerSigmaYVsY->Fill(dTrueY, dShowerUncertaintyY);
	dPluginHist_FCAL_ShowerSigmaYVsZ->Fill(dTrueZ, dShowerUncertaintyY);
	dPluginHist_FCAL_ShowerSigmaYVsE->Fill(dTrueE, dShowerUncertaintyY);
	dPluginHist_FCAL_ShowerSigmaYVsT->Fill(dTrueT, dShowerUncertaintyY);

	//ShowerSigmaZ Dependence
	dPluginHist_FCAL_ShowerSigmaZVsX->Fill(dTrueX, dShowerUncertaintyZ);
	dPluginHist_FCAL_ShowerSigmaZVsY->Fill(dTrueY, dShowerUncertaintyZ);
	dPluginHist_FCAL_ShowerSigmaZVsZ->Fill(dTrueZ, dShowerUncertaintyZ);
	dPluginHist_FCAL_ShowerSigmaZVsE->Fill(dTrueE, dShowerUncertaintyZ);
	dPluginHist_FCAL_ShowerSigmaZVsT->Fill(dTrueT, dShowerUncertaintyZ);

	//ShowerSigmaE Dependence
	dPluginHist_FCAL_ShowerSigmaEVsX->Fill(dTrueX, dShowerUncertaintyE);
	dPluginHist_FCAL_ShowerSigmaEVsY->Fill(dTrueY, dShowerUncertaintyE);
	dPluginHist_FCAL_ShowerSigmaEVsZ->Fill(dTrueZ, dShowerUncertaintyE);
	dPluginHist_FCAL_ShowerSigmaEVsE->Fill(dTrueE, dShowerUncertaintyE);
	dPluginHist_FCAL_ShowerSigmaEVsT->Fill(dTrueT, dShowerUncertaintyE);

	//ShowerSigmaT Dependence
	dPluginHist_FCAL_ShowerSigmaTVsX->Fill(dTrueX, dShowerUncertaintyT);
	dPluginHist_FCAL_ShowerSigmaTVsY->Fill(dTrueY, dShowerUncertaintyT);
	dPluginHist_FCAL_ShowerSigmaTVsZ->Fill(dTrueZ, dShowerUncertaintyT);
	dPluginHist_FCAL_ShowerSigmaTVsE->Fill(dTrueE, dShowerUncertaintyT);
	dPluginHist_FCAL_ShowerSigmaTVsT->Fill(dTrueT, dShowerUncertaintyT);

   return kTRUE;
}

void FCALSelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void FCALSelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	dOutputFile->Write();
}

