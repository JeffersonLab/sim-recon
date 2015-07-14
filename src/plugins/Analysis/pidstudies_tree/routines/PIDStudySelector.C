#define PIDStudySelector_cxx
// The class definition in PIDStudySelector.h has been generated automatically
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
// Root > T->Process("PIDStudySelector.C")
// Root > T->Process("PIDStudySelector.C","some options")
// Root > T->Process("PIDStudySelector.C+")
//

#include "PIDStudySelector.h"
#include <TH2.h>
#include <TStyle.h>

double Convert_RadiansToDegrees(double locRadians){
	return locRadians*180.0/3.141592654;
}

void PIDStudySelector::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();
	TH2F *loc2FHist;
	TH1F *loc1FHist;
	ostringstream locHistName, locHistTitle;
	unsigned int loc_i, loc_j;

	vector<int> locParticleIDVector; //p, pi+, pi-, k+
	vector<int> locParticleIDVector_Results; //p, pi+, pi-, k+, pbar, k-
	vector<string> locParticleSymbolVector, locParticleSymbolVector_Results;
	int locPID;
	string locParticleSymbol;

	locParticleIDVector.push_back(14);  locParticleIDVector.push_back(8);  locParticleIDVector.push_back(9);  locParticleIDVector.push_back(11);
	locParticleIDVector_Results.push_back(14);  locParticleIDVector_Results.push_back(8);  locParticleIDVector_Results.push_back(9);  locParticleIDVector_Results.push_back(11);
	locParticleIDVector_Results.push_back(15);  locParticleIDVector_Results.push_back(12);

	locParticleSymbolVector.push_back("p");  locParticleSymbolVector.push_back("#pi^{+}");  locParticleSymbolVector.push_back("#pi^{-}");  locParticleSymbolVector.push_back("K^{+}");
	locParticleSymbolVector_Results.push_back("p");  locParticleSymbolVector_Results.push_back("#pi^{+}");  locParticleSymbolVector_Results.push_back("#pi^{-}");  locParticleSymbolVector_Results.push_back("K^{+}");
	locParticleSymbolVector_Results.push_back("#bar{p}");  locParticleSymbolVector_Results.push_back("K^{-}");

	int locNumBinsP = 900;
	int locNumBinsTheta = 600;
	int locNumBinsConfidenceLevel = 1000;
	float locRangeMinP = 0.0, locRangeMaxP = 9.0;
	float locRangeMinTheta = 0.0, locRangeMaxTheta = 150.0;
	float locRangeMinConfidenceLevel = 0.0, locRangeMaxConfidenceLevel = 1.0;


	dOutputFile = new TFile("dh_PIDStudies.root", "RECREATE");

	for(loc_i = 0; loc_i < locParticleIDVector.size(); loc_i++){
		locPID = locParticleIDVector[loc_i];
		locParticleSymbol = locParticleSymbolVector[loc_i];

		locHistName.str("");
		locHistName << "dh_GeneratedTracksPVsTheta_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Generated " << locParticleSymbol << "'s;#theta (degrees);p (GeV/c)";
		loc2FHist = new TH2F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsTheta, locRangeMinTheta, locRangeMaxTheta, locNumBinsP, locRangeMinP, locRangeMaxP);

		locHistName.str("");
		locHistName << "dh_NonReconstructedTracksPVsTheta_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Non-Reconstructed " << locParticleSymbol << "'s;#theta (degrees);p (GeV/c)";
		loc2FHist = new TH2F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsTheta, locRangeMinTheta, locRangeMaxTheta, locNumBinsP, locRangeMinP, locRangeMaxP);

		locHistName.str("");
		locHistName << "dh_ReconstructedTracksPVsTheta_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Reconstructed " << locParticleSymbol << "'s;#theta (degrees);p (GeV/c)";
		loc2FHist = new TH2F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsTheta, locRangeMinTheta, locRangeMaxTheta, locNumBinsP, locRangeMinP, locRangeMaxP);

		locHistName.str("");
		locHistName << "dh_IdentifiedTracksPVsTheta_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Identified " << locParticleSymbol << "'s;#theta (degrees);p (GeV/c)";
		loc2FHist = new TH2F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsTheta, locRangeMinTheta, locRangeMaxTheta, locNumBinsP, locRangeMinP, locRangeMaxP);

		locHistName.str("");
		locHistName << "dh_MisIdentifiedTracksPVsTheta_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Misidentified " << locParticleSymbol << "'s;#theta (degrees);p (GeV/c)";
		loc2FHist = new TH2F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsTheta, locRangeMinTheta, locRangeMaxTheta, locNumBinsP, locRangeMinP, locRangeMaxP);

		for(loc_j = 0; loc_j < locParticleIDVector_Results.size(); loc_j++){
			if(loc_j == loc_i)
				continue;
			locHistName.str("");
			locHistName << "dh_MisIdentifiedAsTracksPVsTheta_" << locPID << "_" << locParticleIDVector_Results[loc_j];
			locHistTitle.str("");
			locHistTitle << "MisIdentified " << locParticleSymbol << "as " << locParticleSymbolVector_Results[loc_j] << ";#theta (degrees);p (GeV/c)";
			loc2FHist = new TH2F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsTheta, locRangeMinTheta, locRangeMaxTheta, locNumBinsP, locRangeMinP, locRangeMaxP);
		}

		for(loc_j = 0; loc_j < locParticleIDVector_Results.size(); loc_j++){
			locHistName.str("");
			locHistName << "dh_ConfidenceLevel_ProblemArea_Timing_" << locPID << "_" << locParticleIDVector_Results[loc_j];
			locHistTitle.str("");
			locHistTitle << locParticleSymbol << " as " << locParticleSymbolVector_Results[loc_j] << ", 1.5 < p (GeV/c) < 3.0, 30#circ < #theta < 120#circ;Timing Confidence Level";
			loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);

			locHistName.str("");
			locHistName << "dh_ConfidenceLevel_ProblemArea_DCdEdx_" << locPID << "_" << locParticleIDVector_Results[loc_j];
			locHistTitle.str("");
			locHistTitle << locParticleSymbol << " as " << locParticleSymbolVector_Results[loc_j] << ", 1.5 < p (GeV/c) < 3.0, 30#circ < #theta < 120#circ;DC dEdx Confidence Level";
			loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);
		}

		//CL hists
		locHistName.str("");
		locHistName << "dh_ConfidenceLevel_Overall_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Overall " << locParticleSymbol << " PID;Confidence Level";
		loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);

		locHistName.str("");
		locHistName << "dh_ConfidenceLevel_Tracking_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Tracking " << locParticleSymbol << " PID;Confidence Level";
		loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);

		locHistName.str("");
		locHistName << "dh_ConfidenceLevel_DCdEdx_" << locPID;
		locHistTitle.str("");
		locHistTitle << "DC dEdx " << locParticleSymbol << " PID;Confidence Level";
		loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);

		locHistName.str("");
		locHistName << "dh_ConfidenceLevel_Timing_" << locPID;
		locHistTitle.str("");
		locHistTitle << "Timing " << locParticleSymbol << " PID;Confidence Level";
		loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);

		locHistName.str("");
		locHistName << "dh_ConfidenceLevel_TOFdEdx_" << locPID;
		locHistTitle.str("");
		locHistTitle << "TOFdEdx " << locParticleSymbol << " PID;Confidence Level";
		loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);

		locHistName.str("");
		locHistName << "dh_ConfidenceLevel_BCALdEdx_" << locPID;
		locHistTitle.str("");
		locHistTitle << "BCALdEdx " << locParticleSymbol << " PID;Confidence Level";
		loc1FHist = new TH1F(locHistName.str().c_str(), locHistTitle.str().c_str(), locNumBinsConfidenceLevel, locRangeMinConfidenceLevel, locRangeMaxConfidenceLevel);
	}
}

void PIDStudySelector::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t PIDStudySelector::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either PIDStudySelector::GetEntry() or TBranch::GetEntry()
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
		cout << "entry " << (entry + 1) << " of " << fChain->GetEntries() << endl;

	unsigned int loc_i, loc_j;
	const MCReconstructionStatus *locMCReconstructionStatus;
	const ReconstructedHypothesis *locReconstructedHypothesis;
	float locFOM, locThrownP, locThrownTheta_Degrees;
	Particle_t locThrownID, locPID;
	DLorentzVector locThrownFourMomentum, locThrownSpacetimeVertex, locFourMomentum, locSpacetimeVertex;
	TH2F *loc2FHist;
	TH1F *loc1FHist;
	ostringstream locHistName;

	for(loc_i = 0; loc_i < dMCReconstructionStatusVector.size(); loc_i++){
		locMCReconstructionStatus = dMCReconstructionStatusVector[loc_i];
		locThrownID = locMCReconstructionStatus->dThrownID;

		locThrownFourMomentum = locMCReconstructionStatus->dThrownFourMomentum;
		locThrownSpacetimeVertex = locMCReconstructionStatus->dThrownSpacetimeVertex;

		locThrownP = locThrownFourMomentum.P();
		locThrownTheta_Degrees = Convert_RadiansToDegrees(locThrownFourMomentum.Theta());

		locHistName.str("");
		locHistName << "dh_GeneratedTracksPVsTheta_" << int(locThrownID);
		loc2FHist = (TH2F*)gDirectory->Get(locHistName.str().c_str());
		loc2FHist->Fill(locThrownTheta_Degrees, locThrownP);

		if(locMCReconstructionStatus->dReconstructedHypothesisVector.size() == 0){
			locHistName.str("");
			locHistName << "dh_NonReconstructedTracksPVsTheta_" << int(locThrownID);
			loc2FHist = (TH2F*)gDirectory->Get(locHistName.str().c_str());
			loc2FHist->Fill(locThrownTheta_Degrees, locThrownP);
		}

		for(loc_j = 0; loc_j < locMCReconstructionStatus->dReconstructedHypothesisVector.size(); loc_j++){
			locReconstructedHypothesis = locMCReconstructionStatus->dReconstructedHypothesisVector[loc_j];
			locPID = locReconstructedHypothesis->dPID;
			locFourMomentum = locReconstructedHypothesis->dFourMomentum;
			locSpacetimeVertex = locReconstructedHypothesis->dSpacetimeVertex;

			if(loc_j == 0){
				locHistName.str("");
				locHistName << "dh_ReconstructedTracksPVsTheta_" << int(locThrownID);
				loc2FHist = (TH2F*)gDirectory->Get(locHistName.str().c_str());
				loc2FHist->Fill(locThrownTheta_Degrees, locThrownP);

				if(locThrownID == locPID){
					locHistName.str("");
					locHistName << "dh_IdentifiedTracksPVsTheta_" << int(locThrownID);
					loc2FHist = (TH2F*)gDirectory->Get(locHistName.str().c_str());
					loc2FHist->Fill(locThrownTheta_Degrees, locThrownP);
				}else{
					locHistName.str("");
					locHistName << "dh_MisIdentifiedTracksPVsTheta_" << int(locThrownID);
					loc2FHist = (TH2F*)gDirectory->Get(locHistName.str().c_str());
					loc2FHist->Fill(locThrownTheta_Degrees, locThrownP);

					locHistName.str("");
					locHistName << "dh_MisIdentifiedAsTracksPVsTheta_" << int(locThrownID) << "_" << int(locPID);
					loc2FHist = (TH2F*)gDirectory->Get(locHistName.str().c_str());
					if(loc2FHist != NULL)
						loc2FHist->Fill(locThrownTheta_Degrees, locThrownP);
				}
			}


			if(locThrownID == locPID){
				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_Overall, locReconstructedHypothesis->dNDF_Overall);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_Overall_" << int(locThrownID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				loc1FHist->Fill(locFOM);

				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_Tracking, locReconstructedHypothesis->dNDF_Tracking);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_Tracking_" << int(locThrownID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				loc1FHist->Fill(locFOM);

				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_DCdEdx, locReconstructedHypothesis->dNDF_DCdEdx);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_DCdEdx_" << int(locThrownID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				loc1FHist->Fill(locFOM);

				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_Timing, locReconstructedHypothesis->dNDF_Timing);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_Timing_" << int(locThrownID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				loc1FHist->Fill(locFOM);
/*
				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_TOFdEdx, locReconstructedHypothesis->dNDF_TOFdEdx);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_TOFdEdx_" << int(locThrownID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				loc1FHist->Fill(locFOM);

				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_BCALdEdx, locReconstructedHypothesis->dNDF_BCALdEdx);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_BCALdEdx_" << int(locThrownID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				loc1FHist->Fill(locFOM);
*/
			}

			//problem area
			if((int(locThrownID) == 8) && (locThrownP > 1.5) && (locThrownP < 3.0) && (locThrownTheta_Degrees > 30.0) && (locThrownTheta_Degrees < 120.0)){
				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_Timing, locReconstructedHypothesis->dNDF_Timing);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_ProblemArea_Timing_" << int(locThrownID) << "_" << int(locPID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				if(loc1FHist != NULL)
					loc1FHist->Fill(locFOM);

				locFOM = TMath::Prob(locReconstructedHypothesis->dChiSq_DCdEdx, locReconstructedHypothesis->dNDF_DCdEdx);
				locHistName.str("");
				locHistName << "dh_ConfidenceLevel_ProblemArea_DCdEdx_" << int(locThrownID) << "_" << int(locPID);
				loc1FHist = (TH1F*)gDirectory->Get(locHistName.str().c_str());
				if(loc1FHist != NULL)
					loc1FHist->Fill(locFOM);
			}

		}

	}


	//fraction of generated tracks reconstructed of each type as function of p & theta
	

	//fraction of generated tracks id'd of each type, as each type, as function of p & theta

	//fraction of reconstructed tracks id'd of each type, as each type, as function of p & theta

	//CL of each PID source for each type in bins of p & theta //& % of events failing 1% CL cut

   return kTRUE;
}

void PIDStudySelector::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void PIDStudySelector::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	dOutputFile->Write();
}

