#define Selector_p3pi_cxx
// The class definition in Selector_p3pi.h has been generated automatically
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
// root> T->Process("Selector_p3pi.C")
// root> T->Process("Selector_p3pi.C","some options")
// root> T->Process("Selector_p3pi.C+")
//

#include "Selector_p3pi.h"
#include <TH2.h>
#include <TStyle.h>


void Selector_p3pi::Begin(TTree * /*tree*/)
{
   // The Begin() function is called at the start of the query.
   // When running with PROOF Begin() is only called on the client.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

	//Create ROOT File (if using PROOF, do differently)
   dFile = new TFile("p3pi_hists.root", "RECREATE");

	//Create histograms (if using PROOF, create in SlaveBegin() instead!)
	string locHistTitle = ";#gamma#gamma Invariant Mass (GeV/c^{2})";
	dHist_Pi0Mass = new TH1I("Pi0Mass", locHistTitle.c_str(), 600, 0.0, 0.3);

	locHistTitle = ";#gammap#rightarrow#pi^{#plus}#pi^{#minus}#gamma#gamma Missing Mass Squared (GeV/c^{2})^{2}";
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", locHistTitle.c_str(), 600, -0.06, 0.06);

	locHistTitle = ";#pi^{#plus}#pi^{#minus}#gamma#gamma Invariant Mass (GeV/c^{2})";
	dHist_OmegaMass = new TH1I("OmegaMass", locHistTitle.c_str(), 600, 0.5, 1.1);
}

void Selector_p3pi::SlaveBegin(TTree * /*tree*/)
{
   // The SlaveBegin() function is called after the Begin() function.
   // When running with PROOF SlaveBegin() is called on each slave server.
   // The tree argument is deprecated (on PROOF 0 is passed).

   TString option = GetOption();

}

Bool_t Selector_p3pi::Process(Long64_t entry)
{
   // The Process() function is called for each entry in the tree (or possibly
   // keyed object in the case of PROOF) to be processed. The entry argument
   // specifies which entry in the currently loaded tree is to be processed.
   // It can be passed to either Selector_p3pi::GetEntry() or TBranch::GetEntry()
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

	TLorentzVector locTargetP4(0.0, 0.0, 0.0, 0.938272046);

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < NumCombos; ++loc_i)
	{
		// Is used to mark when combos are cut
		if(IsComboCut[loc_i]) // Is false initially
			continue; // Combo has been cut previously

		// Particle info is split between combo-dependent and combo-independent
		// For combo-dependent (e.g. PID, kinfit p4), use the branches starting with the particle name (<Name>)
			// e.g. PiPlus__P4_KinFit
		// For combo-independent (e.g. measured p4, dE/dx), use the branches starting with either: 
			// "ChargedHypo," "NeutralShower," or "Beam"
			// However, these are arrays. The array index that you need is given by the branches:
			// "<Name>__ChargedIndex," "<Name>__ShowerIndex," or "ComboBeam__BeamIndex"

		// Get Measured Neutral P4's: Combo-dependent (P4 defined by vertex position)
		TLorentzVector& locPhoton1P4 = *((TLorentzVector*)Photon1__P4_Measured->At(loc_i));
		TLorentzVector& locPhoton2P4 = *((TLorentzVector*)Photon2__P4_Measured->At(loc_i));

		// Get Measured Charged P4's: Combo-independent
		Int_t locPiPlusIndex = PiPlus__ChargedIndex[loc_i];
		TLorentzVector& locPiPlusP4 = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locPiPlusIndex));

		Int_t locPiMinusIndex = PiMinus__ChargedIndex[loc_i];
		TLorentzVector& locPiMinusP4 = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locPiMinusIndex));

		Int_t locProtonIndex = Proton__ChargedIndex[loc_i];
		TLorentzVector& locProtonP4 = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locProtonIndex));

		// Get Measured Beam P4: Combo-independent
		Int_t locBeamIndex = ComboBeam__BeamIndex[loc_i];
		TLorentzVector& locBeamP4 = *((TLorentzVector*)ComboBeam__P4_KinFit->At(locBeamIndex));

		// Combine 4-vectors
		TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;
		TLorentzVector locOmegaP4 = locPiPlusP4 + locPiMinusP4 + locPi0P4;
		TLorentzVector locFinalStateP4 = locOmegaP4 + locProtonP4;
		TLorentzVector locMissingP4 = locBeamP4 + locTargetP4 - locFinalStateP4;

		//Histogram pi0 mass & cut
		double locPi0Mass = locPi0P4.M();
		dHist_Pi0Mass->Fill(locPi0Mass);
		if((locPi0Mass < 0.0775209) || (locPi0Mass > 0.188047))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		//Histogram missing mass squared & cut
		double locMissingMassSquared = locMissingP4.M2();
		dHist_MissingMassSquared->Fill(locMissingMassSquared);
		if((locMissingMassSquared < -0.01) || (locMissingMassSquared > 0.005))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		//Histogram omega mass
		double locOmegaMass = locOmegaP4.M();
		dHist_OmegaMass->Fill(locOmegaMass);
	}

   return kTRUE;
}

void Selector_p3pi::SlaveTerminate()
{
   // The SlaveTerminate() function is called after all entries or objects
   // have been processed. When running with PROOF SlaveTerminate() is called
   // on each slave server.

}

void Selector_p3pi::Terminate()
{
   // The Terminate() function is the last function to be called during
   // a query. It always runs on the client, it can be used to present
   // the results graphically or save the results to file.

	//Write and close output file (if using PROOF, do in SlaveTerminate() instead!)
	dFile->Write();
	dFile->Close();
}

