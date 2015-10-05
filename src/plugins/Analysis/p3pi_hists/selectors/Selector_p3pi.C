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

	//Target p4 //don't need to hard-code, info is available in TTree::GetUserInfo()
	dTargetP4.SetPxPyPzE(0.0, 0.0, 0.0, 0.938272046);

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

	/********************************************* SETUP UNIQUENESS TRACKING ********************************************/

	//PREVENT-DOUBLE COUNTING WHEN HISTOGRAMMING
		//Sometimes, some content is the exact same between one combo and the next
			//e.g. maybe two combos have different beam particles, but the same data for the final-state
		//When histogramming, you don't want to double-count when this happens: artificially inflates your signal (or background)
		//So, for each quantity you histogram, keep track of what particles you used (for a given combo)
			//Use the combo-independent particle indices (i.e. the indices to "ChargedHypo," "NeutralShower," and/or "Beam"
		//Then for each combo, just compare to what you used before, and make sure it's unique
		//Tip: Use std::set if possible, it's easier/faster to search

	//BE CAREFUL: This is very tricky, especially if you have a reaction with the same PID appearing in different steps
		//e.g. g p -> pi+ pi- omega p,   omega -> pi+ pi- pi0
		//Here, beware combos that are identical except for swapped pions (e.g. pi+_1 <--> pi+_2)

	//Pi0
	set<set<Int_t> > locParticlesUsedSoFar_Pi0Mass; //1st dimension: Combo. 2nd: particles used

	//Missing Mass
	//since missing mass uses charged, neutral and beam, need to use multiple containers
	vector<set<Int_t> > locParticlesUsedSoFar_MissingMass_Charged; //1st dimension: Combo. 2nd: particles used
	vector<set<Int_t> > locParticlesUsedSoFar_MissingMass_Showers; //1st dimension: Combo. 2nd: particles used
	vector<Int_t> locParticlesUsedSoFar_MissingMass_Beam; //1st dimension: Combo

	//Omega
	//since omega uses both charged and neutral, need to use multiple containers
	vector<set<Int_t> > locParticlesUsedSoFar_OmegaMass_Charged; //1st dimension: Combo. 2nd: particles used
	vector<set<Int_t> > locParticlesUsedSoFar_OmegaMass_Showers; //1st dimension: Combo. 2nd: particles used

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	for(UInt_t loc_i = 0; loc_i < NumCombos; ++loc_i)
	{
		// Is used to mark when combos are cut
		if(IsComboCut[loc_i]) // Is false initially
			continue; // Combo has been cut previously

		/***************************************** READ/SETUP DATA FOR THIS COMBO ****************************************/

		// Particle info is split between combo-dependent and combo-independent
		// For combo-dependent (e.g. PID, kinfit p4), use the branches starting with the particle name (<Name>)
			// e.g. PiPlus__P4_KinFit
		// For combo-independent (e.g. measured p4, dE/dx), use the branches starting with either: 
			// "ChargedHypo," "NeutralShower," or "Beam"
			// However, these are arrays. The array index that you need is given by the branches:
			// "<Name>__ChargedIndex," "<Name>__ShowerIndex," or "ComboBeam__BeamIndex"

		// Get particle indices: These point from combo-particle to combo-independent data
		Int_t locPhoton1Index = Photon1__ShowerIndex[loc_i];
		Int_t locPhoton2Index = Photon2__ShowerIndex[loc_i];
		Int_t locPiPlusIndex = PiPlus__ChargedIndex[loc_i];
		Int_t locPiMinusIndex = PiMinus__ChargedIndex[loc_i];
		Int_t locProtonIndex = Proton__ChargedIndex[loc_i];
		Int_t locBeamIndex = ComboBeam__BeamIndex[loc_i];

		// Get Measured Neutral P4's: Combo-dependent (P4 defined by combo-dependent vertex position)
		TLorentzVector& locPhoton1P4 = *((TLorentzVector*)Photon1__P4_Measured->At(loc_i));
		TLorentzVector& locPhoton2P4 = *((TLorentzVector*)Photon2__P4_Measured->At(loc_i));

		// Get Measured Charged P4's: Combo-independent
		TLorentzVector& locPiPlusP4 = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locPiPlusIndex));
		TLorentzVector& locPiMinusP4 = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locPiMinusIndex));
		TLorentzVector& locProtonP4 = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locProtonIndex));

		// Get Measured Beam P4: Combo-independent
		TLorentzVector& locBeamP4 = *((TLorentzVector*)ComboBeam__P4_KinFit->At(locBeamIndex));

		// Combine 4-vectors
		TLorentzVector locPi0P4 = locPhoton1P4 + locPhoton2P4;
		TLorentzVector locOmegaP4 = locPiPlusP4 + locPiMinusP4 + locPi0P4;
		TLorentzVector locFinalStateP4 = locOmegaP4 + locProtonP4;
		TLorentzVector locMissingP4 = locBeamP4 + dTargetP4 - locFinalStateP4;

		/****************************************************** PI0 ******************************************************/

		//Mass
		double locPi0Mass = locPi0P4.M();

		//Build the set of photons used for the pi0 mass
		//use set in case particles are in opposite order between combos
			// that can't happen in this case, but could in other cases, so safer to just do it all the time
		set<Int_t> locParticlesUsed_Pi0Mass;
		locParticlesUsed_Pi0Mass.insert(locPhoton1Index);
		locParticlesUsed_Pi0Mass.insert(locPhoton2Index);

		//compare to what's been used so far
		if(locParticlesUsedSoFar_Pi0Mass.find(locParticlesUsed_Pi0Mass) == locParticlesUsedSoFar_Pi0Mass.end())
		{
			//unique pi0 combo: histogram it, and register this combo of particles
			dHist_Pi0Mass->Fill(locPi0Mass);
			locParticlesUsedSoFar_Pi0Mass.insert(locParticlesUsed_Pi0Mass);
		}

		//Cut pi0 mass (+/- 3 sigma)
		if((locPi0Mass < 0.0775209) || (locPi0Mass > 0.188047))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		/********************************************* MISSING MASS SQUARED **********************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4.M2();

		//Build the set of final state charged particles used for the missing mass squared
			//Will also use the pi0 set created earlier
		set<Int_t> locParticlesUsed_MissingMass_Charged;
		locParticlesUsed_MissingMass_Charged.insert(locPiPlusIndex);
		locParticlesUsed_MissingMass_Charged.insert(locPiMinusIndex);
		locParticlesUsed_MissingMass_Charged.insert(locProtonIndex);

		//compare to what's been used so far
		bool locUniqueComboFlag = true; //will mark false if otherwise
		for(size_t loc_j = 0; loc_j < locParticlesUsedSoFar_MissingMass_Charged.size(); ++loc_j)
		{
			if(locParticlesUsed_MissingMass_Charged != locParticlesUsedSoFar_MissingMass_Charged[loc_j])
				continue; //doesn't match this combo: unique so far

			//charged matches this combo, check neutrals
			if(locParticlesUsed_Pi0Mass != locParticlesUsedSoFar_MissingMass_Showers[loc_j])
				continue; //doesn't match this combo: unique so far

			//neutrals matches this combo, check beam
			if(locBeamIndex != locParticlesUsedSoFar_MissingMass_Beam[loc_j])
				continue; //doesn't match this combo: unique so far

			//matches this combo, not unique
			locUniqueComboFlag = false;
			break; //no need to keep searching
		}

		if(locUniqueComboFlag)
		{
			//histogram it
			dHist_MissingMassSquared->Fill(locMissingMassSquared);

			//register this combo of particles
			locParticlesUsedSoFar_MissingMass_Charged.push_back(locParticlesUsed_MissingMass_Charged);
			locParticlesUsedSoFar_MissingMass_Showers.push_back(locParticlesUsed_Pi0Mass);
			locParticlesUsedSoFar_MissingMass_Beam.push_back(locBeamIndex);
		}

		//Cut
		if((locMissingMassSquared < -0.01) || (locMissingMassSquared > 0.005))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		/***************************************************** OMEGA *****************************************************/

		//Mass
		double locOmegaMass = locOmegaP4.M();

		//Build the set of charged particles used for the omega mass
			//Will also use the pi0 set created earlier
		set<Int_t> locParticlesUsed_OmegaMass_Charged;
		locParticlesUsed_OmegaMass_Charged.insert(locPiPlusIndex);
		locParticlesUsed_OmegaMass_Charged.insert(locPiMinusIndex);

		//compare to what's been used so far
		locUniqueComboFlag = true; //will mark false if otherwise
		for(size_t loc_j = 0; loc_j < locParticlesUsedSoFar_OmegaMass_Charged.size(); ++loc_j)
		{
			if(locParticlesUsed_OmegaMass_Charged != locParticlesUsedSoFar_OmegaMass_Charged[loc_j])
				continue; //doesn't match this combo: unique so far

			//charged matches this combo, check neutrals
			if(locParticlesUsed_Pi0Mass != locParticlesUsedSoFar_OmegaMass_Showers[loc_j])
				continue; //doesn't match this combo: unique so far

			//matches this combo, not unique
			locUniqueComboFlag = false;
			break; //no need to keep searching
		}

		if(locUniqueComboFlag)
		{
			//histogram it
			dHist_OmegaMass->Fill(locOmegaMass);

			//register this combo of particles
			locParticlesUsedSoFar_OmegaMass_Charged.push_back(locParticlesUsed_OmegaMass_Charged);
			locParticlesUsedSoFar_OmegaMass_Showers.push_back(locParticlesUsed_Pi0Mass);
		}
	} //end combo loop

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

