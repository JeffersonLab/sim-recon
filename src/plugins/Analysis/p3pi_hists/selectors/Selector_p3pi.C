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
	dTargetP4.SetPxPyPzE(0.0, 0.0, 0.0, ParticleMass(Proton));

	//Create ROOT File (if using PROOF, do differently)
   dFile = new TFile("p3pi_hists.root", "RECREATE");

	//Create histograms (if using PROOF, create in SlaveBegin() instead!)

	string locHistTitle = ";#gamma#gamma Invariant Mass (GeV/c^{2})";
	dHist_Pi0Mass_Measured = new TH1I("Pi0Mass_Measured", locHistTitle.c_str(), 600, 0.0, 0.3);

	locHistTitle = ";#gamma#gamma Invariant Mass (GeV/c^{2})";
	dHist_Pi0Mass_KinFit = new TH1I("Pi0Mass_KinFit", locHistTitle.c_str(), 600, 0.0, 0.3);

	locHistTitle = ";#gammap#rightarrowp#pi^{#plus}#pi^{#minus}#gamma#gamma Missing Mass Squared (GeV/c^{2})^{2}";
	dHist_MissingMassSquared = new TH1I("MissingMassSquared", locHistTitle.c_str(), 600, -0.06, 0.06);

	locHistTitle = ";#pi^{#plus}#pi^{#minus}#gamma#gamma Invariant Mass (GeV/c^{2})";
	dHist_OmegaMass_Measured = new TH1I("OmegaMass_Measured", locHistTitle.c_str(), 600, 0.5, 1.1);

	locHistTitle = ";#pi^{#plus}#pi^{#minus}#gamma#gamma Invariant Mass (GeV/c^{2})";
	dHist_OmegaMass_KinFit = new TH1I("OmegaMass_KinFit", locHistTitle.c_str(), 600, 0.5, 1.1);
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

	//BE CAREFUL: This is very tricky, especially if you have a reaction with the same PID appearing in different steps
		//e.g. g p -> pi+ pi- omega p,   omega -> pi+ pi- pi0
		//Here, beware combos that are identical except for swapped pions (e.g. pi+_1 <--> pi+_2)

	//BE CAREFUL: If using charged PIDs for which hypos are not created by default (e.g. e+, e-), you must be extra cautious
		//e.g. g p -> e+ e- pi+ pi- (Yeah, I need a better example)
		//In this case, the e+(-) and pi+(-) will come from the same set of charged hypos
		//So you need to check uniqueness based on PID as well

	//In general: Could multiple particles with the same PID: Use a set (easier, faster to search)
	//In general: Multiple PIDs, so multiple sets: Contain within a map
	//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Pi0Mass; 
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaMass;

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
		// If using charged PIDs for which hypos are not created by default (e.g. e+, e-), beware!
			// The energy in the "P4_Measured" will be computed with a different mass than the one you're using
			// So you'll need to recompute it yourself.  However, the "P4_KinFit" will be fine. 

		// Get particle indices: These point from combo-particle to combo-independent data
		Int_t locPhoton1Index = Photon1__ShowerIndex[loc_i];
		Int_t locPhoton2Index = Photon2__ShowerIndex[loc_i];
		Int_t locPiPlusIndex = PiPlus__ChargedIndex[loc_i];
		Int_t locPiMinusIndex = PiMinus__ChargedIndex[loc_i];
		Int_t locProtonIndex = Proton__ChargedIndex[loc_i];
		Int_t locBeamIndex = ComboBeam__BeamIndex[loc_i];

		// Get Measured Neutral P4's: Combo-dependent (P4 defined by combo-dependent vertex position)
		TLorentzVector& locPhoton1P4_Measured = *((TLorentzVector*)Photon1__P4_Measured->At(loc_i));
		TLorentzVector& locPhoton2P4_Measured = *((TLorentzVector*)Photon2__P4_Measured->At(loc_i));

		// Get KinFit Neutral P4's: Combo-dependent
		TLorentzVector& locPhoton1P4_KinFit = *((TLorentzVector*)Photon1__P4_KinFit->At(loc_i));
		TLorentzVector& locPhoton2P4_KinFit = *((TLorentzVector*)Photon2__P4_KinFit->At(loc_i));

		// Get Measured Charged P4's: Combo-independent
		TLorentzVector& locPiPlusP4_Measured = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locPiPlusIndex));
		TLorentzVector& locPiMinusP4_Measured = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locPiMinusIndex));
		TLorentzVector& locProtonP4_Measured = *((TLorentzVector*)ChargedHypo__P4_Measured->At(locProtonIndex));

		// Get KinFit Charged P4's: Combo-dependent
		TLorentzVector& locPiPlusP4_KinFit = *((TLorentzVector*)PiPlus__P4_KinFit->At(loc_i));
		TLorentzVector& locPiMinusP4_KinFit = *((TLorentzVector*)PiMinus__P4_KinFit->At(loc_i));
		TLorentzVector& locProtonP4_KinFit = *((TLorentzVector*)Proton__P4_KinFit->At(loc_i));

		// Get Measured Beam P4: Combo-independent
		TLorentzVector& locBeamP4_Measured = *((TLorentzVector*)Beam__P4_Measured->At(locBeamIndex));

		// Get KinFit Beam P4: Combo-dependent
		TLorentzVector& locBeamP4_KinFit = *((TLorentzVector*)ComboBeam__P4_KinFit->At(locBeamIndex));

		// Combine 4-vectors
		TLorentzVector locPi0P4_Measured = locPhoton1P4_Measured + locPhoton2P4_Measured;
		TLorentzVector locPi0P4_KinFit = locPhoton1P4_KinFit + locPhoton2P4_KinFit;

		TLorentzVector locOmegaP4_Measured = locPiPlusP4_Measured + locPiMinusP4_Measured + locPi0P4_Measured;
		TLorentzVector locOmegaP4_KinFit = locPiPlusP4_KinFit + locPiMinusP4_KinFit + locPi0P4_KinFit;

		TLorentzVector locFinalStateP4_Measured = locOmegaP4_Measured + locProtonP4_Measured;
		TLorentzVector locMissingP4_Measured = locBeamP4_Measured + dTargetP4 - locFinalStateP4_Measured;

		/****************************************************** PI0 ******************************************************/

		//Mass
		double locPi0Mass_Measured = locPi0P4_Measured.M();
		double locPi0Mass_KinFit = locPi0P4_KinFit.M();

		//Build the map of particles used for the pi0 mass
		map<Particle_t, set<Int_t> > locUsedThisCombo_Pi0Mass;
		locUsedThisCombo_Pi0Mass[Gamma].insert(locPhoton1Index);
		locUsedThisCombo_Pi0Mass[Gamma].insert(locPhoton2Index);

		//compare to what's been used so far
		if(locUsedSoFar_Pi0Mass.find(locUsedThisCombo_Pi0Mass) == locUsedSoFar_Pi0Mass.end())
		{
			//unique pi0 combo: histogram it, and register this combo of particles
			dHist_Pi0Mass_Measured->Fill(locPi0Mass_Measured);
			dHist_Pi0Mass_KinFit->Fill(locPi0Mass_KinFit);
			locUsedSoFar_Pi0Mass.insert(locUsedThisCombo_Pi0Mass);
		}

		//Cut pi0 mass (+/- 3 sigma)
		if((locPi0Mass_KinFit < 0.0775209) || (locPi0Mass_KinFit > 0.188047))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		/********************************************* MISSING MASS SQUARED **********************************************/

		//Missing Mass Squared
		double locMissingMassSquared = locMissingP4_Measured.M2();

		//Build the map of particles used for the missing mass
			//For beam: Don't want to group with final-state photons. Instead use "Unknown" PID (not ideal, but it's easy). 
		map<Particle_t, set<Int_t> > locUsedThisCombo_MissingMass;
		locUsedThisCombo_MissingMass[Gamma] = locUsedThisCombo_Pi0Mass[Gamma];
		locUsedThisCombo_MissingMass[PiPlus].insert(locPiPlusIndex);
		locUsedThisCombo_MissingMass[PiMinus].insert(locPiMinusIndex);
		locUsedThisCombo_MissingMass[Proton].insert(locProtonIndex);
		locUsedThisCombo_MissingMass[Unknown].insert(locBeamIndex); //beam

		//compare to what's been used so far
		if(locUsedSoFar_MissingMass.find(locUsedThisCombo_MissingMass) == locUsedSoFar_MissingMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_MissingMassSquared->Fill(locMissingMassSquared);
			locUsedSoFar_MissingMass.insert(locUsedThisCombo_MissingMass);
		}

		//Cut
		if((locMissingMassSquared < -0.01) || (locMissingMassSquared > 0.005))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		/***************************************************** OMEGA *****************************************************/

		//Mass
		double locOmegaMass_Measured = locOmegaP4_Measured.M();
		double locOmegaMass_KinFit = locOmegaP4_KinFit.M();

		//Build the map of particles used for the omega mass
		map<Particle_t, set<Int_t> > locUsedThisCombo_OmegaMass;
		locUsedThisCombo_OmegaMass[Gamma] = locUsedThisCombo_Pi0Mass[Gamma];
		locUsedThisCombo_OmegaMass[PiPlus].insert(locPiPlusIndex);
		locUsedThisCombo_OmegaMass[PiMinus].insert(locPiMinusIndex);

		//compare to what's been used so far
		if(locUsedSoFar_OmegaMass.find(locUsedThisCombo_OmegaMass) == locUsedSoFar_OmegaMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			dHist_OmegaMass_Measured->Fill(locOmegaMass_Measured);
			dHist_OmegaMass_KinFit->Fill(locOmegaMass_KinFit);
			locUsedSoFar_OmegaMass.insert(locUsedThisCombo_OmegaMass);
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

