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

	//Mass peaks
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

	//Beam energy
	dHist_BeamEnergy = new TH1I("BeamEnergy", ";Beam Energy (GeV)", 600, 0.0, 6.0);
	dHist_MandelstamT = new TH1I("MandelstamT", ";t (GeV^{2})", 200, -2.0, 0.0);

	//Extra particles
	dHist_NumExtraTracks = new TH1I("NumExtraTracks", ";# Extra Tracks", 5, -0.5, 4.5);
	locHistTitle = ";Extra #gamma#gamma Invariant Mass (GeV/c^{2})";
	dHist_ExtraPi0InvariantMass = new TH1I("ExtraPi0InvariantMass", locHistTitle.c_str(), 300, 0.0, 0.3);

	//Angles
	dHist_OmegaPsi = new TH1I("OmegaPsi", ";#psi#circ", 180, -180.0, 180.0);
	dHist_OmegaCosTheta = new TH1I("OmegaCosTheta", ";cos(#theta)", 100, -1.0, 1.0);
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

	//In general: Could have multiple particles with the same PID: Use a set (easier, faster to search)
	//In general: Multiple PIDs, so multiple sets: Contain within a map
	//Multiple combos: Contain maps within a set (easier, faster to search)
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_Pi0Mass; 
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_MissingMass;
	set<map<Particle_t, set<Int_t> > > locUsedSoFar_OmegaMass;
	set<Int_t> locUsedSoFar_BeamEnergy;
	set<Int_t> locUsedSoFar_MandelstamT;
	set<set<Int_t> > locUsedSoFar_ExtraPi0;
	bool locNumExtraTracksFilledFlag = false;

	/************************************************* LOOP OVER COMBOS *************************************************/

	//Loop over combos
	int locNumSurvivingCombos = 0;
	for(UInt_t loc_i = 0; loc_i < NumCombos; ++loc_i)
	{
		// Is used to mark when combos are cut
		if(IsComboCut[loc_i]) // Is false initially
			continue; // Combo has been cut previously

		/***************************************** READ/SETUP DATA FOR THIS COMBO ****************************************/

		// Particle info is split between combo-dependent and combo-independent
		// For combo-dependent (e.g. PID, kinfit p4), use the branches starting with the particle name (<Name>)
			// e.g. PiPlus__P4_KinFit
		// For combo-independent (e.g. measured p4, dE/dx, etc.), use the branches starting with either: 
			// "ChargedHypo," "NeutralShower," or "Beam"
			// However, these are arrays. The array index that you need is given by the branches:
			// "<Name>__ChargedIndex," "<Name>__ShowerIndex," or "ComboBeam__BeamIndex"
		// If using charged PIDs for which hypos are not created by default (e.g. e+, e-), beware! (pi+/-, k+/-, and p are fine)
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
		TLorentzVector locFinalStateP4_KinFit = locOmegaP4_KinFit + locProtonP4_KinFit;

		TLorentzVector locInitialStateP4_Measured = locBeamP4_Measured + dTargetP4;
		TLorentzVector locInitialStateP4_KinFit = locBeamP4_KinFit + dTargetP4;

		TLorentzVector locMissingP4_Measured = locInitialStateP4_Measured - locFinalStateP4_Measured;

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
//		if((locPi0Mass_Measured < 0.0775209) || (locPi0Mass_Measured > 0.188047))
//			continue; //could also mark combo as cut, then save cut results to a new TTree
		if((locPi0Mass_KinFit < 0.102284) || (locPi0Mass_KinFit > 0.167278))
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
		if((locMissingMassSquared < -0.007) || (locMissingMassSquared > 0.005))
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
		bool locOmegaUniqueFlag = false; //will check later for asymmetry
		if(locUsedSoFar_OmegaMass.find(locUsedThisCombo_OmegaMass) == locUsedSoFar_OmegaMass.end())
		{
			//unique missing mass combo: histogram it, and register this combo of particles
			locOmegaUniqueFlag = true;
			dHist_OmegaMass_Measured->Fill(locOmegaMass_Measured);
			dHist_OmegaMass_KinFit->Fill(locOmegaMass_KinFit);
			locUsedSoFar_OmegaMass.insert(locUsedThisCombo_OmegaMass);
		}

		//Cut
		if((locOmegaMass_KinFit < 0.72) || (locOmegaMass_KinFit > 0.84))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		/************************************************ BEAM ENERGY, T *************************************************/

		//Histogram beam energy (if haven't already)
		if(locUsedSoFar_BeamEnergy.find(locBeamIndex) == locUsedSoFar_BeamEnergy.end())
		{
			dHist_BeamEnergy->Fill(locBeamP4_KinFit.E());
			locUsedSoFar_BeamEnergy.insert(locBeamIndex);
		}
		if((locBeamP4_KinFit.E() < 2.5) || (locBeamP4_KinFit.E() > 3.0))
			continue; //could also mark combo as cut, then save cut results to a new TTree

		double locT = (locProtonP4_KinFit - dTargetP4).Mag2();
		if(locUsedSoFar_MandelstamT.find(locProtonIndex) == locUsedSoFar_MandelstamT.end())
		{
			dHist_MandelstamT->Fill(locT);
			locUsedSoFar_MandelstamT.insert(locProtonIndex);
		}

		if((fabs(locT) < 0.02) || (fabs(locT) > 0.3))
			continue;

		/*********************************************** EXTRA PARTICLES *************************************************/

		//loop through all hypotheses, find how many physical tracks there are
		set<Int_t> locFoundTrackIDs; //keep track of unique track ids
		for(size_t loc_j = 0; loc_j < NumChargedHypos; ++loc_j)
			locFoundTrackIDs.insert(ChargedHypo__TrackID[loc_j]);

		int locNumExtraTracks = locFoundTrackIDs.size() - 3; //proton, pi+, pi- used
		//Fill the histogram if it hasn't already been filled for this event (quantity is combo-independent)
		if(!locNumExtraTracksFilledFlag)
		{
			dHist_NumExtraTracks->Fill(locNumExtraTracks);
			locNumExtraTracksFilledFlag = true;
		}

		//cut, requiring no extra tracks in the event
		if(locNumExtraTracks > 0)
			continue;

		//Loop through showers, see if there are any additional pi0s
		TVector3 locProductionVertex = ((TLorentzVector*)ComboBeam__X4_Measured->At(loc_i))->Vect();
		for(Int_t loc_j = 0; loc_j < Int_t(NumNeutralShowers); ++loc_j)
		{
			if((loc_j == locPhoton1Index) || (loc_j == locPhoton2Index))
				continue; //don't choose a shower that's already in the combo

			//Construct extra photon 1 p4
			TLorentzVector& locShower1X4 = *((TLorentzVector*)NeutralShower__X4_Shower->At(loc_j));
			double locShower1Energy = (NeutralShower__Energy_BCAL[loc_j] > 0.0) ? NeutralShower__Energy_BCAL[loc_j] : NeutralShower__Energy_FCAL[loc_j];
			TVector3 locExtraPhoton1P3 = locShower1X4.Vect() - locProductionVertex;
			locExtraPhoton1P3.SetMag(locShower1Energy);
			TLorentzVector locExtraPhoton1P4(locExtraPhoton1P3, locShower1Energy);

			//need 2 photons for a pi0
			for(Int_t loc_k = loc_j + 1; loc_k < Int_t(NumNeutralShowers); ++loc_k)
			{
				if((loc_k == locPhoton1Index) || (loc_k == locPhoton2Index))
					continue; //don't choose a shower that's already in the combo

				//Construct extra photon 2 p4
				TLorentzVector& locShower2X4 = *((TLorentzVector*)NeutralShower__X4_Shower->At(loc_k));
				double locShower2Energy = (NeutralShower__Energy_BCAL[loc_k] > 0.0) ? NeutralShower__Energy_BCAL[loc_k] : NeutralShower__Energy_FCAL[loc_k];
				TVector3 locExtraPhoton2P3 = locShower2X4.Vect() - locProductionVertex;
				locExtraPhoton2P3.SetMag(locShower2Energy);
				TLorentzVector locExtraPhoton2P4(locExtraPhoton2P3, locShower2Energy);

				TLorentzVector locExtraPi0P4 = locExtraPhoton1P4 + locExtraPhoton2P4;
				double locExtraPi0Mass = locExtraPi0P4.M();

				//see if this combo has been histogrammed yet
				set<Int_t> locUsedThisCombo_ExtraPi0;
				locUsedThisCombo_ExtraPi0.insert(loc_j);
				locUsedThisCombo_ExtraPi0.insert(loc_k);
				if(locUsedSoFar_ExtraPi0.find(locUsedThisCombo_ExtraPi0) == locUsedSoFar_ExtraPi0.end())
				{
					//it has not: histogram and register
					dHist_ExtraPi0InvariantMass->Fill(locExtraPi0Mass);
					locUsedSoFar_ExtraPi0.insert(locUsedThisCombo_ExtraPi0);
				}
			}
		}

		++locNumSurvivingCombos;
		if(locNumSurvivingCombos > 1)
			cout << "# combos = " << locNumSurvivingCombos << endl;

		/************************************************ OMEGA ASYMMETRY ************************************************/

		//Polarization plane:
			//Beam is in the lab z-direction
			//(FYI) Circularly polarized photon beam: Polarization rotates through the plane perpendicular to the direction of the photon: The XY Plane
			//The polarization vector is perpendicular to the direction of the photon
			//Linearly polarized photon beam: Polarization is confined to a plane along the direction of the photon
				//Plane defined by z-direction & some angle phi. Thus, polarization vector defined by phi.
				//PARA: Polarization plane parallel to the floor: The XZ plane. Polarization Vector = +/- x-axis
				//PERP: Polarization plane perpendicular to the floor: The YZ plane. Polarization Vector = +/- y-axis
			//Here: Assume that the beam polarization plane is parallel to the floor (Run 3185) (xz plane, choose +x-axis)

		//Production CM frame: The center-of-mass frame of the production step. Here: g, p -> omega, p
			//In general, the beam energy is measured more accurately than the combination of all of the final-state particles
			//So define the production CM frame using the initial state
		TVector3 locBoostVector_ProdCM = -1.0*(locInitialStateP4_KinFit.BoostVector()); //negative due to coordinate system convention

		//boost beam & proton to production CM frame
		TLorentzVector locBeamP4_ProdCM(locBeamP4_KinFit);
		locBeamP4_ProdCM.Boost(locBoostVector_ProdCM);
		TLorentzVector locProtonP4_ProdCM(locProtonP4_KinFit);
		locProtonP4_ProdCM.Boost(locBoostVector_ProdCM);

		//Production plane:
			//The production plane is the plane containing the produced particles. Here: Defined by the proton and the omega
			//However, when you boost to the production CM frame, the production plane is no longer well defined: the particles are back-to-back
			//So, by convention, define the production plane in the production CM frame by the beam and the vector meson.

		//Production CM frame axes: "HELICITY SYSTEM"
			//The z-axis is defined as the direction of the meson (omega): z = Omega
			//The y-axis is defined by the vector cross product: y = Beam X Omega
			//The x-axis is defined by the vector cross product: x = y cross z
			//However, the proton momentum is in general better known than the omega momentum, so use it instead (they are back-to-back)
				//z = -1 * Proton
				//y = -1 * (Beam X Proton)
				//x = y cross z
			//Thus the production plane in the production frame is the XZ plane, and the normal vector is the Y-axis

		//Define production CM frame helicity axes
		TVector3 locHelicityZAxis_ProdCM = -1.0*locProtonP4_ProdCM.Vect().Unit();
		TVector3 locHelicityYAxis_ProdCM = -1.0*locBeamP4_ProdCM.Vect().Cross(locProtonP4_ProdCM.Vect()).Unit();
		TVector3 locHelicityXAxis_ProdCM = locHelicityYAxis_ProdCM.Cross(locHelicityZAxis_ProdCM).Unit();

		//Since the beam is in PARA configuration (Run 3185), the polarization vector is along the lab x-axis
			//Since the boost is in the z-direction, this vector is the same in the production CM frame
		TVector3 locPolUnit(1.0, 0.0, 0.0);

		//In the production CM frame, locPHI is the angle between the polarization vector and the production plane
		double locCosPHI = locBeamP4_ProdCM.Vect().Unit().Dot(locPolUnit.Cross(locHelicityYAxis_ProdCM));
		double locPHI = acos(locCosPHI); //reports phi between 0 and pi: sign ambiguity
		//Resolve the sign ambiguity
		double locSinPHI = locPolUnit.Dot(locHelicityYAxis_ProdCM);
		if(locSinPHI < 0.0)
			locPHI *= -1.0;

		//Now, we need the theta, phi angles between the omega decay plane and the production plane
		//The omega decay plane is defined by decay products in the omega CM frame
			//2 particles (vectors) define a plane.  
			//However, to conserve momentum, the third particle cannot be out of that plane (so must also be in it)
			//So, use the pi+ and the pi- to define the plane (pi0 measurement has less resolution)
		//By the way, for rho decays, the theta & phi angles are those of the pi+ in the rho CM frame, with respect to the helicity axes

		//boost pi+/- to omega CM frame
		TVector3 locBoostVector_OmegaCM = -1.0*(locOmegaP4_KinFit.BoostVector()); //negative due to coordinate system convention
		TLorentzVector locBeamP4_OmegaCM(locBeamP4_KinFit);
		locBeamP4_OmegaCM.Boost(locBoostVector_OmegaCM);
		TLorentzVector locProtonP4_OmegaCM(locProtonP4_KinFit);
		locProtonP4_OmegaCM.Boost(locBoostVector_OmegaCM);
		TLorentzVector locPiPlusP4_OmegaCM(locPiPlusP4_KinFit);
		locPiPlusP4_OmegaCM.Boost(locBoostVector_OmegaCM);
		TLorentzVector locPiMinusP4_OmegaCM(locPiMinusP4_KinFit);
		locPiMinusP4_OmegaCM.Boost(locBoostVector_OmegaCM);

		//Define omega CM frame helicity axes
			//These are defined the same way as before, but with the boost, the direction of the x & y axes has changed
		TVector3 locHelicityZAxis_OmegaCM = -1.0*locProtonP4_OmegaCM.Vect().Unit();
		TVector3 locHelicityYAxis_OmegaCM = -1.0*locBeamP4_OmegaCM.Vect().Cross(locProtonP4_OmegaCM.Vect()).Unit();
		TVector3 locHelicityXAxis_OmegaCM = locHelicityYAxis_OmegaCM.Cross(locHelicityZAxis_OmegaCM).Unit();

		//Compute the normal vector to the omega decay plane (pi+ x pi-)
		TVector3 locOmegaNormal = (locPiPlusP4_OmegaCM.Vect().Cross(locPiMinusP4_OmegaCM.Vect()));

		//Compute the theta angle to the omega decay plane
		double locCosTheta = locOmegaNormal.Dot(locHelicityZAxis_OmegaCM)/locOmegaNormal.Mag();
		double lcoTheta = acos(locCosTheta);

		//Compute the phi angle to the omega decay plane
		TVector3 locZCrossOmegaNormal = locHelicityZAxis_OmegaCM.Cross(locOmegaNormal);
		double locZCrossOmegaNormalMag = locZCrossOmegaNormal.Mag();
		double locCosPhi = locHelicityYAxis_OmegaCM.Dot(locZCrossOmegaNormal)/locZCrossOmegaNormalMag;
		double locPhi = acos(locCosPhi); //reports phi between 0 and pi: sign ambiguity
		//Resolve the sign ambiguity
		double locSinPhi = -1.0*locHelicityXAxis_OmegaCM.Dot(locZCrossOmegaNormal)/locZCrossOmegaNormalMag;
		if(locSinPhi < 0.0)
			locPhi *= -1.0;

		//Compute the "psi" angle: works at forward angles
		double locPsi = locPhi - locPHI;
		while(locPsi < -1.0*TMath::Pi())
			locPsi += 2.0*TMath::Pi();
		while(locPsi > TMath::Pi())
			locPsi -= 2.0*TMath::Pi();

		//result is defined by omega, only histogram if omega is unique
		if(locOmegaUniqueFlag)
		{
			dHist_OmegaPsi->Fill(180.0*locPsi/TMath::Pi());
			dHist_OmegaCosTheta->Fill(locCosTheta);
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

