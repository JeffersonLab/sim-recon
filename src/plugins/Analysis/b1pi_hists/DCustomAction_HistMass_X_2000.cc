// $Id$
//
//    File: DCustomAction_HistMass_X_2000.cc
// Created: Mon Dec  2 12:18:54 EST 2013
// Creator: pmatt (on Darwin pmattLaptop 10.8.0 i386)
//

#include "DCustomAction_HistMass_X_2000.h"

void DCustomAction_HistMass_X_2000::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 

	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Optional: Useful utility functions.
		locEventLoop->GetSingle(dAnalysisUtilities);

		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		// Optional: Create a ROOT subfolder.
			//If another thread has already created the folder, it just changes to it. 
		// CreateAndChangeTo_Directory("MyDirName", "MyDirTitle");
			//make sub-directory content here
		// gDirectory->cd(".."); //return to the action directory

		//	(Optional) Example: Create a histogram.
		string locHistTitle = string(";") + ParticleName_ROOT(PiMinus) + ParticleName_ROOT(PiPlus) + ParticleName_ROOT(omega) + string(" Invariant Mass (GeV/c^{2});# Combos / 2 MeV/c^{2}");
		if(gDirectory->Get("InvariantMass") == NULL) //check to see if already created by another thread
			dMassHist = new TH1I("InvariantMass", locHistTitle.c_str(), 500, 1.5, 2.5);
		else //already created by another thread
			dMassHist = static_cast<TH1I*>(gDirectory->Get("InvariantMass"));
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_HistMass_X_2000::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Optional: check whether the user wanted to use the kinematic fit results when performing this action
//	bool locUseKinFitResultsFlag = Get_UseKinFitResultsFlag();

	const DParticleComboStep* locParticleComboStep0 = locParticleCombo->Get_ParticleComboStep(0);
	const DParticleComboStep* locParticleComboStep1 = locParticleCombo->Get_ParticleComboStep(1);
	const DParticleComboStep* locParticleComboStep2 = locParticleCombo->Get_ParticleComboStep(2);

	//first get measured particle objects, and check to see if combination of particles is unique
	const DKinematicData* locPiMinus1 = locParticleComboStep0->Get_FinalParticle_Measured(1);
	const DKinematicData* locPiPlus1 = locParticleComboStep0->Get_FinalParticle_Measured(2);

	const DKinematicData* locPiPlus2 = locParticleComboStep1->Get_FinalParticle_Measured(0);
	const DKinematicData* locPiMinus2 = locParticleComboStep1->Get_FinalParticle_Measured(1);

	const DKinematicData* locPhoton1 = locParticleComboStep2->Get_FinalParticle_Measured(0);
	const DKinematicData* locPhoton2 = locParticleComboStep2->Get_FinalParticle_Measured(1);

	set<const DKinematicData*> locCurrentParticles;
	locCurrentParticles.insert(locPiMinus1);
	locCurrentParticles.insert(locPiMinus2);
	locCurrentParticles.insert(locPiPlus1);
	locCurrentParticles.insert(locPiPlus2);
	locCurrentParticles.insert(locPhoton1);
	locCurrentParticles.insert(locPhoton2);

	//if new event: clear past particles, else check if duplicate
	if(Get_NumPreviousParticleCombos() == 0)
		dPastParticles.clear();
	else //have had previous combos for this event, check to make sure particles used to compute this quantity aren't duplicate
	{
		for(size_t loc_i = 0; loc_i < dPastParticles.size(); ++loc_i)
		{
			if(locCurrentParticles != dPastParticles[loc_i])
				continue;
			return true; //duplicate combo of particles, don't fill histogram!
		}
	}
	dPastParticles.push_back(locCurrentParticles);

	DLorentzVector locP4;
//	if(!locUseKinFitResultsFlag || (locParticleCombo->Get_KinFitResults() == NULL)) //measured or kinfit failed to converge
	{
		locP4 += locPiMinus1->lorentzMomentum() + locPiMinus2->lorentzMomentum();
		locP4 += locPiPlus1->lorentzMomentum() + locPiPlus2->lorentzMomentum();
		locP4 += locPhoton1->lorentzMomentum() + locPhoton2->lorentzMomentum();
	}
/*	else //kinfit: get kinfit objects first
	{
		locPiMinus1 = locParticleComboStep0->Get_FinalParticle(1);
		locPiPlus1 = locParticleComboStep0->Get_FinalParticle(2);

		locPiPlus2 = locParticleComboStep1->Get_FinalParticle(0);
		locPiMinus2 = locParticleComboStep1->Get_FinalParticle(1);
		const DKinematicData* locPiZero = locParticleComboStep1->Get_FinalParticle(2);

		locP4 += locPiMinus1->lorentzMomentum() + locPiMinus2->lorentzMomentum();
		locP4 += locPiPlus1->lorentzMomentum() + locPiPlus2->lorentzMomentum();
		locP4 += locPiZero->lorentzMomentum();
	}*/
	double locInvariantMass = locP4.M();

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		// Fill any histograms here
		dMassHist->Fill(locInvariantMass);
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

