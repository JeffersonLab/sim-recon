// $Id$
//
//    File: DCustomAction_HistMass_b1_1235.cc
// Created: Mon Dec  2 12:18:47 EST 2013
// Creator: pmatt (on Darwin pmattLaptop 10.8.0 i386)
//

#include "DCustomAction_HistMass_b1_1235.h"

void DCustomAction_HistMass_b1_1235::Initialize(JEventLoop* locEventLoop)
{
	//Optional: Create histograms and/or modify member variables.
	//Create any histograms/trees/etc. within a ROOT lock. 
		//This is so that when running multithreaded, only one thread is writing to the ROOT file at a time. 

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
		string locHistTitle = string(";") + ParticleName_ROOT(PiPlus) + ParticleName_ROOT(omega) + string(" Invariant Mass (GeV/c^{2});# Combos / 2 MeV/c^{2}");
		if(gDirectory->Get("InvariantMass") == NULL) //check to see if already created by another thread
			dMassHist = new TH1D("InvariantMass", locHistTitle.c_str(), 600, 0.6, 1.8);
		else //already created by another thread
			dMassHist = static_cast<TH1D*>(gDirectory->Get("InvariantMass"));
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_HistMass_b1_1235::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	//Optional: check whether the user wanted to use the kinematic fit results when performing this action
	bool locUseKinFitResultsFlag = Get_UseKinFitResultsFlag();

	//Get the DParticleCombo objects for which this action has been previously executed.
		//This is useful for determining whether filling a histogram will result in double-counting. 
	deque<pair<const DParticleCombo*, bool> > locPreviousParticleCombos;
	Get_PreviousParticleCombos(locPreviousParticleCombos);

	const DParticleComboStep* locParticleComboStep0 = locParticleCombo->Get_ParticleComboStep(0);
	const DParticleComboStep* locParticleComboStep1 = locParticleCombo->Get_ParticleComboStep(1);
	const DParticleComboStep* locParticleComboStep2 = locParticleCombo->Get_ParticleComboStep(2);

	const DKinematicData* locPiPlus1 = locUseKinFitResultsFlag ? locParticleComboStep0->Get_FinalParticle(2) : locParticleComboStep0->Get_FinalParticle_Measured(2);

	const DKinematicData* locPiPlus2 = locUseKinFitResultsFlag ? locParticleComboStep1->Get_FinalParticle(0) : locParticleComboStep1->Get_FinalParticle_Measured(0);
	const DKinematicData* locPiMinus2 = locUseKinFitResultsFlag ? locParticleComboStep1->Get_FinalParticle(1) : locParticleComboStep1->Get_FinalParticle_Measured(1);

	const DKinematicData* locPhoton1 = locUseKinFitResultsFlag ? locParticleComboStep2->Get_FinalParticle(0) : locParticleComboStep2->Get_FinalParticle_Measured(0);
	const DKinematicData* locPhoton2 = locUseKinFitResultsFlag ? locParticleComboStep2->Get_FinalParticle(1) : locParticleComboStep2->Get_FinalParticle_Measured(1);

	set<const DKinematicData*> locCurrentParticles;
	locCurrentParticles.insert(locPiMinus2);
	locCurrentParticles.insert(locPiPlus1);
	locCurrentParticles.insert(locPiPlus2);
	locCurrentParticles.insert(locPhoton1);
	locCurrentParticles.insert(locPhoton2);

	//if new event: clear past particles, else check if duplicate
	if(locPreviousParticleCombos.empty())
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
	if(!locUseKinFitResultsFlag) //measured
	{
		locP4 += locPiMinus2->lorentzMomentum();
		locP4 += locPiPlus1->lorentzMomentum() + locPiPlus2->lorentzMomentum();
		locP4 += locPhoton1->lorentzMomentum() + locPhoton2->lorentzMomentum();
	}
	else //kinfit
	{
		locP4 += locPiMinus2->lorentzMomentum();
		locP4 += locPiPlus1->lorentzMomentum() + locPiPlus2->lorentzMomentum();
		const DKinematicData* locPiZero = locParticleComboStep1->Get_FinalParticle(2);
		locP4 += locPiZero->lorentzMomentum();
	}
	double locInvariantMass = locP4.M();

	//Optional: Fill histograms
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Fill any histograms here
		dMassHist->Fill(locInvariantMass);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
