// $Id$
//
//    File: DCustomAction_CutExtraPi0.cc
// Created: Sun Jun 28 15:10:49 EDT 2015
// Creator: pmatt (on Darwin Pauls-MacBook-Pro-2.local 13.4.0 i386)
//

#include "DCustomAction_CutExtraPi0.h"

void DCustomAction_CutExtraPi0::Initialize(JEventLoop* locEventLoop)
{
	//CREATE THE HISTOGRAMS
	//Since we are creating histograms, the contents of gDirectory will be modified: must use JANA-wide ROOT lock
	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		// Optional: Useful utility functions.
		locEventLoop->GetSingle(dAnalysisUtilities);

		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dHist_Pi0InvariantMass = GetOrCreate_Histogram<TH1I>("InvariantMass_Pi0", ";#gamma#gamma Invariant Mass (GeV/c^{2})", 600, 0.0, 0.3);
	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_CutExtraPi0::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	if(Get_NumPreviousParticleCombos() == 0)
		dPreviousSourceObjects.clear();

	vector<const DNeutralParticle*> locUnusedNeutralParticles;
	dAnalysisUtilities->Get_UnusedNeutralParticles(locEventLoop, locParticleCombo, locUnusedNeutralParticles);

	bool locPassesCutFlag = true;
	vector<double> locInvariantMasses; //collect them all, then fill at once (MUCH fewer locks)
	for(size_t loc_i = 0; loc_i < locUnusedNeutralParticles.size(); ++loc_i)
	{
		const DNeutralParticleHypothesis* locPhoton1 = locUnusedNeutralParticles[loc_i]->Get_Hypothesis(Gamma);
		for(size_t loc_j = loc_i + 1; loc_j < locUnusedNeutralParticles.size(); ++loc_j)
		{
			const DNeutralParticleHypothesis* locPhoton2 = locUnusedNeutralParticles[loc_j]->Get_Hypothesis(Gamma);
			DLorentzVector locP4 = locPhoton1->lorentzMomentum() + locPhoton2->lorentzMomentum();

			double locInvariantMass = locP4.M();
			if((locInvariantMass >= dLowMassCut) && (locInvariantMass <= dHighMassCut))
				locPassesCutFlag = false;

			//check to see if this combination has already been histogrammed (e.g. on previous combo)
			set<const DNeutralParticleHypothesis*> locPhotonSet;
			locPhotonSet.insert(locPhoton1);
			locPhotonSet.insert(locPhoton2);
			if(dPreviousSourceObjects.find(locPhotonSet) != dPreviousSourceObjects.end())
				continue; //dupe: already histed!
			dPreviousSourceObjects.insert(locPhotonSet);

			//save to vector to fill hist later
			locInvariantMasses.push_back(locInvariantMass);
		}
	}

	//FILL HISTOGRAMS
	//Since we are filling histograms local to this action, it will not interfere with other ROOT operations: can use action-wide ROOT lock
	//Note, the mutex is unique to this DReaction + action_string combo: actions of same class with different hists will have a different mutex
	Lock_Action(); //ACQUIRE ROOT LOCK!!
	{
		for(size_t loc_i = 0; loc_i < locInvariantMasses.size(); ++loc_i)
			dHist_Pi0InvariantMass->Fill(locInvariantMasses[loc_i]);
	}
	Unlock_Action(); //RELEASE ROOT LOCK!!

	return locPassesCutFlag; //return false if you want to use this action to apply a cut (and it fails the cut!)
}

