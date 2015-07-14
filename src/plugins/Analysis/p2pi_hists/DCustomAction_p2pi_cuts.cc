// $Id$
//
//    File: DCustomAction_p2pi_hists.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2pi_cuts.h"

void DCustomAction_p2pi_cuts::Initialize(JEventLoop* locEventLoop)
{
	// check if a particle is missing
	Get_Reaction()->Get_MissingPID(dMissingPID);
}

bool DCustomAction_p2pi_cuts::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	// should only have one reaction step
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);

	// get beam photon energy and final state particles
        deque<const DKinematicData*> locParticles;
        if(!Get_UseKinFitResultsFlag()) { //measured
                locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	}
	else {
		locParticleComboStep->Get_FinalParticles(locParticles);
	}

	// calculate 2pi P4
	DLorentzVector locP4_2pi;
	for(size_t loc_i = 0; loc_i < 3; ++loc_i) {
		if(locParticles[loc_i] == NULL) continue; // missing proton
                if(locParticles[loc_i]->PID() == PiPlus || locParticles[loc_i]->PID() == PiMinus)
                        locP4_2pi += locParticles[loc_i]->lorentzMomentum();
        }

	// calculate missing P4
	DLorentzVector locMissingP4 = dAnalysisUtilities->Calc_MissingP4(locParticleCombo,Get_UseKinFitResultsFlag());

	if(dMissingPID != Proton) {
		
		if(fabs(locMissingP4.E()) > 0.2)
			return false;
		
		if(locP4_2pi.M() < 0.6 || locP4_2pi.M() > 0.9)
			return false;				
	}

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
