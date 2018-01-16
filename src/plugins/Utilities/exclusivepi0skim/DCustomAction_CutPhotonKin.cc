// $Id$
//
//    File: DCustomAction_CutPhotonKin.cc
// Created: Fri Jul 18 12:51:03 EDT 2014
// Creator: jrsteven (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#include "DCustomAction_CutPhotonKin.h"

void DCustomAction_CutPhotonKin::Initialize(JEventLoop* locEventLoop)
{
	
}

bool DCustomAction_CutPhotonKin::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
    auto locParticles = locParticleCombo->Get_FinalParticles_Measured(Get_Reaction());

	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i)
	{
		if(ParticleCharge(locParticles[loc_i]->PID()) == 0)
		{
			const DNeutralParticleHypothesis* locNeutralParticleHypothesis = static_cast<const DNeutralParticleHypothesis*>(locParticles[loc_i]);
			const DNeutralShower* locNeutralShower = locNeutralParticleHypothesis->Get_NeutralShower();

			// make BCAL cut on photon here
			if(locNeutralShower->dDetectorSystem != SYS_BCAL) return false;
            // Require each photon to have at least 500 MeV of energy
            if(locNeutralShower->dEnergy < 0.5)  return false;
		}
	}

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
