// $Id$
//
//    File: DCustomAction_p3pi_Pi0Cuts.cc
// Created: Thu Jan 22 11:19:47 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p3pi_Pi0Cuts.h"

void DCustomAction_p3pi_Pi0Cuts::Initialize(JEventLoop* locEventLoop)
{
	
}

bool DCustomAction_p3pi_Pi0Cuts::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(2);
        if(locParticleComboStep->Get_InitialParticleID() != Pi0)
                return false;
	
	// get final state particles
        deque<const DKinematicData*> locParticles;
        if(!Get_UseKinFitResultsFlag()) //measured
                locParticleComboStep->Get_FinalParticles_Measured(locParticles);
        else
                locParticleComboStep->Get_FinalParticles(locParticles);

	int nFCAL = 0;
	
        // loop final state particles
        for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {
                if(locParticleComboStep->Is_FinalParticleMissing(loc_i)) continue;
		
		// get shower object
		const DNeutralShower* locNeutralShower = static_cast<const DNeutralShower*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
		if(locNeutralShower == NULL) 
			continue;
		
		// count # of FCAL photons and set separate thresholds on FCAL and BCAL energies
		if(locNeutralShower->dDetectorSystem == SYS_FCAL) {
			nFCAL++;
		}
	}

	// require 1 or 2 photons in FCAL (specify parameter in DReaction setup)
	if(nFCAL != dMinFCAL) 
		return false;

	// passed all cuts
	return true; 
}
