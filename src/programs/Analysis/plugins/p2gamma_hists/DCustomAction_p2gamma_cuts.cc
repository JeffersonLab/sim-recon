// $Id$
//
//    File: DCustomAction_p2gamma_cuts.cc
// Created: Wed Jan 21 16:53:41 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p2gamma_cuts.h"

void DCustomAction_p2gamma_cuts::Initialize(JEventLoop* locEventLoop)
{

}

bool DCustomAction_p2gamma_cuts::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(0);

        // get beam photon energy and final state particles
	const DKinematicData* locBeamPhoton = NULL;
        deque<const DKinematicData*> locParticles;
        if(!Get_UseKinFitResultsFlag()) { //measured
		locBeamPhoton = locParticleComboStep->Get_InitialParticle_Measured();
                locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	}
	else {
		locBeamPhoton = locParticleComboStep->Get_InitialParticle();
		locParticleComboStep->Get_FinalParticles(locParticles);
	}
       
	// calculate missing mass
	DLorentzVector locMissingP4; 
	DLorentzVector locProtonP4Init(0,0,0,0.938);
	locMissingP4 += locProtonP4Init;
	locMissingP4 += locBeamPhoton->lorentzMomentum();
	DLorentzVector locSumInitP4 = locMissingP4;

	DLorentzVector locProtonP4;

	// calculate missing mass
	DLorentzVector loc2g_P4;
	for(size_t loc_i = 0; loc_i < 3; ++loc_i) {
		locMissingP4 -= locParticles[loc_i]->lorentzMomentum();
		
		if(locParticles[loc_i]->PID() == Gamma) 
			loc2g_P4 += locParticles[loc_i]->lorentzMomentum();
		else
			locProtonP4 += locParticles[loc_i]->lorentzMomentum();
	}

	double locDeltaPhi = (locProtonP4.Phi() - loc2g_P4.Phi())*180./TMath::Pi();
	if(locDeltaPhi > 360.) locDeltaPhi -= 360.;
	if(locDeltaPhi < 0.) locDeltaPhi += 360.;

	// require proton and pi0 are back-to-back
	if(locDeltaPhi < 175. || locDeltaPhi > 185.) 
		return false;
	
	// for pi0 candidates require recoil proton
	if(loc2g_P4.M() < 0.10 || loc2g_P4.M() > 0.16 || fabs(locMissingP4.M2()) < 0.05 || fabs(locMissingP4.E()) < 0.5) 
		return false;

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
