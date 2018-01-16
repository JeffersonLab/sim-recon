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
        auto locParticles = Get_UseKinFitResultsFlag() ? locParticleCombo->Get_FinalParticles(Get_Reaction(), false, false) : locParticleCombo->Get_FinalParticles_Measured(Get_Reaction());
        if(!Get_UseKinFitResultsFlag()) { //measured
		locBeamPhoton = locParticleComboStep->Get_InitialParticle_Measured();
	}
	else {
		locBeamPhoton = locParticleComboStep->Get_InitialParticle();
	}
       
	// calculate missing mass
	DLorentzVector locMissingP4; 
	DLorentzVector locProtonP4Init(0,0,0,0.938);
	locMissingP4 += locProtonP4Init;
	locMissingP4 += locBeamPhoton->lorentzMomentum();
	DLorentzVector locSumInitP4 = locMissingP4;

	DLorentzVector locProtonP4;
	DLorentzVector locPiPlusP4;
	DLorentzVector locPiMinusP4;
	DLorentzVector locCandidateP4;

	// calculate missing mass
	DLorentzVector loc2g_P4;
	for(size_t loc_i = 0; loc_i < 5; ++loc_i) {
		locMissingP4 -= locParticles[loc_i]->lorentzMomentum();
		
		if(locParticles[loc_i]->PID() == Gamma) 
			loc2g_P4 += locParticles[loc_i]->lorentzMomentum();
		else if(locParticles[loc_i]->PID() == PiPlus)
			locPiPlusP4 += locParticles[loc_i]->lorentzMomentum();
		else if(locParticles[loc_i]->PID() == PiMinus)
			locPiMinusP4 += locParticles[loc_i]->lorentzMomentum();
		else
			locProtonP4 += locParticles[loc_i]->lorentzMomentum();
	}

	locCandidateP4 = loc2g_P4 + locPiPlusP4 + locPiMinusP4;

	double locDeltaPhi = (locProtonP4.Phi() - locCandidateP4.Phi())*180./TMath::Pi();
	if(locDeltaPhi > 360.) locDeltaPhi -= 360.;
	if(locDeltaPhi < 0.) locDeltaPhi += 360.;

	// require proton and pi0 are back-to-back
	if(locDeltaPhi < 160. || locDeltaPhi > 200.) 
		return false;
	// exclusive kinematics cuts
	if(fabs(locMissingP4.M2()) > 0.1 || fabs(locMissingP4.E()) > 1.0) 
		return false;

	return true; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
