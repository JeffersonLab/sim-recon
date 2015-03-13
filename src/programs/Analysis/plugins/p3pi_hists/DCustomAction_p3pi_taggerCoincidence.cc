// $Id$
//
//    File: DCustomAction_p3pi_taggerCoincidence.cc
// Created: Thu Jan 22 08:06:18 EST 2015
// Creator: jrsteven (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "DCustomAction_p3pi_taggerCoincidence.h"

void DCustomAction_p3pi_taggerCoincidence::Initialize(JEventLoop* locEventLoop)
{

	// get PID algos
	const DParticleID* locParticleID = NULL;
        locEventLoop->GetSingle(locParticleID);
	dParticleID = locParticleID;

	japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
	{
		//Required: Create a folder in the ROOT output file that will contain all of the output ROOT objects (if any) for this action.
			//If another thread has already created the folder, it just changes to it. 
		CreateAndChangeTo_ActionDirectory();

		dMatch_E_DeltaT_All = GetOrCreate_Histogram<TH2I>("Match_E_DeltaT_All", "Match Charged Track - TAGGER: Energy vs #Delta t; #Delta t; Energy", 100, -50, 50, 200, 2., 12.);
		dMatch_E_DeltaT_SC = GetOrCreate_Histogram<TH2I>("Match_E_DeltaT_SC", "Match SC - TAGGER: Energy vs #Delta t; #Delta t; Energy", 100, -50, 50, 200, 2., 12.);

	}
	japp->RootUnLock(); //RELEASE ROOT LOCK!!
}

bool DCustomAction_p3pi_taggerCoincidence::Perform_Action(JEventLoop* locEventLoop, const DParticleCombo* locParticleCombo)
{
	const DDetectorMatches* locDetectorMatches = NULL;
        locEventLoop->GetSingle(locDetectorMatches);
		
	const DParticleComboStep* locParticleComboStep = locParticleCombo->Get_ParticleComboStep(1);
	if(locParticleComboStep->Get_InitialParticleID() != omega)
		return false;

	// get beam photon time
	const DKinematicData* locBeamPhoton =  locParticleCombo->Get_ParticleComboStep(0)->Get_InitialParticle();
	double locBeamPhotonTime = locBeamPhoton->time();
	double locBeamPhotonEnergy = locBeamPhoton->energy();

	// get final state particles
	deque<const DKinematicData*> locParticles;
	if(!Get_UseKinFitResultsFlag()) //measured
		locParticleComboStep->Get_FinalParticles_Measured(locParticles);
	else
		locParticleComboStep->Get_FinalParticles(locParticles);

	bool locGoodTime = false;

	// loop final state particles
	for(size_t loc_i = 0; loc_i < locParticles.size(); ++loc_i) {
		Particle_t locPID = locParticleComboStep->Get_FinalParticleID(loc_i);
		if(locPID != PiPlus && locPID != PiMinus) continue;

		// get time based track
		const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboStep->Get_FinalParticle_SourceObject(loc_i));
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_BestFOM();
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);
		
		double locTime = locChargedTrackHypothesis->time();
		double locDeltaT = locTime - locBeamPhotonTime;

		japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
		{
			dMatch_E_DeltaT_All->Fill(locDeltaT, locBeamPhotonEnergy);
		}
		japp->RootUnLock(); //RELEASE ROOT LOCK!!

		// get match to SC
		DSCHitMatchParams locSCHitMatchParams;
		bool foundSC = dParticleID->Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams);
                if(foundSC){
                        double locSCTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime;
			double locDeltaT = locSCTime - locBeamPhotonTime;
			
			japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
			{
				dMatch_E_DeltaT_SC->Fill(locDeltaT, locBeamPhotonEnergy);
			}
			japp->RootUnLock(); //RELEASE ROOT LOCK!!

			// reject out of time beam photons
			if(fabs(locDeltaT) < dDeltaTMax) {
				locGoodTime = true;
			}
		}
	}

	return locGoodTime; //return false if you want to use this action to apply a cut (and it fails the cut!)
}
