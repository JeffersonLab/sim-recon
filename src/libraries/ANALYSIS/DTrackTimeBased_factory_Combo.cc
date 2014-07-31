// $Id$
//
//    File: DTrackTimeBased_factory_Combo.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#include "DTrackTimeBased_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DTrackTimeBased_factory_Combo::init(void)
{
	MAX_dReferenceTrajectoryPoolSize = 10;

	//remember, charge sign could have flipped during track reconstruction
	deque<pair<Particle_t, bool> > locPIDDeque; //bool is true/false if should/shouldn't reswim

	//order of preference:
		//very similar mass & charge (e.g. pi -> mu, e) //no reswim (flight time virtually identical)
		//slightly similar mass & charge (e.g. pi -> K) //reswim: recalculate flight time
		//similar mass, different charge (e.g. pi+ -> pi-) //reswim! (different charge)
		//if to/from proton/deuteron: always reswim (significant mass difference)
		//don't need to list every PID: if none of the below PIDs are available, will choose the one with the best PID FOM and reswim

	//pi+
	locPIDDeque.resize(4);
	locPIDDeque[0] = pair<Particle_t, bool>(KPlus, true);  locPIDDeque[1] = pair<Particle_t, bool>(PiMinus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(KMinus, true);  locPIDDeque[3] = pair<Particle_t, bool>(Proton, true);
	dParticleIDsToTry[PiPlus] = locPIDDeque;

	//pi-
	locPIDDeque.resize(3);
	locPIDDeque[0] = pair<Particle_t, bool>(KMinus, true);  locPIDDeque[1] = pair<Particle_t, bool>(PiPlus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(KPlus, true);
	dParticleIDsToTry[PiMinus] = locPIDDeque;

	//K+
	locPIDDeque.resize(4);
	locPIDDeque[0] = pair<Particle_t, bool>(PiPlus, true);  locPIDDeque[1] = pair<Particle_t, bool>(KMinus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(PiMinus, true);  locPIDDeque[3] = pair<Particle_t, bool>(Proton, true);
	dParticleIDsToTry[KPlus] = locPIDDeque;

	//K-
	locPIDDeque.resize(4);
	locPIDDeque[0] = pair<Particle_t, bool>(PiMinus, true);  locPIDDeque[1] = pair<Particle_t, bool>(KPlus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(PiPlus, true);  locPIDDeque[3] = pair<Particle_t, bool>(Proton, true);
	dParticleIDsToTry[KMinus] = locPIDDeque;

	//p
	locPIDDeque.resize(3);
	locPIDDeque[0] = pair<Particle_t, bool>(KPlus, true);  locPIDDeque[1] = pair<Particle_t, bool>(PiPlus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(KMinus, true);
	dParticleIDsToTry[Proton] = locPIDDeque;

	//pbar
	locPIDDeque.resize(4);
	locPIDDeque[0] = pair<Particle_t, bool>(KMinus, true);  locPIDDeque[1] = pair<Particle_t, bool>(Proton, true);
	locPIDDeque[2] = pair<Particle_t, bool>(PiMinus, true);  locPIDDeque[3] = pair<Particle_t, bool>(KPlus, true);
	dParticleIDsToTry[AntiProton] = locPIDDeque;

	//e+
	locPIDDeque.resize(4);
	locPIDDeque[0] = pair<Particle_t, bool>(PiPlus, false);  locPIDDeque[1] = pair<Particle_t, bool>(KPlus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(PiMinus, true);  locPIDDeque[3] = pair<Particle_t, bool>(KMinus, true);
	dParticleIDsToTry[Positron] = locPIDDeque;

	//e-
	locPIDDeque.resize(4);
	locPIDDeque[0] = pair<Particle_t, bool>(PiMinus, false);  locPIDDeque[1] = pair<Particle_t, bool>(KMinus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(PiPlus, true);  locPIDDeque[3] = pair<Particle_t, bool>(KPlus, true);
	dParticleIDsToTry[Electron] = locPIDDeque;

	//mu+
	dParticleIDsToTry[MuonPlus] = dParticleIDsToTry[Positron];

	//mu-
	dParticleIDsToTry[MuonMinus] = dParticleIDsToTry[Electron];

	//D
	locPIDDeque.resize(3);
	locPIDDeque[0] = pair<Particle_t, bool>(Proton, true);  locPIDDeque[1] = pair<Particle_t, bool>(KPlus, true);
	locPIDDeque[2] = pair<Particle_t, bool>(PiPlus, true);
	dParticleIDsToTry[Deuteron] = locPIDDeque;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory_Combo::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	dGeometry = locApplication ? locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber()) : NULL;
	dMagneticFieldMap = locApplication ? locApplication->GetBfield() : NULL;
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackTimeBased_factory_Combo::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dReferenceTrajectoryPool_All.size() > MAX_dReferenceTrajectoryPoolSize){
		for(size_t loc_i = MAX_dReferenceTrajectoryPoolSize; loc_i < dReferenceTrajectoryPool_All.size(); ++loc_i)
			delete dReferenceTrajectoryPool_All[loc_i];
		dReferenceTrajectoryPool_All.resize(MAX_dReferenceTrajectoryPoolSize);
	}
	dReferenceTrajectoryPool_Available = dReferenceTrajectoryPool_All;

 	vector<const DParticleComboBlueprint*> locParticleComboBlueprints;
	locEventLoop->Get(locParticleComboBlueprints);

	//Create New DTrackTimeBased Objects (as needed) //for PIDs not available in the original objects
		//if DTrackTimeBased doesn't exist, get it for the closest pid
	map<pair<const DChargedTrack*, Particle_t>, const DTrackTimeBased*> locCreatedTrackMap; //if a DChargedTrack was used a given PID, just copy the result
	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	{
		const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_k)
			{
				if((!locParticleComboBlueprintStep->Is_FinalParticleDetected(loc_k)) || (!locParticleComboBlueprintStep->Is_FinalParticleCharged(loc_k)))
					continue;

				const DChargedTrack* locChargedTrack = static_cast<const DChargedTrack*>(locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_k));

				//check if pid already exists for this track
				Particle_t locPID = locParticleComboBlueprintStep->Get_FinalParticleID(loc_k);
				if(locChargedTrack->Get_Hypothesis(locPID) != NULL)
					continue; //already exists, don't need to create a new one

				//check to see if new pid already created for this track
				pair<const DChargedTrack*, Particle_t> locTrackPIDPair(locChargedTrack, locPID);
				if(locCreatedTrackMap.find(locTrackPIDPair) != locCreatedTrackMap.end())
					continue; //new one already created for this pair

				//create the DTrackTimeBased for the given PID
				DTrackTimeBased* locTrackTimeBased = Create_TrackTimeBased(locChargedTrack, locPID);
				if(locTrackTimeBased == NULL)
					continue;
				locTrackTimeBased->AddAssociatedObject(locChargedTrack);
				_data.push_back(locTrackTimeBased);
				locCreatedTrackMap[locTrackPIDPair] = locTrackTimeBased;
			}
		}
	}

	return NOERROR;
}

DTrackTimeBased* DTrackTimeBased_factory_Combo::Create_TrackTimeBased(const DChargedTrack* locChargedTrack, Particle_t locDesiredPID)
{
	if(dParticleIDsToTry.find(locDesiredPID) == dParticleIDsToTry.end())
		return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), locDesiredPID, true);

	for(size_t loc_i = 0; loc_i < dParticleIDsToTry[locDesiredPID].size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(dParticleIDsToTry[locDesiredPID][loc_i].first);
		if(locChargedTrackHypothesis != NULL)
			return Convert_ChargedTrack(locChargedTrackHypothesis, locDesiredPID, dParticleIDsToTry[locDesiredPID][loc_i].second);
	}
	return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), locDesiredPID, true);
}

DReferenceTrajectory* DTrackTimeBased_factory_Combo::Get_ReferenceTrajectoryResource(void)
{
	DReferenceTrajectory* locReferenceTrajectory = NULL;
	if(dReferenceTrajectoryPool_Available.empty())
	{
		locReferenceTrajectory = new DReferenceTrajectory(dMagneticFieldMap);
		dReferenceTrajectoryPool_All.push_back(locReferenceTrajectory);
	}
	else
	{
		locReferenceTrajectory = dReferenceTrajectoryPool_Available.back();
		locReferenceTrajectory->Reset();
		dReferenceTrajectoryPool_Available.pop_back();
	}
	return locReferenceTrajectory;
}

DTrackTimeBased* DTrackTimeBased_factory_Combo::Convert_ChargedTrack(const DChargedTrackHypothesis* locChargedTrackHypothesis, Particle_t locNewPID, bool locSwimFlag)
{
	DTrackTimeBased* locTrackTimeBased = new DTrackTimeBased();

 	const DTrackTimeBased* locOriginalTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingleT(locOriginalTrackTimeBased);
	*locTrackTimeBased = *locOriginalTrackTimeBased;

	locTrackTimeBased->setMass(ParticleMass(locNewPID));
	locTrackTimeBased->setPID(locNewPID);
	locTrackTimeBased->setCharge(ParticleCharge(locNewPID));

	// swim the DReferenceTrajectory object if necessary
	if(locSwimFlag)
	{
		DReferenceTrajectory *locReferenceTrajectory = Get_ReferenceTrajectoryResource();
		locReferenceTrajectory->SetMass(locTrackTimeBased->mass());
		locReferenceTrajectory->SetDGeometry(dGeometry);
		locReferenceTrajectory->Swim(locTrackTimeBased->position(), locTrackTimeBased->momentum(), locTrackTimeBased->charge());
		locTrackTimeBased->rt = locReferenceTrajectory;
	}
	else
		locTrackTimeBased->rt = NULL;

	locTrackTimeBased->AddAssociatedObject(locOriginalTrackTimeBased);
	locTrackTimeBased->AddAssociatedObject(locChargedTrackHypothesis);

	return locTrackTimeBased;
}

//------------------
// erun
//------------------
jerror_t DTrackTimeBased_factory_Combo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackTimeBased_factory_Combo::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReferenceTrajectoryPool_All.size(); ++loc_i)
		delete dReferenceTrajectoryPool_All[loc_i];
	return NOERROR;
}


