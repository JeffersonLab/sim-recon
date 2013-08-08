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
	MAX_dReferenceTrajectoryPoolSize = 5;

	//remember, charge sign could have flipped during track reconstruction
	deque<Particle_t> locPIDDeque;

	locPIDDeque.resize(4);
	locPIDDeque[0] = KPlus;  locPIDDeque[1] = Proton;  locPIDDeque[2] = PiMinus;  locPIDDeque[3] = KMinus;
	dParticleIDsToTry[PiPlus] = locPIDDeque;

	locPIDDeque.resize(3);
	locPIDDeque[0] = KMinus;  locPIDDeque[1] = PiPlus;  locPIDDeque[2] = KPlus;
	dParticleIDsToTry[PiMinus] = locPIDDeque;

	locPIDDeque.resize(4);
	locPIDDeque[0] = PiPlus;  locPIDDeque[1] = Proton;  locPIDDeque[2] = KMinus;  locPIDDeque[3] = PiMinus;
	dParticleIDsToTry[KPlus] = locPIDDeque;

	locPIDDeque.resize(4);
	locPIDDeque[0] = PiMinus;  locPIDDeque[1] = KPlus;  locPIDDeque[2] = PiPlus;  locPIDDeque[3] = Proton;
	dParticleIDsToTry[KMinus] = locPIDDeque;

	locPIDDeque.resize(3);
	locPIDDeque[0] = KPlus;  locPIDDeque[1] = PiPlus;  locPIDDeque[2] = KMinus;
	dParticleIDsToTry[Proton] = locPIDDeque;

	locPIDDeque.resize(4);
	locPIDDeque[0] = KMinus;  locPIDDeque[1] = PiMinus;  locPIDDeque[2] = Proton;  locPIDDeque[3] = KPlus;
	dParticleIDsToTry[AntiProton] = locPIDDeque;

	locPIDDeque.resize(5);
	locPIDDeque[0] = PiPlus;  locPIDDeque[1] = KPlus;  locPIDDeque[2] = Proton;  locPIDDeque[3] = PiMinus;  locPIDDeque[4] = KMinus;
	dParticleIDsToTry[Positron] = locPIDDeque;

	locPIDDeque.resize(4);
	locPIDDeque[0] = PiMinus;  locPIDDeque[1] = KMinus;  locPIDDeque[2] = PiPlus;  locPIDDeque[3] = KPlus;
	dParticleIDsToTry[Electron] = locPIDDeque;

	dParticleIDsToTry[MuonPlus] = dParticleIDsToTry[Positron];
	dParticleIDsToTry[MuonMinus] = dParticleIDsToTry[Electron];

	locPIDDeque.resize(3);
	locPIDDeque[0] = Proton;  locPIDDeque[1] = KPlus;  locPIDDeque[2] = PiPlus;
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
	//for each charged track, get the dchargedtrackhypothesis corresponding to the Analysis pid
	//if it doesn't exist, get it for the closest pid, and recalculate the PID FOM
	map<pair<const DChargedTrack*, Particle_t>, const DTrackTimeBased*> locCreatedTrackMap; //if a DChargedTrack was used a given PID, just copy the result
	pair<const DChargedTrack*, Particle_t> locTrackPIDPair;
	Particle_t locPID;
	DTrackTimeBased* locTrackTimeBased;
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
				locPID = locParticleComboBlueprintStep->Get_FinalParticleID(loc_k);
				if(locChargedTrack->Get_Hypothesis(locPID) != NULL)
					continue; //already exists, don't need to create a new one

				//check to see if new pid already created for this track
				locTrackPIDPair.first = locChargedTrack;
				locTrackPIDPair.second = locPID;
				if(locCreatedTrackMap.find(locTrackPIDPair) != locCreatedTrackMap.end())
					continue; //new one already created for this pair

				//create the DTrackTimeBased for the given PID
				locTrackTimeBased = Create_TrackTimeBased(locChargedTrack, locPID);
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
		return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), locDesiredPID);

	for(size_t loc_i = 0; loc_i < dParticleIDsToTry[locDesiredPID].size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(dParticleIDsToTry[locDesiredPID][loc_i]);
		if(locChargedTrackHypothesis != NULL)
			return Convert_ChargedTrack(locChargedTrackHypothesis, locDesiredPID);
	}
	return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), locDesiredPID);
}

DReferenceTrajectory* DTrackTimeBased_factory_Combo::Get_ReferenceTrajectoryResource(void)
{
	DReferenceTrajectory* locReferenceTrajectory;
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

DTrackTimeBased* DTrackTimeBased_factory_Combo::Convert_ChargedTrack(const DChargedTrackHypothesis* locChargedTrackHypothesis, Particle_t locNewPID)
{
	DTrackTimeBased* locTrackTimeBased = new DTrackTimeBased();

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locChargedTrackHypothesis->GetT(locTrackTimeBasedVector);
	*locTrackTimeBased = *(locTrackTimeBasedVector[0]);
	locTrackTimeBased->AddAssociatedObject(locTrackTimeBasedVector[0]);

	locTrackTimeBased->setMass(ParticleMass(locNewPID));
	locTrackTimeBased->setPID(locNewPID);
	locTrackTimeBased->setCharge(ParticleCharge(locNewPID));

	// configure the DReferenceTrajectory object
	DReferenceTrajectory *locReferenceTrajectory = Get_ReferenceTrajectoryResource();
	locReferenceTrajectory->SetMass(locTrackTimeBased->mass());
	locReferenceTrajectory->SetDGeometry(dGeometry);
	locReferenceTrajectory->Swim(locTrackTimeBased->position(), locTrackTimeBased->momentum(), locTrackTimeBased->charge());
	locTrackTimeBased->rt = locReferenceTrajectory;
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


