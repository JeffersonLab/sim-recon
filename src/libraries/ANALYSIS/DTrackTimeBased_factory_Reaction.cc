// $Id$
//
//    File: DTrackTimeBased_factory_Reaction.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#include "DTrackTimeBased_factory_Reaction.h"

//------------------
// init
//------------------
jerror_t DTrackTimeBased_factory_Reaction::init(void)
{
	MAX_dReferenceTrajectoryPoolSize = 5;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrackTimeBased_factory_Reaction::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	dGeometry = locApplication ? locApplication->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber()) : NULL;
	dMagneticFieldMap = locApplication ? locApplication->GetBfield() : NULL;
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrackTimeBased_factory_Reaction::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
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
	for(int loc_i = 0; loc_i < int(locParticleComboBlueprints.size()); ++loc_i)
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
				locTrackTimeBased->AddAssociatedObject(locChargedTrack);
				_data.push_back(locTrackTimeBased);
				locCreatedTrackMap[locTrackPIDPair] = locTrackTimeBased;
			}
		}
	}

	return NOERROR;
}

DTrackTimeBased* DTrackTimeBased_factory_Reaction::Create_TrackTimeBased(const DChargedTrack* locChargedTrack, Particle_t locDesiredPID)
{
	//this assumes that the track reconstruction cannot get the incorrect charge
	const DChargedTrackHypothesis* locChargedTrackHypothesis = NULL;
	switch (locDesiredPID)
	{
		case PiPlus:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, PiPlus);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, PiPlus);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), PiPlus);
		case PiMinus:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, PiMinus);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), PiMinus);
		case KPlus:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, KPlus);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, KPlus);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), KPlus);
		case KMinus:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, KMinus);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), KMinus);
		case Proton:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, Proton);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, Proton);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), Proton);
		case AntiProton:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, AntiProton);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, AntiProton);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), AntiProton);
		case Positron:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, Positron);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, Positron);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, Positron);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), Positron);
		case Electron:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, Electron);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, Electron);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), Electron);
		case MuonPlus:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, MuonPlus);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KPlus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, MuonPlus);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(Proton);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, MuonPlus);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), MuonPlus);
		case MuonMinus:
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(PiMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, MuonMinus);
			locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(KMinus);
			if(locChargedTrackHypothesis != NULL)
				return Convert_ChargedTrack(locChargedTrackHypothesis, MuonMinus);
			return Convert_ChargedTrack(locChargedTrack->Get_BestFOM(), MuonMinus);
		default:
			return NULL;
	}
}

DReferenceTrajectory* DTrackTimeBased_factory_Reaction::Get_ReferenceTrajectoryResource(void)
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

DTrackTimeBased* DTrackTimeBased_factory_Reaction::Convert_ChargedTrack(const DChargedTrackHypothesis* locChargedTrackHypothesis, Particle_t locNewPID)
{
	DTrackTimeBased* locTrackTimeBased = new DTrackTimeBased();

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locChargedTrackHypothesis->GetT(locTrackTimeBasedVector);
	*locTrackTimeBased = *(locTrackTimeBasedVector[0]);
	locTrackTimeBased->AddAssociatedObject(locTrackTimeBasedVector[0]);

	locTrackTimeBased->setMass(ParticleMass(locNewPID));
	locTrackTimeBased->setPID(locNewPID);

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
jerror_t DTrackTimeBased_factory_Reaction::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTrackTimeBased_factory_Reaction::fini(void)
{
	for(size_t loc_i = 0; loc_i < dReferenceTrajectoryPool_All.size(); ++loc_i)
		delete dReferenceTrajectoryPool_All[loc_i];
	return NOERROR;
}


