// $Id$
//
//    File: DChargedTrackHypothesis_factory_Combo.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DChargedTrackHypothesis_factory_Combo.h"
using namespace std;
using namespace jana;

//------------------
// init
//------------------
jerror_t DChargedTrackHypothesis_factory_Combo::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrackHypothesis_factory_Combo::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses); //make sure that brun() is called for the default factory!!!
	dChargedTrackHypothesisFactory = static_cast<DChargedTrackHypothesis_factory*>(locEventLoop->GetFactory("DChargedTrackHypothesis"));
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrackHypothesis_factory_Combo::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DChargedTrackHypothesis_factory_Combo::evnt()");
#endif

 	vector<const DParticleComboBlueprint*> locParticleComboBlueprints;
	locEventLoop->Get(locParticleComboBlueprints);

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector, "Combo");

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Combo");

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches, "Combo");

	map<const DEventRFBunch*, deque<pair<const DChargedTrack*, Particle_t> > > locCreatedParticleMap; //don't create if already done!

	DChargedTrackHypothesis* locNewChargedTrackHypothesis;
	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	{
		const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];

		//select the corresponding rf bunch
		const DEventRFBunch* locEventRFBunch = NULL;
		for(size_t loc_j = 0; loc_j < locEventRFBunches.size(); ++loc_j)
		{
			if(!locEventRFBunches[loc_j]->IsAssociated(locParticleComboBlueprint))
				continue;
			locEventRFBunch = locEventRFBunches[loc_j];
			break;
		}
		if(locEventRFBunch == NULL)
		{
			cout << "SOMETHING IS VERY WRONG IN DParticleCombo_factory_PreKinFit.cc" << endl;
			abort();
		}

		if(locCreatedParticleMap.find(locEventRFBunch) == locCreatedParticleMap.end())
			locCreatedParticleMap[locEventRFBunch] = deque<pair<const DChargedTrack*, Particle_t> >();

		for(size_t loc_j = 0; loc_j < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_k)
			{
				const JObject* locSourceObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_k);
				if(locSourceObject == NULL)
					continue; //missing or decaying
				const DChargedTrack* locChargedTrack = dynamic_cast<const DChargedTrack*>(locSourceObject);
				if(locChargedTrack == NULL)
					continue; //neutral
				Particle_t locPID = locParticleComboBlueprintStep->Get_FinalParticleID(loc_k);

				//see if already created for this pid/track/rf bunch pair
				deque<pair<const DChargedTrack*, Particle_t> > locParticlePairs = locCreatedParticleMap[locEventRFBunch];
				bool locAlreadyCreatedFlag = false;
				for(size_t loc_l = 0; loc_l < locParticlePairs.size(); ++loc_l)
				{
					if((locParticlePairs[loc_l].first != locChargedTrack) || (locParticlePairs[loc_l].second != locPID))
						continue;
					locAlreadyCreatedFlag = true; //e.g., several combos have the same rf bunch
					break;
				}
				if(locAlreadyCreatedFlag)
					continue;

				//see if DChargedTrackHypothesis with the desired PID was created by the default factory
				const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
				if(locChargedTrackHypothesis != NULL)
				{
					//yes it was: create new object with same PID (so that is registered with the combo factory, and because rf bunch could be different)
					const DTrackTimeBased* locTrackTimeBased = NULL;
					locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);
					locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBased, locDetectorMatches, locEventRFBunch);
					locNewChargedTrackHypothesis->AddAssociatedObject(locEventRFBunch);
					locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
					_data.push_back(locNewChargedTrackHypothesis);
					locCreatedParticleMap[locEventRFBunch].push_back(pair<const DChargedTrack*, Particle_t>(locChargedTrack, locPID));
					continue;
				}

				//no DChargedTrackHypothesis with this PID: get track info from DTrackTimeBased "Combo"-tag objects
				for(size_t loc_l = 0; loc_l < locTrackTimeBasedVector.size(); ++loc_l)
				{
					const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_l];
					if(locTrackTimeBased->PID() != locPID)
						continue;
					const DChargedTrack* locTimeBasedSourceChargedTrack = NULL;
					locTrackTimeBased->GetSingleT(locTimeBasedSourceChargedTrack);
					if(locTimeBasedSourceChargedTrack != locChargedTrack)
						continue;
					//correct DTrackTimeBased grabbed for this source object: create new DChargedTrackHypothesis object
					locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBased, locDetectorMatches, locEventRFBunch);
					locNewChargedTrackHypothesis->AddAssociatedObject(locEventRFBunch);
					locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
					_data.push_back(locNewChargedTrackHypothesis);
					locCreatedParticleMap[locEventRFBunch].push_back(pair<const DChargedTrack*, Particle_t>(locChargedTrack, locPID));
					break;
				}
			}
		}
	}

cout << "Event, # charged = " << locEventLoop->GetJEvent().GetEventNumber() << ", " << _data.size() << endl;
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DChargedTrackHypothesis_factory_Combo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTrackHypothesis_factory_Combo::fini(void)
{
	return NOERROR;
}


