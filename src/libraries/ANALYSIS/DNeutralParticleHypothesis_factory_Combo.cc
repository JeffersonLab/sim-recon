// $Id$
//
//    File: DNeutralParticleHypothesis_factory_Combo.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DNeutralParticleHypothesis_factory_Combo.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses); //make sure that brun() is called for the default factory!!!
	dNeutralParticleHypothesisFactory = static_cast<DNeutralParticleHypothesis_factory*>(locEventLoop->GetFactory("DNeutralParticleHypothesis"));
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DParticleComboBlueprint*> locParticleComboBlueprints;
	locEventLoop->Get(locParticleComboBlueprints);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Combo");

	map<const DEventRFBunch*, deque<pair<const DNeutralShower*, Particle_t> > > locCreatedParticleMap; //don't create if already done!

	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	{
		const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];

		//select the corresponding rf bunch
		const DEventRFBunch* locEventRFBunch = locEventRFBunches[0];
/*
		const DEventRFBunch* locEventRFBunch = NULL;
		for(size_t loc_j = 0; loc_j < locEventRFBunches.size(); ++loc_j)
		{
			vector<const DParticleComboBlueprint*> locParticleComboBlueprints_Bunch;
			locEventRFBunches[loc_j]->Get(locParticleComboBlueprints_Bunch);
			bool locMatchFoundFlag = false;
			for(size_t loc_k = 0; loc_k < locParticleComboBlueprints_Bunch.size(); ++loc_k)
			{
				if(locParticleComboBlueprints_Bunch[loc_k] != locParticleComboBlueprint)
					continue;
				locEventRFBunch = locEventRFBunches[loc_j];
				locMatchFoundFlag = true;
				break;
			}
			if(locMatchFoundFlag)
				break;
		}
*/
		if(locCreatedParticleMap.find(locEventRFBunch) == locCreatedParticleMap.end())
			locCreatedParticleMap[locEventRFBunch] = deque<pair<const DNeutralShower*, Particle_t> >();

		//find the combo showers
		for(size_t loc_j = 0; loc_j < locParticleComboBlueprint->Get_NumParticleComboBlueprintSteps(); ++loc_j)
		{
			const DParticleComboBlueprintStep* locParticleComboBlueprintStep = locParticleComboBlueprint->Get_ParticleComboBlueprintStep(loc_j);

			for(size_t loc_k = 0; loc_k < locParticleComboBlueprintStep->Get_NumFinalParticleSourceObjects(); ++loc_k)
			{
				const JObject* locSourceObject = locParticleComboBlueprintStep->Get_FinalParticle_SourceObject(loc_k);
				if(locSourceObject == NULL)
					continue; //missing or decaying
				const DNeutralShower* locNeutralShower = dynamic_cast<const DNeutralShower*>(locSourceObject);
				if(locNeutralShower == NULL)
					continue; //charged
				Particle_t locPID = locParticleComboBlueprintStep->Get_FinalParticleID(loc_k);

				//see if already created for this pid/track/rf bunch pair
				deque<pair<const DNeutralShower*, Particle_t> > locParticlePairs = locCreatedParticleMap[locEventRFBunch];
				bool locAlreadyCreatedFlag = false;
				for(size_t loc_l = 0; loc_l < locParticlePairs.size(); ++loc_l)
				{
					if((locParticlePairs[loc_l].first != locNeutralShower) || (locParticlePairs[loc_l].second != locPID))
						continue;
					locAlreadyCreatedFlag = true; //e.g., several combos have the same rf bunch
					break;
				}
				if(locAlreadyCreatedFlag)
					continue;

				//create the objects
				DNeutralParticleHypothesis* locNeutralParticleHypothesis = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locNeutralShower, locPID, locEventRFBunch);
				if(locNeutralParticleHypothesis == NULL)
					continue;
				locNeutralParticleHypothesis->AddAssociatedObject(locEventRFBunch);
				_data.push_back(locNeutralParticleHypothesis);
				locCreatedParticleMap[locEventRFBunch].push_back(pair<const DNeutralShower*, Particle_t>(locNeutralShower, locPID));
			}
		}
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::fini(void)
{
	return NOERROR;
}

