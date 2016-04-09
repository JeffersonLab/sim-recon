// $Id$
//
//    File: DNeutralParticleHypothesis_factory_Combo.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#ifdef VTRACE
#include "vt_user.h"
#endif

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
	dShowerSelectionTag = "PreSelect";

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses); //make sure that brun() is called for the default factory!!!
	dNeutralParticleHypothesisFactory = static_cast<DNeutralParticleHypothesis_factory*>(locEventLoop->GetFactory("DNeutralParticleHypothesis"));

	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	dReactions.clear();
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	for(size_t loc_i = 0; loc_i < locFactories.size(); ++loc_i)
	{
		JFactory<DReaction>* locFactory = dynamic_cast<JFactory<DReaction>* >(locFactories[loc_i]);
		if(locFactory == NULL)
			continue;
		if(string(locFactory->Tag()) == "Thrown")
			continue;
		// Found a factory producing DReactions. The reaction objects are
		// produced at the init stage and are persistent through all event
		// processing so we can grab the list here and append it to our
		// overall list.
		vector<const DReaction*> locReactionsSubset;
		locFactory->Get(locReactionsSubset);
		dReactions.insert(dReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}

	//build list of all needed neutral PIDs
	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		deque<Particle_t> locDetectedNeutralPIDs;
		dReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedNeutralPIDs, 2, false);
		for(size_t loc_j = 0; loc_j < locDetectedNeutralPIDs.size(); ++loc_j)
			dNeutralPIDs[locDetectedNeutralPIDs[loc_j]].push_back(dReactions[loc_i]);
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DNeutralParticleHypothesis_factory_Combo::evnt()");
#endif

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Combo");

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles, dShowerSelectionTag.c_str());

	vector<const DNeutralParticleHypothesis*> locOrigNeutralParticleHypotheses;
	locEventLoop->Get(locOrigNeutralParticleHypotheses);

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	//pre-sort rf bunches
	//for each PID, what are the rf bunches that I need? determined via links to reactions
	map<const DReaction*, set<const DEventRFBunch*> > locRFBunchReactionMap;
	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locEventRFBunches.size(); ++loc_j)
		{
			if(locEventRFBunches[loc_j]->IsAssociated(dReactions[loc_i]))
				locRFBunchReactionMap[dReactions[loc_i]].insert(locEventRFBunches[loc_j]);
		}
	}
	map<Particle_t, set<const DEventRFBunch*> > locRFBunchPIDMap;
	Build_RFPIDMap(dNeutralPIDs, locRFBunchReactionMap, locRFBunchPIDMap);

	//assume that: for each RF bunch, must calculate each PID FOM for each shower
		//is true unless: somehow for a given combination of the above, it always fails an invariant mass cut
		//however, it doesn't take much memory if an extra object is created
		//and it's WAY faster than looping over EVERY blueprint and checking what the exceptions are
			//scales much better this way

	for(size_t loc_j = 0; loc_j < locNeutralParticles.size(); ++loc_j)
	{
		map<Particle_t, vector<const DReaction*> >::const_iterator locPIDIterator = dNeutralPIDs.begin();
		for(; locPIDIterator != dNeutralPIDs.end(); ++locPIDIterator)
		{
			Particle_t locPID = locPIDIterator->first;
			vector<const DReaction*> locReactions = locPIDIterator->second;

			//loop over rf bunches for this PID
			set<const DEventRFBunch*> locRFBunches = locRFBunchPIDMap[locPID];
			set<const DEventRFBunch*>::const_iterator locRFIterator = locRFBunches.begin();
			for(; locRFIterator != locRFBunches.end(); ++locRFIterator)
			{
				//create
				const DNeutralShower* locNeutralShower = locNeutralParticles[loc_j]->dNeutralShower;
				DNeutralParticleHypothesis* locNeutralParticleHypothesis = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locNeutralShower, locPID, *locRFIterator, locVertex);
				if(locNeutralParticleHypothesis == NULL)
					continue;

				locNeutralParticleHypothesis->AddAssociatedObject(locNeutralParticles[loc_j]);
				locNeutralParticleHypothesis->AddAssociatedObject(*locRFIterator);
				for(size_t loc_k = 0; loc_k < locReactions.size(); ++loc_k)
					locNeutralParticleHypothesis->AddAssociatedObject(locReactions[loc_k]);

				const DNeutralParticleHypothesis* locOrigNeutralParticleHypothesis = locNeutralParticles[loc_j]->Get_Hypothesis(locPID);
				if(locOrigNeutralParticleHypothesis != NULL)
					locNeutralParticleHypothesis->AddAssociatedObject(locOrigNeutralParticleHypothesis);

				_data.push_back(locNeutralParticleHypothesis);
			}
		}
	}

	return NOERROR;
}

void DNeutralParticleHypothesis_factory_Combo::Build_RFPIDMap(map<Particle_t, vector<const DReaction*> >& locPIDMap, map<const DReaction*, set<const DEventRFBunch*> >& locRFBunchReactionMap, map<Particle_t, set<const DEventRFBunch*> >& locRFBunchPIDMap)
{
	map<Particle_t, vector<const DReaction*> >::const_iterator locIterator = locPIDMap.begin();
	for(; locIterator != locPIDMap.end(); ++locIterator)
	{
		Particle_t locPID = locIterator->first;
		vector<const DReaction*> locReactions = locIterator->second;
		for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
		{
			set<const DEventRFBunch*>& locRFBunches = locRFBunchReactionMap[locReactions[loc_i]];
			locRFBunchPIDMap[locPID].insert(locRFBunches.begin(), locRFBunches.end());
		}
	}
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

