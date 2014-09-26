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
jerror_t DNeutralParticleHypothesis_factory_Combo::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses); //make sure that brun() is called for the default factory!!!
	dNeutralParticleHypothesisFactory = static_cast<DNeutralParticleHypothesis_factory*>(locEventLoop->GetFactory("DNeutralParticleHypothesis"));

	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	vector<const DReaction*> locReactions;
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
		locReactions.insert(locReactions.end(), locReactionsSubset.begin(), locReactionsSubset.end());
	}

	//build list of all needed neutral PIDs
	for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
	{
		deque<Particle_t> locDetectedNeutralPIDs;
		locReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedNeutralPIDs, 2, false);
		for(size_t loc_j = 0; loc_j < locDetectedNeutralPIDs.size(); ++loc_j)
			dNeutralPIDs.insert(locDetectedNeutralPIDs[loc_j]);
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticleHypothesis_factory_Combo::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DNeutralParticleHypothesis_factory_Combo::evnt()");
#endif

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Combo");

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers, dShowerSelectionTag.c_str());

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	//assume that: for each RF bunch, must calculate each PID FOM for each shower
		//is true unless: somehow for a given combination of the above, it always fails an invariant mass cut
		//however, it doesn't take much memory if an extra object is created
		//and it's WAY faster than looping over EVERY blueprint and checking what the exceptions are
			//scales much better this way

	set<Particle_t>::iterator locIterator = dNeutralPIDs.begin();
	for(; locIterator != dNeutralPIDs.end(); ++locIterator)
	{
		for(size_t loc_i = 0; loc_i < locEventRFBunches.size(); ++loc_i)
		{
			for(size_t loc_j = 0; loc_j < locNeutralShowers.size(); ++loc_j)
			{
				//create the objects
				DNeutralParticleHypothesis* locNeutralParticleHypothesis = dNeutralParticleHypothesisFactory->Create_DNeutralParticleHypothesis(locNeutralShowers[loc_j], *locIterator, locEventRFBunches[loc_i], locVertex);
				if(locNeutralParticleHypothesis == NULL)
					continue;
				locNeutralParticleHypothesis->AddAssociatedObject(locEventRFBunches[loc_i]);
				_data.push_back(locNeutralParticleHypothesis);
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

