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
	dTrackSelectionTag = "PreSelect";
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrackHypothesis_factory_Combo::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses); //make sure that brun() is called for the default factory!!!
	dChargedTrackHypothesisFactory = static_cast<DChargedTrackHypothesis_factory*>(locEventLoop->GetFactory("DChargedTrackHypothesis"));

	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);

	// Get DReactions:
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	dReactions.clear();
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

	//Get Needed PIDs
	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		//locChargeFlag: 0/1/2/3/4 for all, charged, neutral, q+, q- particles
		deque<Particle_t> locDetectedPIDs;
		dReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedPIDs, 3, false); //q+
		for(size_t loc_j = 0; loc_j < locDetectedPIDs.size(); ++loc_j)
			dPositivelyChargedPIDs[locDetectedPIDs[loc_j]].push_back(dReactions[loc_i]);

		locDetectedPIDs.clear();
		dReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedPIDs, 4, false); //q-
		for(size_t loc_j = 0; loc_j < locDetectedPIDs.size(); ++loc_j)
			dNegativelyChargedPIDs[locDetectedPIDs[loc_j]].push_back(dReactions[loc_i]);
	}

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrackHypothesis_factory_Combo::evnt(jana::JEventLoop* locEventLoop, uint64_t eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DChargedTrackHypothesis_factory_Combo::evnt()");
#endif

	//assume that: for each RF bunch, must calculate each PID FOM for each "good" charged track (& time-based "combo" track)
		//definition of "good" is different from DReaction to reaction
		//very few rf bunches: don't worry about only selecting bunches corresponding to the given DReaction: just loop over all of them 

	//is true unless: for a given PID, a track always fails an invariant mass cut
		//certainly possible, maybe even likely (e.g. testing a fast pion as a proton)
		//however, it doesn't take much memory if an extra object is created
		//and it's WAY faster than looping over EVERY blueprint and checking what the exceptions are
			//scales much better this way

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector, "Combo");

 	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Combo");

	locEventLoop->GetSingle(dDetectorMatches, "Combo");

	//pre-sort time-based tracks
	dTimeBasedSourceMap.clear();
	for(size_t loc_l = 0; loc_l < locTrackTimeBasedVector.size(); ++loc_l)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_l];
		const DChargedTrack* locChargedTrack = NULL;
		locTrackTimeBased->GetSingleT(locChargedTrack);
		dTimeBasedSourceMap[pair<const DChargedTrack*, Particle_t>(locChargedTrack, locTrackTimeBased->PID())] = locTrackTimeBased;
	}

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
	Build_RFPIDMap(dPositivelyChargedPIDs, locRFBunchReactionMap, locRFBunchPIDMap);
	Build_RFPIDMap(dNegativelyChargedPIDs, locRFBunchReactionMap, locRFBunchPIDMap);

	//do it
	for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j)
	{
		if(locChargedTracks[loc_j]->Contains_Charge(1))
			Create_TrackHypos(locEventLoop, locChargedTracks[loc_j], dPositivelyChargedPIDs, locRFBunchPIDMap);
		if(locChargedTracks[loc_j]->Contains_Charge(-1))
			Create_TrackHypos(locEventLoop, locChargedTracks[loc_j], dNegativelyChargedPIDs, locRFBunchPIDMap);
	}

	return NOERROR;
}

void DChargedTrackHypothesis_factory_Combo::Build_RFPIDMap(map<Particle_t, vector<const DReaction*> >& locPIDMap, map<const DReaction*, set<const DEventRFBunch*> >& locRFBunchReactionMap, map<Particle_t, set<const DEventRFBunch*> >& locRFBunchPIDMap)
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

void DChargedTrackHypothesis_factory_Combo::Create_TrackHypos(JEventLoop* locEventLoop, const DChargedTrack* locChargedTrack, map<Particle_t, vector<const DReaction*> >& locPIDMap, map<Particle_t, set<const DEventRFBunch*> >& locRFBunchPIDMap)
{
	map<Particle_t, vector<const DReaction*> >::const_iterator locPIDIterator = locPIDMap.begin();
	for(; locPIDIterator != locPIDMap.end(); ++locPIDIterator)
	{
		Particle_t locPID = locPIDIterator->first;
		vector<const DReaction*> locReactions = locPIDIterator->second;

		//loop over rf bunches for this PID
		set<const DEventRFBunch*> locRFBunches = locRFBunchPIDMap[locPID];
		set<const DEventRFBunch*>::const_iterator locRFIterator = locRFBunches.begin();
		for(; locRFIterator != locRFBunches.end(); ++locRFIterator)
		{
			DChargedTrackHypothesis* locNewHypothesis = Create_TrackHypo(locEventLoop, *locRFIterator, locChargedTrack, locPID);
			if(locNewHypothesis == NULL)
				continue;
			for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
				locNewHypothesis->AddAssociatedObject(locReactions[loc_i]);
		}
	}
}

DChargedTrackHypothesis* DChargedTrackHypothesis_factory_Combo::Create_TrackHypo(JEventLoop* locEventLoop, const DEventRFBunch* locEventRFBunch, const DChargedTrack* locChargedTrack, Particle_t locPID)
{
	//see if DChargedTrackHypothesis with the desired PID was created by the default factory, AND it passed the PreSelect cuts
	const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
	if(locChargedTrackHypothesis != NULL)
	{
		//it was/did. create new object with same PID (so that is registered with the combo factory, and because rf bunch could be different)
		const DTrackTimeBased* locTrackTimeBased = NULL;
		locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);
		DChargedTrackHypothesis* locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBased, dDetectorMatches, locEventRFBunch);

		locNewChargedTrackHypothesis->AddAssociatedObject(locEventRFBunch);
		locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
		locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrackHypothesis);
		_data.push_back(locNewChargedTrackHypothesis);
		return locNewChargedTrackHypothesis;
	}

	//no DChargedTrackHypothesis with this PID:
		//either it failed the PreSelect cuts, or was never reconstructed in the first place
	//see if in DTrackTimeBased "Combo"-tag objects
	pair<const DChargedTrack*, Particle_t> locTrackPair(locChargedTrack, locPID);
	map<pair<const DChargedTrack*, Particle_t>, const DTrackTimeBased*>::iterator locIterator = dTimeBasedSourceMap.find(locTrackPair);
	if(locIterator == dTimeBasedSourceMap.end())
		return NULL; //it failed the PreSelect cuts: bad track
	const DTrackTimeBased* locTrackTimeBased = locIterator->second;

	//correct DTrackTimeBased grabbed for this source object: create new DChargedTrackHypothesis object
	DChargedTrackHypothesis* locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBased, dDetectorMatches, locEventRFBunch);

	locNewChargedTrackHypothesis->AddAssociatedObject(locEventRFBunch);
	locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
	_data.push_back(locNewChargedTrackHypothesis);
	return locNewChargedTrackHypothesis;
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
