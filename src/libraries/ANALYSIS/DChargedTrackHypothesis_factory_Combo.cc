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
	dMinTrackingFOM = pair<bool, double>(false, -1.0);
	dHasDetectorMatchFlag = pair<bool, bool>(false, false);
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

	if(gPARMS->Exists("COMBO:MIN_TRACKING_FOM"))
	{
		dMinTrackingFOM.first = true;
		gPARMS->GetParameter("COMBO:MIN_TRACKING_FOM", dMinTrackingFOM.second);
	}

	if(gPARMS->Exists("COMBO:HAS_DETECTOR_MATCH_FLAG"))
	{
		dHasDetectorMatchFlag.first = true;
		gPARMS->GetParameter("COMBO:HAS_DETECTOR_MATCH_FLAG", dHasDetectorMatchFlag.second);
	}

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
			dPositivelyChargedPIDs[dReactions[loc_i]].insert(locDetectedPIDs[loc_j]);

		locDetectedPIDs.clear();
		dReactions[loc_i]->Get_DetectedFinalPIDs(locDetectedPIDs, 4, false); //q+
		for(size_t loc_j = 0; loc_j < locDetectedPIDs.size(); ++loc_j)
			dNegativelyChargedPIDs[dReactions[loc_i]].insert(locDetectedPIDs[loc_j]);
	}

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

	//assume that: for each RF bunch, must calculate each PID FOM for each "good" charged track (& time-based "combo" track)
		//definition of "good" is different from DReaction to reaction
		//very few rf bunches: don't worry about only selecting bunches corresponding to the given DReaction: just loop over all of them 

	//is true unless: for a given PID, a track always fails an invariant mass cut
		//certainly possible, maybe even likely (e.g. testing a fast pion as a proton)
		//however, it doesn't take much memory if an extra object is created
		//and it's WAY faster than looping over EVERY blueprint and checking what the exceptions are
			//scales much better this way

	dCreatedParticleMap.clear();
	dTimeBasedSourceMap.clear();

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector, "Combo");

 	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	vector<const DEventRFBunch*> locEventRFBunches;
	locEventLoop->Get(locEventRFBunches, "Combo");

	locEventLoop->GetSingle(dDetectorMatches, "Combo");

	//pre-sort time-based tracks
	for(size_t loc_l = 0; loc_l < locTrackTimeBasedVector.size(); ++loc_l)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_l];
		const DChargedTrack* locChargedTrack = NULL;
		locTrackTimeBased->GetSingleT(locChargedTrack);
		dTimeBasedSourceMap[pair<const DChargedTrack*, Particle_t>(locChargedTrack, locTrackTimeBased->PID())] = locTrackTimeBased;
	}

	//do it
	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locEventRFBunches.size(); ++loc_j)
		{
			for(size_t loc_k = 0; loc_k < locChargedTracks.size(); ++loc_k)
			{
				//q+
				if(locChargedTracks[loc_k]->Contains_Charge(1))
				{
					set<Particle_t>& locPIDs = dPositivelyChargedPIDs[dReactions[loc_i]];
					Create_PIDsAsNeeded(locEventLoop, dReactions[loc_i], locEventRFBunches[loc_j], locChargedTracks[loc_k], locPIDs);
				}
				//q-
				if(locChargedTracks[loc_k]->Contains_Charge(-1))
				{
					set<Particle_t>& locPIDs = dNegativelyChargedPIDs[dReactions[loc_i]];
					Create_PIDsAsNeeded(locEventLoop, dReactions[loc_i], locEventRFBunches[loc_j], locChargedTracks[loc_k], locPIDs);
				}
			}
		}
	}

	return NOERROR;
}

void DChargedTrackHypothesis_factory_Combo::Create_PIDsAsNeeded(JEventLoop* locEventLoop, const DReaction* locReaction, const DEventRFBunch* locEventRFBunch, const DChargedTrack* locChargedTrack, set<Particle_t>& locPIDs)
{
	set<Particle_t>::iterator locIterator = locPIDs.begin();
	for(; locIterator != locPIDs.end(); ++locIterator)
	{
		Particle_t locPID = *locIterator;

		//if already created, don't create new object
		set<Particle_t>& locCreatedPIDs = dCreatedParticleMap[locEventRFBunch][locChargedTrack];
		if(locCreatedPIDs.find(locPID) != locCreatedPIDs.end())
			continue; //already created

		//see if DChargedTrackHypothesis with the desired PID was created by the default factory
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
		if(locChargedTrackHypothesis != NULL)
		{
			//it was. create new object it if it's a good track
			if(!Cut_HasDetectorMatch(locReaction, locChargedTrackHypothesis))
				continue;
			if(!Cut_TrackingFOM(locReaction, locChargedTrackHypothesis))
				continue;

			//it is. create new object with same PID (so that is registered with the combo factory, and because rf bunch could be different)
			const DTrackTimeBased* locTrackTimeBased = NULL;
			locChargedTrackHypothesis->GetSingleT(locTrackTimeBased);
			DChargedTrackHypothesis* locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBased, dDetectorMatches, locEventRFBunch);
			locNewChargedTrackHypothesis->AddAssociatedObject(locEventRFBunch);
			locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
			_data.push_back(locNewChargedTrackHypothesis);
			locCreatedPIDs.insert(locPID);
			continue;
		}

		//no DChargedTrackHypothesis with this PID: get track info from DTrackTimeBased "Combo"-tag objects
		pair<const DChargedTrack*, Particle_t> locTrackPair(locChargedTrack, locPID);
		map<pair<const DChargedTrack*, Particle_t>, const DTrackTimeBased*>::iterator locIterator = dTimeBasedSourceMap.find(locTrackPair);
		if(locIterator == dTimeBasedSourceMap.end())
			continue; //bad track
		const DTrackTimeBased* locTrackTimeBased = locIterator->second;
		if(!Cut_TrackingFOM(locReaction, locTrackTimeBased))
			continue; //bad track

		//correct DTrackTimeBased grabbed for this source object: create new DChargedTrackHypothesis object
		DChargedTrackHypothesis* locNewChargedTrackHypothesis = dChargedTrackHypothesisFactory->Create_ChargedTrackHypothesis(locEventLoop, locTrackTimeBased, dDetectorMatches, locEventRFBunch);

		//create it if it's a good track
		if(!Cut_HasDetectorMatch(locReaction, locNewChargedTrackHypothesis))
		{
			delete locNewChargedTrackHypothesis;
			continue; //bad track
		}

		locNewChargedTrackHypothesis->AddAssociatedObject(locEventRFBunch);
		locNewChargedTrackHypothesis->AddAssociatedObject(locChargedTrack);
		_data.push_back(locNewChargedTrackHypothesis);
		locCreatedPIDs.insert(locPID);
	}
}

bool DChargedTrackHypothesis_factory_Combo::Cut_HasDetectorMatch(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locHasDetectorMatchFlag = dHasDetectorMatchFlag.first ? dHasDetectorMatchFlag : locReaction->Get_HasDetectorMatchFlag();
	if((!locHasDetectorMatchFlag.first) || (!locHasDetectorMatchFlag.second))
		return true;

	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
	return dDetectorMatches->Get_IsMatchedToHit(locTrackTimeBased);
}

bool DChargedTrackHypothesis_factory_Combo::Cut_TrackingFOM(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locMinTrackingFOM = dMinTrackingFOM.first ? dMinTrackingFOM : locReaction->Get_MinTrackingFOM();
	if(!locMinTrackingFOM.first)
		return true;
	double locFOM = TMath::Prob(locChargedTrackHypothesis->dChiSq_Track, locChargedTrackHypothesis->dNDF_Track);
	return ((locChargedTrackHypothesis->dNDF_Track == 0) ? true : (locFOM >= locMinTrackingFOM.second));
}

bool DChargedTrackHypothesis_factory_Combo::Cut_TrackingFOM(const DReaction* locReaction, const DTrackTimeBased* locTrackTimeBased) const
{
	pair<bool, double> locMinTrackingFOM = dMinTrackingFOM.first ? dMinTrackingFOM : locReaction->Get_MinTrackingFOM();
	if(!locMinTrackingFOM.first)
		return true;
	return ((locTrackTimeBased->Ndof == 0) ? true : (locTrackTimeBased->FOM >= locMinTrackingFOM.second));
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


