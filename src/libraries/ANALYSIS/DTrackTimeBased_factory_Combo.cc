// $Id$
//
//    File: DTrackTimeBased_factory_Combo.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DTrackTimeBased_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DTrackTimeBased_factory_Combo::init(void)
{
	MAX_dReferenceTrajectoryPoolSize = 10;

	dMinTrackingFOM = pair<bool, double>(false, -1.0);
	dHasDetectorMatchFlag = pair<bool, bool>(false, false);
	dMinProtonMomentum = pair<bool, double>(false, -1.0);

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

	if(gPARMS->Exists("COMBO:MIN_PROTON_MOMENTUM"))
	{
		dMinProtonMomentum.first = true;
		gPARMS->GetParameter("COMBO:MIN_PROTON_MOMENTUM", dMinProtonMomentum.second);
	}

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
jerror_t DTrackTimeBased_factory_Combo::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DTrackTimeBased_factory_Combo::evnt()");
#endif

	// delete pool sizes if too large, preventing memory-leakage-like behavor.
	if(dReferenceTrajectoryPool_All.size() > MAX_dReferenceTrajectoryPoolSize)
	{
		for(size_t loc_i = MAX_dReferenceTrajectoryPoolSize; loc_i < dReferenceTrajectoryPool_All.size(); ++loc_i)
			delete dReferenceTrajectoryPool_All[loc_i];
		dReferenceTrajectoryPool_All.resize(MAX_dReferenceTrajectoryPoolSize);
	}
	dReferenceTrajectoryPool_Available = dReferenceTrajectoryPool_All;

	locEventLoop->GetSingle(dDetectorMatches);

 	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks);

	//assume that: for each PID, a hypothesis must exist for each track that contains that charge
		//definition of "good" is different from DReaction to reaction

	//is true unless: for a given PID, a track always fails an invariant mass cut
		//certainly possible, maybe even likely (e.g. testing a fast pion as a proton)
		//however, it doesn't take much memory if an extra object is created, unless you have to swim
			//and even if you swim, it's probably only a few cases (since swimming should have been performed earlier)
		//and it's WAY faster than looping over EVERY blueprint and checking what the exceptions are
			//scales much better this way

	for(size_t loc_i = 0; loc_i < dReactions.size(); ++loc_i)
	{
		for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j)
		{
			//q+
			if(locChargedTracks[loc_j]->Contains_Charge(1))
			{
				set<Particle_t>& locPIDs = dPositivelyChargedPIDs[dReactions[loc_i]];
				Create_PIDsAsNeeded(dReactions[loc_i], locChargedTracks[loc_j], locPIDs);
			}
			//q-
			if(locChargedTracks[loc_j]->Contains_Charge(-1))
			{
				set<Particle_t>& locPIDs = dNegativelyChargedPIDs[dReactions[loc_i]];
				Create_PIDsAsNeeded(dReactions[loc_i], locChargedTracks[loc_j], locPIDs);
			}
		}
	}

	return NOERROR;
}

void DTrackTimeBased_factory_Combo::Create_PIDsAsNeeded(const DReaction* locReaction, const DChargedTrack* locChargedTrack, set<Particle_t>& locPIDs)
{
	set<Particle_t>::iterator locIterator = locPIDs.begin();
	for(; locIterator != locPIDs.end(); ++locIterator)
	{
		Particle_t locPID = *locIterator;
		if(locChargedTrack->Get_Hypothesis(locPID) != NULL)
			continue; //already exists

		bool locReSwimFlag = false;
		const DChargedTrackHypothesis* locChargedTrackHypothesis = Get_ChargedHypothesisToUse(locChargedTrack, locPID, locReSwimFlag);

		//create it if it's a good track
		if((!locReSwimFlag) && (!Cut_HasDetectorMatch(locReaction, locChargedTrackHypothesis)))
			continue;
		if(!Cut_TrackingFOM(locReaction, locChargedTrackHypothesis))
			continue;

		DTrackTimeBased* locTrackTimeBased = Convert_ChargedTrack(locChargedTrackHypothesis, locPID, locReSwimFlag);
		if(locTrackTimeBased == NULL)
			continue;
		locTrackTimeBased->AddAssociatedObject(locChargedTrack);
		_data.push_back(locTrackTimeBased);
	}
}

const DChargedTrackHypothesis* DTrackTimeBased_factory_Combo::Get_ChargedHypothesisToUse(const DChargedTrack* locChargedTrack, Particle_t locDesiredPID, bool& locReSwimFlag)
{
	//pid not found for this track: loop over other possible pids
	for(size_t loc_i = 0; loc_i < dParticleIDsToTry[locDesiredPID].size(); ++loc_i)
	{
		const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(dParticleIDsToTry[locDesiredPID][loc_i].first);
		if(locChargedTrackHypothesis == NULL)
			continue;
		locReSwimFlag = dParticleIDsToTry[locDesiredPID][loc_i].second;
		return locChargedTrackHypothesis;
	}

	//still none found, take the one with the best FOM
	locReSwimFlag = true;
	return locChargedTrack->Get_BestFOM();
}

bool DTrackTimeBased_factory_Combo::Cut_HasDetectorMatch(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locHasDetectorMatchFlag = dHasDetectorMatchFlag.first ? dHasDetectorMatchFlag : locReaction->Get_HasDetectorMatchFlag();
	if((!locHasDetectorMatchFlag.first) || (!locHasDetectorMatchFlag.second))
		return true;

	const DTrackTimeBased* locTrackTimeBased = NULL;
	locChargedTrackHypothesis->GetSingle(locTrackTimeBased);
	return dDetectorMatches->Get_IsMatchedToHit(locTrackTimeBased);
}

bool DTrackTimeBased_factory_Combo::Cut_TrackingFOM(const DReaction* locReaction, const DChargedTrackHypothesis* locChargedTrackHypothesis) const
{
	pair<bool, double> locMinTrackingFOM = dMinTrackingFOM.first ? dMinTrackingFOM : locReaction->Get_MinTrackingFOM();
	if(!locMinTrackingFOM.first)
		return true;
	double locFOM = TMath::Prob(locChargedTrackHypothesis->dChiSq_Track, locChargedTrackHypothesis->dNDF_Track);
	return ((locChargedTrackHypothesis->dNDF_Track == 0) ? true : (locFOM >= locMinTrackingFOM.second));
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
#ifdef VTRACE
	VT_TRACER("DTrackTimeBased_factory_Combo::Convert_ChargedTrack()");
#endif

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


