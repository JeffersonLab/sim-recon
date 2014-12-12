// $Id$
//
//    File: DEventRFBunch_factory_Combo.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt
//

#ifdef VTRACE
#include "vt_user.h"
#endif

#include "DEventRFBunch_factory_Combo.h"

//------------------
// init
//------------------
jerror_t DEventRFBunch_factory_Combo::init(void)
{
	dRFBunchFrequency = 2.004;
	dShowerSelectionTag = "PreSelect";
	dTrackSelectionTag = "PreSelect";
	dMinThrownMatchFOM = 5.73303E-7;
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventRFBunch_factory_Combo::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	gPARMS->SetDefaultParameter("COMBO:TRACK_SELECT_TAG", dTrackSelectionTag);
	gPARMS->SetDefaultParameter("COMBO:SHOWER_SELECT_TAG", dShowerSelectionTag);

	DApplication *locApplication = dynamic_cast<DApplication*> (locEventLoop->GetJApplication());
	DGeometry *locGeometry = locApplication ? locApplication->GetDGeometry(runnumber):NULL;
	locGeometry->GetTargetZ(dTargetCenterZ);

	locEventLoop->GetSingle(dParticleID);

	dEventRFBunchFactory = static_cast<DEventRFBunch_factory*>(locEventLoop->GetFactory("DEventRFBunch"));

	// Get DReactions:
	// Get list of factories and find all the ones producing
	// DReaction objects. (A simpler way to do this would be to
	// just use locEventLoop->Get(...), but then only one plugin could
	// be used at a time.)
	vector<JFactory_base*> locFactories = locEventLoop->GetFactories();
	vector<const DReaction*> locReactions;
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

	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	string locHistName, locHistTitle;
	TH1I* loc1DHist;

	//Create Diagnostic Histograms
	japp->RootWriteLock();
	{
		string locOutputFileName = "hd_root.root";
		if(gPARMS->Exists("OUTPUT_FILENAME"))
			gPARMS->GetParameter("OUTPUT_FILENAME", locOutputFileName);
		TFile* locFile = (TFile*)gROOT->FindObject(locOutputFileName.c_str());
		if(locFile == NULL)
			return NOERROR;
		locFile->cd("");

		for(size_t loc_i = 0; loc_i < locReactions.size(); ++loc_i)
		{
			const DReaction* locReaction = locReactions[loc_i];

			//get to the correct directory
			string locReactionName = locReaction->Get_ReactionName();
			string locDirName = locReactionName;
			string locDirTitle = locReactionName;

			//action directory
			locFile->cd();
			TDirectoryFile* locDirectoryFile = static_cast<TDirectoryFile*>(locFile->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
			locDirectoryFile->cd();

			//pre-combo directory
			locDirName = "Hist_RFSelection";
			locDirectoryFile = static_cast<TDirectoryFile*>(gDirectory->GetDirectory(locDirName.c_str()));
			if(locDirectoryFile == NULL)
				locDirectoryFile = new TDirectoryFile(locDirName.c_str(), locDirTitle.c_str());
			locDirectoryFile->cd();

			// RFTime
			locHistName = "RFTime";
			loc1DHist = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
			if(loc1DHist == NULL)
			{
				locHistTitle = locReactionName + string(";Combo-Selected RF Time (ns)");
				loc1DHist = new TH1I(locHistName.c_str(), locHistTitle.c_str(), 420, -21.0, 21.0);
			}
			dHistMap_RFTime[locReaction] = loc1DHist;

			if(!locMCThrowns.empty())
			{
				// DeltaRFTime
				locHistName = "DeltaRFTime";
				loc1DHist = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				if(loc1DHist == NULL)
				{
					locHistTitle = locReactionName + string(";RF #Deltat_{Selected - True} (ns)");
					loc1DHist = new TH1I(locHistName.c_str(), locHistTitle.c_str(), 220, -11.0, 11.0);
				}
				dHistMap_DeltaRFTime[locReaction] = loc1DHist;

				// DeltaRFTime_TruePID
				locHistName = "DeltaRFTime_TruePID";
				loc1DHist = static_cast<TH1I*>(gDirectory->Get(locHistName.c_str()));
				if(loc1DHist == NULL)
				{
					locHistTitle = locReactionName + string(", True PID;RF #Deltat_{Selected - True} (ns)");
					loc1DHist = new TH1I(locHistName.c_str(), locHistTitle.c_str(), 220, -11.0, 11.0);
				}
				dHistMap_DeltaRFTime_TruePID[locReaction] = loc1DHist;
			}
		}
	}
	japp->RootUnLock(); //unlock

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory_Combo::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
#ifdef VTRACE
	VT_TRACER("DEventRFBunch_factory_Combo::evnt()");
#endif

	//use all particles in the combo to find the vertex
		//ignore detached vertices
		//selection is done by propagating all times to the target center, and then voting
	//for neutrals, use DVertex to calculate times
		//if its the good combo: has good tracks, will be good vertex
		//vertex should be bad only if on the wrong combo anyway, or if all tracks are junk

 	vector<const DParticleComboBlueprint*> locParticleComboBlueprints;
	locEventLoop->Get(locParticleComboBlueprints);

	const DEventRFBunch* locEventRFBunch = NULL;
	locEventLoop->GetSingle(locEventRFBunch);

	double locRFTime, locRFVariance;
	if(locEventRFBunch->dTime == locEventRFBunch->dTime)
	{
		locRFTime = locEventRFBunch->dTime;
		locRFVariance = locEventRFBunch->dTimeVariance;
	}
	else if(!dEventRFBunchFactory->Get_RFTimeGuess(locEventLoop, locRFTime, locRFVariance))
	{
		//no good RF time, set to NaN for all combos
		DEventRFBunch* locNewEventRFBunch = new DEventRFBunch(*locEventRFBunch);
		for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
		{
			const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];
			locNewEventRFBunch->AddAssociatedObject(locParticleComboBlueprint);
			locNewEventRFBunch->AddAssociatedObject(locParticleComboBlueprint->Get_Reaction());
		}
		_data.push_back(locNewEventRFBunch);
		return NOERROR;
	}

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector, "Combo");

 	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers, dShowerSelectionTag.c_str());

	vector<const DEventRFBunch*> locThrownEventRFBunches;
	locEventLoop->Get(locThrownEventRFBunches, "Thrown");
	const DEventRFBunch* locThrownEventRFBunch = locThrownEventRFBunches.empty() ? NULL : locThrownEventRFBunches[0];

	vector<const DMCThrownMatching*> locMCThrownMatchingVector;
	locEventLoop->Get(locMCThrownMatchingVector);
	const DMCThrownMatching* locMCThrownMatching = locMCThrownMatchingVector.empty() ? NULL : locMCThrownMatchingVector[0];

 	vector<const DChargedTrack*> locChargedTracks;
	locEventLoop->Get(locChargedTracks, dTrackSelectionTag.c_str());

	const DVertex* locVertex = NULL;
	locEventLoop->GetSingle(locVertex);

	map<int, DEventRFBunch*> locComboRFBunchMap;

	//pre-sort time-based tracks
	map<pair<const DChargedTrack*, Particle_t>, const DTrackTimeBased*> locTimeBasedSourceMap;
	for(size_t loc_l = 0; loc_l < locTrackTimeBasedVector.size(); ++loc_l)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_l];
		const DChargedTrack* locChargedTrack = NULL;
		locTrackTimeBased->GetSingleT(locChargedTrack);
		locTimeBasedSourceMap[pair<const DChargedTrack*, Particle_t>(locChargedTrack, locTrackTimeBased->PID())] = locTrackTimeBased;
	}

	//pre-calculate propagated times: charged
	map<const DChargedTrackHypothesis*, double> dPropagatedStartTimes_Charged;
	for(size_t loc_i = 0; loc_i < locChargedTracks.size(); ++loc_i)
	{
		const vector<const DChargedTrackHypothesis*>& locChargedTrackHypotheses = locChargedTracks[loc_i]->dChargedTrackHypotheses;
		for(size_t loc_j = 0; loc_j < locChargedTrackHypotheses.size(); ++loc_j)
		{
			double locPropagatedTime = locChargedTrackHypotheses[loc_j]->time() + (dTargetCenterZ - locChargedTrackHypotheses[loc_j]->z())/29.9792458;
			dPropagatedStartTimes_Charged[locChargedTrackHypotheses[loc_j]] = locPropagatedTime;
		}
	}

	//pre-calculate propagated times: time-based
	map<const DTrackTimeBased*, double> dPropagatedStartTimes_TimeBased;
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		double locStartTime = 0.0;
		if(!Get_StartTime(locEventLoop, locTrackTimeBasedVector[loc_i], locStartTime))
			continue;
		double locPropagatedTime = locStartTime + (dTargetCenterZ - locTrackTimeBasedVector[loc_i]->z())/29.9792458;
		dPropagatedStartTimes_TimeBased[locTrackTimeBasedVector[loc_i]] = locPropagatedTime;
	}

	//pre-calculate propagated times: neutrals (ignore non-photons)
	map<const DNeutralShower*, double> dPropagatedStartTimes_Neutral;
	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
	{
		double locStartTime = Calc_StartTime(locNeutralShowers[loc_i], locVertex);
		double locPropagatedTime = locStartTime + (dTargetCenterZ - locVertex->dSpacetimeVertex.Z())/29.9792458;
		dPropagatedStartTimes_Neutral[locNeutralShowers[loc_i]] = locPropagatedTime;
	}

	for(size_t loc_i = 0; loc_i < locParticleComboBlueprints.size(); ++loc_i)
	{
		const DParticleComboBlueprint* locParticleComboBlueprint = locParticleComboBlueprints[loc_i];
		vector<double> locPropagatedTimes;

		//Charged Tracks
		deque<pair<const DChargedTrack*, Particle_t> > locChargedTracks;
		locParticleComboBlueprint->Get_DetectedChargedTrackSourceObjects(locChargedTracks);
		for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j)
		{
			const DChargedTrack* locChargedTrack = locChargedTracks[loc_j].first;
			Particle_t locPID = locChargedTracks[loc_j].second;
			const DChargedTrackHypothesis* locChargedTrackHypothesis = locChargedTrack->Get_Hypothesis(locPID);
			if(locChargedTrackHypothesis != NULL)
			{
				locPropagatedTimes.push_back(dPropagatedStartTimes_Charged[locChargedTrackHypothesis]);
				continue;
			}

			//get from time-based
			pair<const DChargedTrack*, Particle_t> locTrackPair(locChargedTrack, locPID);
			map<pair<const DChargedTrack*, Particle_t>, const DTrackTimeBased*>::iterator locIterator = locTimeBasedSourceMap.find(locTrackPair);
			if(locIterator == locTimeBasedSourceMap.end())
				continue; //bad track

			const DTrackTimeBased* locTrackTimeBased = locIterator->second;
			locPropagatedTimes.push_back(dPropagatedStartTimes_TimeBased[locTrackTimeBased]);
		}

		//Neutrals
		deque<pair<const DNeutralShower*, Particle_t> > locNeutralShowers;
		locParticleComboBlueprint->Get_DetectedNeutralShowerSourceObjects(locNeutralShowers);
		for(size_t loc_j = 0; loc_j < locNeutralShowers.size(); ++loc_j)
		{
			Particle_t locPID = locNeutralShowers[loc_j].second;
			if(locPID != Gamma)
				continue; //other neutrals (e.g. neutron) can't be used to pick the time: their momentum is defined by the time

			const DNeutralShower* locNeutralShower = locNeutralShowers[loc_j].first;
			locPropagatedTimes.push_back(dPropagatedStartTimes_Neutral[locNeutralShower]);
		}

		if(locPropagatedTimes.empty())
		{
			//no timing information somehow: use the pre-existing value
			DEventRFBunch* locNewEventRFBunch = new DEventRFBunch(*locEventRFBunch);
			locNewEventRFBunch->AddAssociatedObject(locParticleComboBlueprint);
			locNewEventRFBunch->AddAssociatedObject(locParticleComboBlueprint->Get_Reaction());
			_data.push_back(locNewEventRFBunch);
			continue;
		}

		// Find # RF Bunch Shifts
		int locNumBunchShifts = Find_BestRFBunchShift(locRFTime, locPropagatedTimes);
		double locNewRFTime = locRFTime + (double)(locNumBunchShifts)*dRFBunchFrequency;

		//Hist
		bool locIsAllTruePID = false;
		if(locThrownEventRFBunch != NULL)
			locIsAllTruePID = Is_AllTruePID(locMCThrownMatching, locParticleComboBlueprint);
		const DReaction* locReaction = locParticleComboBlueprint->Get_Reaction();
		japp->RootWriteLock();
		{
			dHistMap_RFTime[locReaction]->Fill(locNewRFTime);
			if(locThrownEventRFBunch != NULL)
			{
				double locDeltaT = locNewRFTime - locThrownEventRFBunch->dTime;
				dHistMap_DeltaRFTime[locReaction]->Fill(locDeltaT); //diff between selected & true RF times (all combos)
				if(locIsAllTruePID)
					dHistMap_DeltaRFTime_TruePID[locReaction]->Fill(locDeltaT); //diff between selected & true RF times (all combos)
			}
		}
		japp->RootUnLock(); //unlock

		// Create new RF Bunch if doesn't already exist
		if(locComboRFBunchMap.find(locNumBunchShifts) != locComboRFBunchMap.end()) //already created, don't recreate identical object!
		{
			locComboRFBunchMap[locNumBunchShifts]->AddAssociatedObject(locParticleComboBlueprint);
			locComboRFBunchMap[locNumBunchShifts]->AddAssociatedObject(locParticleComboBlueprint->Get_Reaction());
		}
		else
		{
			DEventRFBunch* locNewEventRFBunch = new DEventRFBunch();
			locNewEventRFBunch->dTime = locNewRFTime;
			locNewEventRFBunch->dTimeVariance = locRFVariance;
			locNewEventRFBunch->dNumParticlesVotedForThisTime = locPropagatedTimes.size();
			locNewEventRFBunch->AddAssociatedObject(locParticleComboBlueprint);
			locNewEventRFBunch->AddAssociatedObject(locParticleComboBlueprint->Get_Reaction());
			_data.push_back(locNewEventRFBunch);
			locComboRFBunchMap[locNumBunchShifts] = locNewEventRFBunch;
		}
	}

	return NOERROR;
}

bool DEventRFBunch_factory_Combo::Get_StartTime(JEventLoop* locEventLoop, const DTrackTimeBased* locTrackTimeBased, double& locStartTime)
{
	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	// Use time-based tracking time as initial guess
	locStartTime = 0.0;

	//BCAL
	DShowerMatchParams locBCALShowerMatchParams;
	if(dParticleID->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
	{
		const DBCALShower* locBCALShower = dynamic_cast<const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
		locStartTime = locBCALShower->t - locBCALShowerMatchParams.dFlightTime;
		return true;
	}

	//TOF
	DTOFHitMatchParams locTOFHitMatchParams;
	if(dParticleID->Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
	{
		const DTOFPoint* locTOFPoint = locTOFHitMatchParams.dTOFPoint;
		locStartTime = locTOFPoint->t - locTOFHitMatchParams.dFlightTime;
		return true;
	}

	//SC
	DSCHitMatchParams locSCHitMatchParams;
	if(dParticleID->Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams))
	{
		locStartTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime;
		return true;
	}

	//FCAL
	DShowerMatchParams locFCALShowerMatchParams;
	if(dParticleID->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams))
	{
		const DFCALShower* locFCALShower = dynamic_cast<const DFCALShower*>(locFCALShowerMatchParams.dShowerObject);
		locStartTime = locFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime;
		return true;
	}

	return false;
}

double DEventRFBunch_factory_Combo::Calc_StartTime(const DNeutralShower* locNeutralShower, const DVertex* locVertex)
{
	//doesn't work for neutrons!!

	DVector3 locHitPoint = locNeutralShower->dSpacetimeVertex.Vect();
	DVector3 locPath = locHitPoint - locVertex->dSpacetimeVertex.Vect();
	double locPathLength = locPath.Mag();

	double locFlightTime = locPathLength/29.9792458;
	double locHitTime = locNeutralShower->dSpacetimeVertex.T();
	return locHitTime - locFlightTime;
}


int DEventRFBunch_factory_Combo::Find_BestRFBunchShift(double locRFHitTime, const vector<double>& locTimes)
{
	//then find the #beam buckets the RF time needs to shift to match it
	map<int, unsigned int> locNumRFBucketsShiftedMap; //key is # RF buckets the RF time is shifted to match the track, value is the # tracks at that shift
	int locBestRFBunchShift = 0;
	for(unsigned int loc_i = 0; loc_i < locTimes.size(); ++loc_i)
	{
		double locPropagatedTrackTime = locTimes[loc_i];
		//do manually: tricky to convert int to float...
		int locNumRFBucketsShifted = 0;
		double locTempRFHitTime = locRFHitTime;
		while((locTempRFHitTime - locPropagatedTrackTime) > (0.5*dRFBunchFrequency))
		{
			locTempRFHitTime -= dRFBunchFrequency;
			--locNumRFBucketsShifted;
		}
		while((locTempRFHitTime - locPropagatedTrackTime) < (-0.5*dRFBunchFrequency))
		{
			locTempRFHitTime += dRFBunchFrequency;
			++locNumRFBucketsShifted;
		}

		if(locNumRFBucketsShiftedMap.find(locNumRFBucketsShifted) == locNumRFBucketsShiftedMap.end())
			locNumRFBucketsShiftedMap[locNumRFBucketsShifted] = 1;
		else
			++(locNumRFBucketsShiftedMap[locNumRFBucketsShifted]);

		if(locNumRFBucketsShifted == locBestRFBunchShift)
			continue;

		unsigned int locBestNumTracks = locNumRFBucketsShiftedMap[locBestRFBunchShift];
		unsigned int locNumTracks = locNumRFBucketsShiftedMap[locNumRFBucketsShifted];
		if(locNumTracks > locBestNumTracks)
			locBestRFBunchShift = locNumRFBucketsShifted;
	}
	return locBestRFBunchShift;
}

bool DEventRFBunch_factory_Combo::Is_AllTruePID(const DMCThrownMatching* locMCThrownMatching, const DParticleComboBlueprint* locParticleComboBlueprint)
{
	//Charged
	deque<pair<const DChargedTrack*, Particle_t> > locChargedTracks;
	locParticleComboBlueprint->Get_DetectedChargedTrackSourceObjects(locChargedTracks);
	for(size_t loc_j = 0; loc_j < locChargedTracks.size(); ++loc_j)
	{
		const DChargedTrack* locChargedTrack = locChargedTracks[loc_j].first;
		Particle_t locPID = locChargedTracks[loc_j].second;

		double locMatchFOM = 0.0;
		const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locChargedTrack, locMatchFOM);
		if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
			return false;
		if(((Particle_t)locMCThrown->type) != locPID)
			return false;
	}

	//Neutrals
	deque<pair<const DNeutralShower*, Particle_t> > locNeutralShowers;
	locParticleComboBlueprint->Get_DetectedNeutralShowerSourceObjects(locNeutralShowers);
	for(size_t loc_j = 0; loc_j < locNeutralShowers.size(); ++loc_j)
	{
		Particle_t locPID = locNeutralShowers[loc_j].second;
		if(locPID != Gamma)
			continue; //other neutrals (e.g. neutron) can't be used to pick the time: their momentum is defined by the time

		const DNeutralShower* locNeutralShower = locNeutralShowers[loc_j].first;

		double locMatchFOM = 0.0;
		const DMCThrown* locMCThrown = locMCThrownMatching->Get_MatchingMCThrown(locNeutralShower, locMatchFOM);
		if((locMCThrown == NULL) || (locMatchFOM < dMinThrownMatchFOM))
			return false;
		if(((Particle_t)locMCThrown->type) != locPID)
			return false;
	}

	return true;
}

//------------------
// erun
//------------------
jerror_t DEventRFBunch_factory_Combo::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventRFBunch_factory_Combo::fini(void)
{
	return NOERROR;
}


