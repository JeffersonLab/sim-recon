// $Id$
//
//    File: DEventRFBunch_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DEventRFBunch_factory.h"

using namespace jana;

//------------------
// init
//------------------
jerror_t DEventRFBunch_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventRFBunch_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	dTargetCenter.SetXYZ(0.0, 0.0, 65.0);
	dMinTrackingFOM = 0.001;
	dRFBunchFrequency = 2.004;
	dMinVertexZ = 45.0;
	dMaxVertexZ = 85.0;

	vector<const DParticleID*> locParticleIDVector;
	locEventLoop->Get(locParticleIDVector);
	dPIDAlgorithm = locParticleIDVector[0];

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	// https://halldweb1.jlab.org/wiki/index.php/How_HDGeant_defines_time-zero_for_physics_events

/*
	//COMMENTED UNTIL DRawRFBunch CREATED
	vector<const DRawRFBunch*> locRawRFBunches;
	locEventLoop->Get(locRawRFBunches);
	if(locRawRFBunches.empty())
		return RESOURCE_NOT_AVAILABLE;
	double locRFHitTime = locRawRFBunches[0]->dTime; //get from raw hit when ready
	double locTimeVariance = locRawRFBunches[0]->dTimeVariance;
*/

	double locRFHitTime = 0.0;
	double locTimeVariance = 0.0;

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	//select the best DTrackTimeBased for each track: use best tracking FOM
	map<JObject::oid_t, const DTrackTimeBased*> locBestTrackTimeBasedMap; //lowest tracking chisq/ndf for each candidate id
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		if((locTrackTimeBasedVector[loc_i]->z() > dMaxVertexZ) || (locTrackTimeBasedVector[loc_i]->z() < dMinVertexZ))
			continue; //vertex-z cut: ignore total garbage
		JObject::oid_t locCandidateID = locTrackTimeBasedVector[loc_i]->candidateid;
		if(locBestTrackTimeBasedMap.find(locCandidateID) == locBestTrackTimeBasedMap.end())
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
		else if(locTrackTimeBasedVector[loc_i]->FOM > locBestTrackTimeBasedMap[locCandidateID]->FOM)
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
	}

	//separate the tracks based on high/low (potentially garbage) tracking FOM
	map<JObject::oid_t, const DTrackTimeBased*>::iterator locIterator;
	vector<const DTrackTimeBased*> locTrackTimeBasedVector_OnePerTrack, locTrackTimeBasedVector_OnePerTrack_GoodFOM;
	for(locIterator = locBestTrackTimeBasedMap.begin(); locIterator != locBestTrackTimeBasedMap.end(); ++locIterator)
	{
		const DTrackTimeBased* locTrackTimeBased = locIterator->second;
		locTrackTimeBasedVector_OnePerTrack.push_back(locTrackTimeBased);
		if(locTrackTimeBased->FOM > dMinTrackingFOM)
			locTrackTimeBasedVector_OnePerTrack_GoodFOM.push_back(locTrackTimeBased);
	}

	//project the track time to the beamline, then propagate that time to the center of the target (so that can compare with the RF time)
	//then figure out which RF bunch matches the most # tracks
		//need to try different sources of times: start with best quality
	int locBestRFBunchShift = 0; //default is no shift
	vector<pair<double, double> > locTimeFOMPairs;
	if(Find_TimeFOMPairs_ST(locSCHits, locTrackTimeBasedVector_OnePerTrack_GoodFOM, locTimeFOMPairs)) //good tracks, use ST info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);
	else if(Find_TimeFOMPairs_ST(locSCHits, locTrackTimeBasedVector_OnePerTrack, locTimeFOMPairs)) //potentially bad tracks, use ST info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);
	else if(Find_TimeFOMPairs_T0(locTrackTimeBasedVector_OnePerTrack_GoodFOM, locTimeFOMPairs)) //good tracks, use tracking time info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);
	else if(Find_TimeFOMPairs_T0(locTrackTimeBasedVector_OnePerTrack, locTimeFOMPairs)) //potentially bad tracks, use tracking time info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);

	DEventRFBunch *locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = locRFHitTime + dRFBunchFrequency*double(locBestRFBunchShift);
	locEventRFBunch->dTimeVariance = locTimeVariance;
	_data.push_back(locEventRFBunch);

	return NOERROR;
}

bool DEventRFBunch_factory::Find_TimeFOMPairs_ST(vector<const DSCHit*>& locSCHitVector, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, double> >& locTimeFOMPairs)
{
	locTimeFOMPairs.clear();
	for(unsigned int loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];

		// Match to the start counter using the result of the time-based fit
		double locProjectedSTTime = locTrackTimeBased->t0(); // to reject hits that are not in time with the track
		unsigned int locSCIndex;
		double locPathLength, locFlightTime;
		if(dPIDAlgorithm->MatchToSC(locTrackTimeBased->rt, DTrackFitter::kTimeBased, locSCHitVector, locProjectedSTTime, locSCIndex, locPathLength, locFlightTime) != NOERROR)
			continue;

		double locPropagatedTime = locProjectedSTTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
		locTimeFOMPairs.push_back(pair<double, double>(locPropagatedTime, locTrackTimeBased->FOM));
	}
	return (!locTimeFOMPairs.empty());
}

bool DEventRFBunch_factory::Find_TimeFOMPairs_T0(const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, double> >& locTimeFOMPairs)
{
	locTimeFOMPairs.clear();
	for(unsigned int loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];
		double locPropagatedTime = locTrackTimeBased->t0() + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
		locTimeFOMPairs.push_back(pair<double, double>(locPropagatedTime, locTrackTimeBased->FOM));
	}
	return (!locTimeFOMPairs.empty());
}

int DEventRFBunch_factory::Find_BestRFBunchShift(double locRFHitTime, const vector<pair<double, double> >& locTimeFOMPairs)
{
	//match DTrackTimeBased objects to ST hits, then project the ST hit time to the beamline (using the assumed PID) (ignore tracks with FOM < dMinTrackingFOM)
	//then find the #beam buckets the RF time needs to shift to match it
	map<int, pair<unsigned int, double> > locNumRFBucketsShiftedMap; //key is # RF buckets the RF time is shifted to match the track, value is the # tracks at that shift & the avg FOM of those tracks
	int locBestRFBunchShift = 0;
	for(unsigned int loc_i = 0; loc_i < locTimeFOMPairs.size(); ++loc_i)
	{
		double locPropagatedTrackTime = locTimeFOMPairs[loc_i].first;
		double locTrackFOM = locTimeFOMPairs[loc_i].second;
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
			locNumRFBucketsShiftedMap[locNumRFBucketsShifted] = pair<unsigned int, double>(1, locTrackFOM);
		else
		{
			unsigned int locOldNumTracks = locNumRFBucketsShiftedMap[locNumRFBucketsShifted].first;
			double locOldAverageFOM = locNumRFBucketsShiftedMap[locNumRFBucketsShifted].second;
			++(locNumRFBucketsShiftedMap[locNumRFBucketsShifted].first);
			locNumRFBucketsShiftedMap[locNumRFBucketsShifted].second = (locOldAverageFOM*(double(locOldNumTracks)) + locTrackFOM)/(double(locOldNumTracks + 1));
		}

		if(locNumRFBucketsShifted == locBestRFBunchShift)
			continue;

		unsigned int locBestNumTracks = locNumRFBucketsShiftedMap[locBestRFBunchShift].first;
		unsigned int locNumTracks = locNumRFBucketsShiftedMap[locNumRFBucketsShifted].first;
		if(locNumTracks > locBestNumTracks)
			locBestRFBunchShift = locNumRFBucketsShifted;
		else if(locNumTracks == locBestNumTracks)
		{
			if(locNumRFBucketsShiftedMap[locNumRFBucketsShifted].second > locNumRFBucketsShiftedMap[locBestRFBunchShift].second)
				locBestRFBunchShift = locNumRFBucketsShifted; //avg FOM higher in new bunch
		}
	}
	return locBestRFBunchShift;
}

//------------------
// erun
//------------------
jerror_t DEventRFBunch_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventRFBunch_factory::fini(void)
{
	return NOERROR;
}

