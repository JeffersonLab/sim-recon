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
	DApplication* locApplication = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = locApplication->GetDGeometry(runnumber);

	dMinTrackingFOM = 0.001;
	dRFBunchFrequency = 2.004;

	double locTargetCenterZ;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);

	double locTargetLength;
	locGeometry->GetTargetLength(locTargetLength);
	dMinVertexZ = locTargetCenterZ - 0.5*locTargetLength - 5.0;
	dMaxVertexZ = locTargetCenterZ + 0.5*locTargetLength + 5.0;

	locEventLoop->GetSingle(dParticleID);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	// https://halldweb1.jlab.org/wiki/index.php/How_HDGeant_defines_time-zero_for_physics_events

/*
	//COMMENTED UNTIL DRFTime CREATED
	vector<const DRFTime*> locRFTimes;
	locEventLoop->Get(locRFTimes);
	if(locRFTimes.empty())
		return RESOURCE_NOT_AVAILABLE;
	double locRFHitTime = locRFTimes[0]->dTime; //get from raw hit when ready
	double locTimeVariance = locRFTimes[0]->dTimeVariance;
*/

	//cheat & disable RF bunch selection until timing/etc. issues fixed
	vector<const DBeamPhoton*> locGenBeamPhotons;
	locEventLoop->Get(locGenBeamPhotons, "MCGEN");

	if(!locGenBeamPhotons.empty())
	{
		DEventRFBunch *locEventRFBunch = new DEventRFBunch;
		locEventRFBunch->dTime = locGenBeamPhotons[0]->time();
		locEventRFBunch->dTimeVariance = 0.0;
		locEventRFBunch->dMatchedToTracksFlag = true;
		_data.push_back(locEventRFBunch);
		return NOERROR;
	}

	DEventRFBunch *locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = 0.0;
	locEventRFBunch->dTimeVariance = 0.0;
	locEventRFBunch->dMatchedToTracksFlag = true;
	_data.push_back(locEventRFBunch);
	return NOERROR;

	double locRFHitTime = 0.0;
	double locTimeVariance = 0.0;

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	vector<const DSCHit*> locSCHits;
	locEventLoop->Get(locSCHits);

	vector<const DTOFPoint*> locTOFPoints;
	locEventLoop->Get(locTOFPoints);

	vector<const DBCALShower*> locBCALShowers;
	locEventLoop->Get(locBCALShowers);

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

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
	int locBestRFBunchShift = 0;
	bool locMatchedToTracksFlag = true; //set to false if not
	vector<pair<double, double> > locTimeFOMPairs;
	if(Find_TimeFOMPairs_Hits(locDetectorMatches, locTrackTimeBasedVector_OnePerTrack_GoodFOM, locTimeFOMPairs)) //good tracks, use TOF/BCAL/ST info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);
	else if(Find_TimeFOMPairs_Hits(locDetectorMatches, locTrackTimeBasedVector_OnePerTrack, locTimeFOMPairs)) //potentially bad tracks, use TOF/BCAL/ST info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);
	else if(Find_TimeFOMPairs_T0(locTrackTimeBasedVector_OnePerTrack_GoodFOM, locTimeFOMPairs)) //good tracks, use tracking time info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);
	else if(Find_TimeFOMPairs_T0(locTrackTimeBasedVector_OnePerTrack, locTimeFOMPairs)) //potentially bad tracks, use tracking time info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimeFOMPairs);
	else
		locMatchedToTracksFlag = false; //no confidence in selecting the RF bunch for the event

	locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = (locMatchedToTracksFlag) ? locRFHitTime + dRFBunchFrequency*double(locBestRFBunchShift) : locRFHitTime;
	locEventRFBunch->dTimeVariance = locTimeVariance;
	locEventRFBunch->dMatchedToTracksFlag = locMatchedToTracksFlag;
	_data.push_back(locEventRFBunch);

	return NOERROR;
}

bool DEventRFBunch_factory::Find_TimeFOMPairs_Hits(const DDetectorMatches* locDetectorMatches, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<pair<double, double> >& locTimeFOMPairs)
{
	locTimeFOMPairs.clear();

	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];

		//Use TOF/BCAL time to match to RF bunches:
			//max can be off = 1.002 ns (if off by more, will pick wrong beam bunch (RF frequency = 2.004 ns))

		// TOF
		DTOFHitMatchParams locTOFHitMatchParams;
		if(dParticleID->Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
		{
			//TOF resolution = 80ps, 3sigma = 240ps
			//max PID delta-t = 762ps (assuming no buffer)
				//for pion-proton: delta-t is ~750ps at ~170 MeV/c or lower: cannot have proton this slow anyway
				//for pion-kaon delta-t is ~750ps at ~80 MeV/c or lower: won't reconstruct these anyway, and not likely to be foreward-going anyway
			const DTOFPoint* locTOFPoint = locTOFHitMatchParams.dTOFPoint;
			double locPropagatedTime = locTOFPoint->t - locTOFHitMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
			locTimeFOMPairs.push_back(pair<double, double>(locPropagatedTime, locTrackTimeBased->FOM));
			continue;
		}

		// Else match to BCAL if fast enough (low time resolution for slow particles)
		double locP = locTrackTimeBased->momentum().Mag();
		//at 225 MeV/c, BCAL t-resolution is ~333ps (3sigma = 999ps), BCAL delta-t error is ~40ps: ~1040ps: bad
		//at 250 MeV/c, BCAL t-resolution is ~310ps (3sigma = 920ps), BCAL delta-t error is ~40ps: ~960 ps < 1 ns: OK!!
		if(locP < 0.25)
			continue; //too slow for the BCAL to distinguish RF bunch
		DShowerMatchParams locBCALShowerMatchParams;
		if(dParticleID->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
		{
			const DBCALShower* locBCALShower = dynamic_cast<const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
			double locPropagatedTime = locBCALShower->t - locBCALShowerMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
			locTimeFOMPairs.push_back(pair<double, double>(locPropagatedTime, locTrackTimeBased->FOM));
			continue;
		}
	}

	if(!locTimeFOMPairs.empty())
		return true;

	//no good matches from TOF/BCAL (or too slow for BCAL), try using ST
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];
		DSCHitMatchParams locSCHitMatchParams;
		if(dParticleID->Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams))
		{
			double locPropagatedTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
			locTimeFOMPairs.push_back(pair<double, double>(locPropagatedTime, locTrackTimeBased->FOM));
		}
		// Else no match, don't use track
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

