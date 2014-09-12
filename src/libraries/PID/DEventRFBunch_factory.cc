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

	dMinTrackingFOM = 5.73303E-7; //5 sigma
	dRFBunchFrequency = 2.004;

	double locTargetCenterZ;
	locGeometry->GetTargetZ(locTargetCenterZ);
	dTargetCenter.SetXYZ(0.0, 0.0, locTargetCenterZ);

	locEventLoop->GetSingle(dParticleID);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DRFTime*> locRFTimes;
	locEventLoop->Get(locRFTimes);
	if(locRFTimes.empty())
		return RESOURCE_UNAVAILABLE;

	double locRFHitTime = locRFTimes[0]->dTime;
	double locTimeVariance = locRFTimes[0]->dTimeVariance;

	//preferentially:
	//use SC hits on tracks with good tracking FOM
		//if no SC hits, use other timing systems
	//if no good tracks (or none with matched hits), use SC hits on all tracks
		//if no SC hits, use other timing systems
	//if no tracks (or none with matched hits), use neutral showers assuming PIDs = photon
	//if no showers ... set time to the original RF time and set the good-flag to false

	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector);

	//select the best DTrackTimeBased for each track: use best tracking FOM
	map<JObject::oid_t, const DTrackTimeBased*> locBestTrackTimeBasedMap; //lowest tracking chisq/ndf for each candidate id
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		JObject::oid_t locCandidateID = locTrackTimeBasedVector[loc_i]->candidateid;
		if(locBestTrackTimeBasedMap.find(locCandidateID) == locBestTrackTimeBasedMap.end())
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
		else if(locTrackTimeBasedVector[loc_i]->FOM > locBestTrackTimeBasedMap[locCandidateID]->FOM)
			locBestTrackTimeBasedMap[locCandidateID] = locTrackTimeBasedVector[loc_i];
	}

	//separate the tracks based on high/low tracking FOM
	map<JObject::oid_t, const DTrackTimeBased*>::iterator locIterator;
	vector<const DTrackTimeBased*> locTrackTimeBasedVector_OnePerTrack, locTrackTimeBasedVector_OnePerTrack_GoodFOM;
	for(locIterator = locBestTrackTimeBasedMap.begin(); locIterator != locBestTrackTimeBasedMap.end(); ++locIterator)
	{
		const DTrackTimeBased* locTrackTimeBased = locIterator->second;
		locTrackTimeBasedVector_OnePerTrack.push_back(locTrackTimeBased);
		if(locTrackTimeBased->FOM >= dMinTrackingFOM)
			locTrackTimeBasedVector_OnePerTrack_GoodFOM.push_back(locTrackTimeBased);
	}

	const DDetectorMatches* locDetectorMatches = NULL;
	locEventLoop->GetSingle(locDetectorMatches);

	//then figure out which RF bunch matches the most # tracks
		//need to try different sources of times: start with best quality
	int locBestRFBunchShift = 0;
	bool locMatchedToTracksFlag = true; //set to false if not
	vector<double> locTimes;
	if(Find_TrackTimes(locDetectorMatches, locTrackTimeBasedVector_OnePerTrack_GoodFOM, locTimes)) //good tracks, use TOF/BCAL/ST info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimes);
	else if(Find_TrackTimes(locDetectorMatches, locTrackTimeBasedVector_OnePerTrack, locTimes)) //potentially bad tracks, use TOF/BCAL/ST info
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimes);
	else if(Find_NeutralTimes(locEventLoop, locTimes))//use neutral showers
		locBestRFBunchShift = Find_BestRFBunchShift(locRFHitTime, locTimes);
	else
		locMatchedToTracksFlag = false; //no confidence in selecting the RF bunch for the event
	
	DEventRFBunch *locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = locRFHitTime + dRFBunchFrequency*double(locBestRFBunchShift);
	locEventRFBunch->dTimeVariance = locTimeVariance;
	locEventRFBunch->dMatchedToTracksFlag = locMatchedToTracksFlag;
	_data.push_back(locEventRFBunch);

	return NOERROR;
}

bool DEventRFBunch_factory::Find_NeutralTimes(JEventLoop* locEventLoop, vector<double>& locTimes)
{
	locTimes.clear();

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses);

	for(size_t loc_i = 0; loc_i < locNeutralParticleHypotheses.size(); ++loc_i)
	{
		if(locNeutralParticleHypotheses[loc_i]->PID() == Gamma)
			locTimes.push_back(locNeutralParticleHypotheses[loc_i]->time());
	}

	return (!locTimes.empty());
}

bool DEventRFBunch_factory::Find_TrackTimes(const DDetectorMatches* locDetectorMatches, const vector<const DTrackTimeBased*>& locTrackTimeBasedVector, vector<double>& locTimes)
{
	locTimes.clear();

	//SC
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];
		DSCHitMatchParams locSCHitMatchParams;
		if(!dParticleID->Get_BestSCMatchParams(locTrackTimeBased, locDetectorMatches, locSCHitMatchParams))
			continue;

		double locPropagatedTime = locSCHitMatchParams.dHitTime - locSCHitMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
		locTimes.push_back(locPropagatedTime);
	}

	if(!locTimes.empty())
		return true;

	//None of the input tracks have SC hits: try other timing systems
		//Order of preference: BCAL, TOF, FCAL
	for(size_t loc_i = 0; loc_i < locTrackTimeBasedVector.size(); ++loc_i)
	{
		const DTrackTimeBased* locTrackTimeBased = locTrackTimeBasedVector[loc_i];

		DShowerMatchParams locBCALShowerMatchParams;
		if(dParticleID->Get_BestBCALMatchParams(locTrackTimeBased, locDetectorMatches, locBCALShowerMatchParams))
		{
			const DBCALShower* locBCALShower = dynamic_cast<const DBCALShower*>(locBCALShowerMatchParams.dShowerObject);
			double locPropagatedTime = locBCALShower->t - locBCALShowerMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
			locTimes.push_back(locPropagatedTime);
			continue;
		}

		DTOFHitMatchParams locTOFHitMatchParams;
		if(dParticleID->Get_BestTOFMatchParams(locTrackTimeBased, locDetectorMatches, locTOFHitMatchParams))
		{
			const DTOFPoint* locTOFPoint = locTOFHitMatchParams.dTOFPoint;
			double locPropagatedTime = locTOFPoint->t - locTOFHitMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
			locTimes.push_back(locPropagatedTime);
			continue;
		}

		DShowerMatchParams locFCALShowerMatchParams;
		if(dParticleID->Get_BestFCALMatchParams(locTrackTimeBased, locDetectorMatches, locFCALShowerMatchParams))
		{
			const DFCALShower* locFCALShower = dynamic_cast<const DFCALShower*>(locFCALShowerMatchParams.dShowerObject);
			double locPropagatedTime = locFCALShower->getTime() - locFCALShowerMatchParams.dFlightTime + (dTargetCenter.Z() - locTrackTimeBased->z())/SPEED_OF_LIGHT;
			locTimes.push_back(locPropagatedTime);
			continue;
		}
	}

	return (!locTimes.empty());
}

int DEventRFBunch_factory::Find_BestRFBunchShift(double locRFHitTime, const vector<double>& locTimes)
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
		{
			unsigned int locOldNumTracks = locNumRFBucketsShiftedMap[locNumRFBucketsShifted];
			++(locNumRFBucketsShiftedMap[locNumRFBucketsShifted]);
		}

		if(locNumRFBucketsShifted == locBestRFBunchShift)
			continue;

		unsigned int locBestNumTracks = locNumRFBucketsShiftedMap[locBestRFBunchShift];
		unsigned int locNumTracks = locNumRFBucketsShiftedMap[locNumRFBucketsShifted];
		if(locNumTracks > locBestNumTracks)
			locBestRFBunchShift = locNumRFBucketsShifted;
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

