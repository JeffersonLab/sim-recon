// $Id$
//
//    File: DChargedTrack_factory_PreSelectTimeCalib.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DChargedTrack_factory_PreSelectTimeCalib.h"

//------------------
// init
//------------------
jerror_t DChargedTrack_factory_PreSelectTimeCalib::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrack_factory_PreSelectTimeCalib::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrack_factory_PreSelectTimeCalib::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DTrackWireBased*> locTrackWireBasedVector;
	locEventLoop->Get(locTrackWireBasedVector);

	//arrange by track
	map<JObject::oid_t, vector<const DTrackWireBased*> > locWireBasedMap;
	for(size_t loc_i = 0; loc_i < locTrackWireBasedVector.size(); ++loc_i)
	{
		const DTrackWireBased* locTrackWireBased = locTrackWireBasedVector[loc_i];
		locWireBasedMap[locTrackWireBased->candidateid].push_back(locTrackWireBased);
	}

	//create charged tracks
	map<JObject::oid_t, vector<const DTrackWireBased*> >::iterator locIterator = locWireBasedMap.begin();
	for(; locIterator != locWireBasedMap.end(); ++locIterator)
	{
		DChargedTrack* locChargedTrack = new DChargedTrack();
		vector<const DTrackWireBased*>& locWireBasedTracks = locIterator->second;
		for(size_t loc_i = 0; loc_i < locWireBasedTracks.size(); ++loc_i)
		{
			DChargedTrackHypothesis* locChargedTrackHypothesis = new DChargedTrackHypothesis();
			DKinematicData *locKinematicData = locChargedTrackHypothesis;
			*locKinematicData = *(static_cast<const DKinematicData*>(locWireBasedTracks[loc_i]));

			locChargedTrack->dChargedTrackHypotheses.push_back(locChargedTrackHypothesis);
		}

		_data.push_back(locChargedTrack);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DChargedTrack_factory_PreSelectTimeCalib::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTrack_factory_PreSelectTimeCalib::fini(void)
{
	return NOERROR;
}
