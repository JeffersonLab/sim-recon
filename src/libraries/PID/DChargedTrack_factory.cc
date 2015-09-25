// $Id$
//
//    File: DChargedTrack_factory.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DChargedTrack_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DChargedTrack_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrack_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrack_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses);

	map<JObject::oid_t, vector<const DChargedTrackHypothesis*> > locHypothesesByTrackID;
	for(size_t loc_i = 0; loc_i < locChargedTrackHypotheses.size(); loc_i++)
		locHypothesesByTrackID[locChargedTrackHypotheses[loc_i]->candidateid].push_back(locChargedTrackHypotheses[loc_i]);

	map<JObject::oid_t, vector<const DChargedTrackHypothesis*> >::iterator locIterator = locHypothesesByTrackID.begin();
	for(; locIterator != locHypothesesByTrackID.end(); ++locIterator)
	{
		DChargedTrack* locChargedTrack = new DChargedTrack();
		locChargedTrack->candidateid = locIterator->first;
		locChargedTrack->dChargedTrackHypotheses = locIterator->second;
		_data.push_back(locChargedTrack);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DChargedTrack_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTrack_factory::fini(void)
{
	return NOERROR;
}


