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

bool Compare_ChargedTrackHypotheses_FOM(const DChargedTrackHypothesis *locTrack1, const DChargedTrackHypothesis *locTrack2){
	return (locTrack1->dFOM > locTrack2->dFOM);
};

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
	unsigned int loc_i, loc_j;
	const DChargedTrackHypothesis* locChargedTrackHypothesis;
	DChargedTrack *locChargedTrack;
	bool locIDMatchFlag;

	vector<const DChargedTrackHypothesis*> locChargedTrackHypotheses;
	locEventLoop->Get(locChargedTrackHypotheses);

	for(loc_i = 0; loc_i < locChargedTrackHypotheses.size(); loc_i++){
		locChargedTrackHypothesis = locChargedTrackHypotheses[loc_i];
		locIDMatchFlag = false;
		for (loc_j = 0; loc_j < _data.size(); loc_j++){
			if(locChargedTrackHypothesis->candidateid == _data[loc_j]->dChargedTrackHypotheses[0]->candidateid){
				_data[loc_j]->dChargedTrackHypotheses.push_back(locChargedTrackHypothesis);
				locIDMatchFlag = true;
				break;
			}

		}
		if(locIDMatchFlag == true)
			continue;
		locChargedTrack = new DChargedTrack();
		locChargedTrack->dChargedTrackHypotheses.push_back(locChargedTrackHypothesis);
		_data.push_back(locChargedTrack);
	}

	for(loc_i = 0; loc_i < _data.size(); loc_i++)
      sort(_data[loc_i]->dChargedTrackHypotheses.begin(), _data[loc_i]->dChargedTrackHypotheses.end(), Compare_ChargedTrackHypotheses_FOM);

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


