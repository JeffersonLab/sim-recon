// $Id$
//
//    File: DNeutralTrack_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DNeutralTrack_factory.h"
using namespace jana;


bool Compare_NeutralTrackHypotheses_FOM(const DNeutralTrackHypothesis *locTrack1, const DNeutralTrackHypothesis *locTrack2){
	return (locTrack1->dFOM > locTrack2->dFOM);
};

//------------------
// init
//------------------
jerror_t DNeutralTrack_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralTrack_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralTrack_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	unsigned int loc_i, loc_j;
	const DNeutralTrackHypothesis* locNeutralTrackHypothesis;
	DNeutralTrack *locNeutralTrack;
	bool locIDMatchFlag;
	vector<const DNeutralShowerCandidate*> locNeutralShowerCandidates;
	vector<const DNeutralShowerCandidate*> locNeutralShowerCandidates_Stored;

	vector<const DNeutralTrackHypothesis*> locNeutralTrackHypotheses;
	locEventLoop->Get(locNeutralTrackHypotheses);

	for(loc_i = 0; loc_i < locNeutralTrackHypotheses.size(); loc_i++){
		locNeutralTrackHypothesis = locNeutralTrackHypotheses[loc_i];
		locNeutralTrackHypothesis->GetT(locNeutralShowerCandidates);
		locIDMatchFlag = false;
		for (loc_j = 0; loc_j < _data.size(); loc_j++){
			_data[loc_j]->dNeutralTrackHypotheses[0]->GetT(locNeutralShowerCandidates_Stored);
			if(locNeutralShowerCandidates[0]->id == locNeutralShowerCandidates_Stored[0]->id){
				_data[loc_j]->dNeutralTrackHypotheses.push_back(locNeutralTrackHypothesis);
				locIDMatchFlag = true;
				break;
			}
		}
		if(locIDMatchFlag == true)
			continue;
		locNeutralTrack = new DNeutralTrack();
		locNeutralTrack->dNeutralTrackHypotheses.push_back(locNeutralTrackHypothesis);
		locNeutralTrack->AddAssociatedObject(locNeutralShowerCandidates[0]);
		_data.push_back(locNeutralTrack);
	}

	for(loc_i = 0; loc_i < _data.size(); loc_i++)
      sort(_data[loc_i]->dNeutralTrackHypotheses.begin(), _data[loc_i]->dNeutralTrackHypotheses.end(), Compare_NeutralTrackHypotheses_FOM);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralTrack_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralTrack_factory::fini(void)
{
	return NOERROR;
}


