// $Id$
//
//    File: DChargedTrackHypothesis_factory_Reaction.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DChargedTrackHypothesis_factory_Reaction.h"
using namespace std;
using namespace jana;

//------------------
// init
//------------------
jerror_t DChargedTrackHypothesis_factory_Reaction::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTrackHypothesis_factory_Reaction::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTrackHypothesis_factory_Reaction::evnt(jana::JEventLoop* locEventLoop, int eventnumber)
{
	DChargedTrackHypothesis_factory* locChargedTrackHypothesisFactory = static_cast<DChargedTrackHypothesis_factory*>(locEventLoop->GetFactory("DChargedTrackHypothesis"));

 	vector<const DTrackTimeBased*> locTrackTimeBasedVector;
	locEventLoop->Get(locTrackTimeBasedVector, "Reaction");

 	vector<DChargedTrackHypothesis*> locChargedTrackHypotheses;
	jerror_t locJError = locChargedTrackHypothesisFactory->Get_ChargedTrackHypotheses(locEventLoop, locTrackTimeBasedVector, locChargedTrackHypotheses);
	if(locJError != NOERROR)
		return locJError;

	for(size_t loc_i = 0; loc_i < locChargedTrackHypotheses.size(); ++loc_i)
	{
	 	vector<const DChargedTrack*> locChargedTracks;
		locTrackTimeBasedVector[loc_i]->GetT(locChargedTracks);
		locChargedTrackHypotheses[loc_i]->AddAssociatedObject(locChargedTracks[0]);
		_data.push_back(locChargedTrackHypotheses[loc_i]);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DChargedTrackHypothesis_factory_Reaction::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTrackHypothesis_factory_Reaction::fini(void)
{
	return NOERROR;
}


