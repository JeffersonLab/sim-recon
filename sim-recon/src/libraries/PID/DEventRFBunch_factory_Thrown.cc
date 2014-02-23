// $Id$
//
//    File: DEventRFBunch_factory_Thrown.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DEventRFBunch_factory_Thrown.h"
#include <deque>

using namespace jana;

//------------------
// init
//------------------
jerror_t DEventRFBunch_factory_Thrown::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DEventRFBunch_factory_Thrown::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DEventRFBunch_factory_Thrown::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);
	if(locMCThrowns.empty())
		return NOERROR;

	//t = 0 at target center: https://halldweb1.jlab.org/wiki/index.php/How_HDGeant_defines_time-zero_for_physics_events
	DEventRFBunch *locEventRFBunch = new DEventRFBunch;
	locEventRFBunch->dTime = 0.0;
	locEventRFBunch->dTimeVariance = 0.0;
	locEventRFBunch->dMatchedToTracksFlag = true;
	_data.push_back(locEventRFBunch);
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DEventRFBunch_factory_Thrown::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DEventRFBunch_factory_Thrown::fini(void)
{
	return NOERROR;
}

