// $Id$
//
//    File: DMCThrown_factory_Decaying.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DMCThrown_factory_Decaying.h"

using namespace jana;

//------------------
// init
//------------------
jerror_t DMCThrown_factory_Decaying::init(void)
{
	SetFactoryFlag(NOT_OBJECT_OWNER);
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCThrown_factory_Decaying::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCThrown_factory_Decaying::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	_data.clear();

	deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locThrownSteps;
	dAnalysisUtilities->Get_ThrownParticleSteps(locEventLoop, locThrownSteps);

	if(locThrownSteps.empty())
		return NOERROR;

	for(size_t loc_i = 1; loc_i < locThrownSteps.size(); ++loc_i)
		_data.push_back(const_cast<DMCThrown*>(locThrownSteps[loc_i].first));

	return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DMCThrown_factory_Decaying::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DMCThrown_factory_Decaying::fini(void)
{
	return NOERROR;
}

