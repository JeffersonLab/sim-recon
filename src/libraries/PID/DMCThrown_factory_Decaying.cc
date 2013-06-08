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
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCThrown_factory_Decaying::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCThrown_factory_Decaying::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		if(locMCThrowns[loc_i]->PID() == Unknown)
			continue;
		if(Is_FinalStateParticle(locMCThrowns[loc_i]->PID()) == 1)
			continue;
		DMCThrown* locMCThrown = new DMCThrown(*locMCThrowns[loc_i]);
		_data.push_back(locMCThrown);
	}

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

