// $Id$
//
//    File: DMCThrown_factory_FinalState.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DMCThrown_factory_FinalState.h"

using namespace jana;

//------------------
// init
//------------------
jerror_t DMCThrown_factory_FinalState::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCThrown_factory_FinalState::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCThrown_factory_FinalState::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DMCThrown*> locMCThrowns;
	locEventLoop->Get(locMCThrowns);

	for(size_t loc_i = 0; loc_i < locMCThrowns.size(); ++loc_i)
	{
		int locIsFinalStateParticleInt = Is_FinalStateParticle(locMCThrowns[loc_i]->PID());
		if(locIsFinalStateParticleInt != 1)
			continue;

		//find parent id, make sure not also a final state particle
		int locParentID = locMCThrowns[loc_i]->parentid;
		int locIsParentFinalStateInt = 0;
		bool locIsParentFoundFlag = false;
		if(locParentID != 0) //parent = 0 means initial production particle
		{
			for(size_t loc_j = 0; loc_j < locMCThrowns.size(); ++loc_j)
			{
				if(locMCThrowns[loc_j]->myid != locParentID)
					continue;
				locIsParentFinalStateInt = Is_FinalStateParticle(locMCThrowns[loc_j]->PID());
				locIsParentFoundFlag = true;
				break;
			}
		}
		if((locIsParentFinalStateInt == 1) || ((!locIsParentFoundFlag) && (locParentID != 0)))
			continue; //don't save: a decay product of a final state particle (e.g. mu+ from pi+ decay) //OR the parent is lost

		DMCThrown* locMCThrown = new DMCThrown(*locMCThrowns[loc_i]);
		_data.push_back(locMCThrown);
	}

	return NOERROR;
}


//------------------
// erun
//------------------
jerror_t DMCThrown_factory_FinalState::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DMCThrown_factory_FinalState::fini(void)
{
	return NOERROR;
}

