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
	SetFactoryFlag(NOT_OBJECT_OWNER);
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DMCThrown_factory_FinalState::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	locEventLoop->GetSingle(dAnalysisUtilities);
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DMCThrown_factory_FinalState::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	_data.clear();

	deque<pair<const DMCThrown*, deque<const DMCThrown*> > > locThrownSteps;
	dAnalysisUtilities->Get_ThrownParticleSteps(locEventLoop, locThrownSteps);

	if(locThrownSteps.empty())
		return NOERROR;

	for(size_t loc_i = 0; loc_i < locThrownSteps.size(); ++loc_i)
	{
		deque<const DMCThrown*>& locParticles = locThrownSteps[loc_i].second;
		for(size_t loc_j = 0; loc_j < locParticles.size(); ++loc_j)
		{
			if(Is_FinalStateParticle(locParticles[loc_j]->PID()))
				_data.push_back(const_cast<DMCThrown*>(locParticles[loc_j]));
		}
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

