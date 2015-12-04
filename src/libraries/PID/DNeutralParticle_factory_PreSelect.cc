// $Id$
//
//    File: DNeutralParticle_factory_PreSelect.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include "DNeutralParticle_factory_PreSelect.h"

//------------------
// init
//------------------
jerror_t DNeutralParticle_factory_PreSelect::init(void)
{
	//Setting this flag makes it so that JANA does not delete the objects in _data.  This factory will manage this memory. 
		//This is because some/all of these pointers are just copied from earlier objects, and should not be deleted.  
	SetFactoryFlag(NOT_OBJECT_OWNER);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticle_factory_PreSelect::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticle_factory_PreSelect::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	//Clear objects from last event
	_data.clear();

	vector<const DNeutralParticle*> locNeutralParticles;
	locEventLoop->Get(locNeutralParticles);

	vector<const DNeutralShower*> locNeutralShowers;
	locEventLoop->Get(locNeutralShowers, "PreSelect");

	set<const DNeutralShower*> locNeutralShowerSet;
	for(size_t loc_i = 0; loc_i < locNeutralShowers.size(); ++loc_i)
		locNeutralShowerSet.insert(locNeutralShowers[loc_i]);

	//if neutral shower was good, keep particle, else ignore it
	for(size_t loc_i = 0; loc_i < locNeutralParticles.size(); ++loc_i)
	{
		if(locNeutralShowerSet.find(locNeutralParticles[loc_i]->dNeutralShower) != locNeutralShowerSet.end())
			_data.push_back(const_cast<DNeutralParticle*>(locNeutralParticles[loc_i]));
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticle_factory_PreSelect::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticle_factory_PreSelect::fini(void)
{
	return NOERROR;
}


