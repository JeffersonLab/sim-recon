// $Id$
//
//    File: DNeutralParticle_factory.cc
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DNeutralParticle_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DNeutralParticle_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DNeutralParticle_factory::brun(jana::JEventLoop *locEventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DNeutralParticle_factory::evnt(jana::JEventLoop *locEventLoop, int eventnumber)
{
	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses);

	map<JObject::oid_t, vector<const DNeutralParticleHypothesis*> > locHypothesesByShowerID;
	for(size_t loc_i = 0; loc_i < locNeutralParticleHypotheses.size(); loc_i++)
		locHypothesesByShowerID[locNeutralParticleHypotheses[loc_i]->dNeutralShowerID].push_back(locNeutralParticleHypotheses[loc_i]);

	map<JObject::oid_t, vector<const DNeutralParticleHypothesis*> >::iterator locIterator = locHypothesesByShowerID.begin();
	for(; locIterator != locHypothesesByShowerID.end(); ++locIterator)
	{
		DNeutralParticle* locNeutralParticle = new DNeutralParticle();
		locNeutralParticle->dNeutralParticleHypotheses = locIterator->second;
		locIterator->second[0]->GetSingle(locNeutralParticle->dNeutralShower);
		locNeutralParticle->AddAssociatedObject(locNeutralParticle->dNeutralShower);
		_data.push_back(locNeutralParticle);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DNeutralParticle_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DNeutralParticle_factory::fini(void)
{
	return NOERROR;
}


