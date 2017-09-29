// $Id$
//
//    File: DBeamPhoton_factory_TRUTH.cc
// Created: Thu Dec  3 17:27:55 EST 2009
// Creator: staylor (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DBeamPhoton_factory_TRUTH.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DBeamPhoton_factory_TRUTH::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBeamPhoton_factory_TRUTH::brun(jana::JEventLoop *locEventLoop, int32_t runnumber)
{
	DApplication* dapp = dynamic_cast<DApplication*>(locEventLoop->GetJApplication());
	DGeometry* locGeometry = dapp->GetDGeometry(locEventLoop->GetJEvent().GetRunNumber());
	dTargetCenterZ = 0.0;
	locGeometry->GetTargetZ(dTargetCenterZ);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBeamPhoton_factory_TRUTH::evnt(jana::JEventLoop *locEventLoop, uint64_t eventnumber)
{
	DVector3 pos(0.0, 0.0, dTargetCenterZ);

	vector<const DTAGMHit*> tagm_hits;
	locEventLoop->Get(tagm_hits, "TRUTH");
	for (unsigned int ih=0; ih < tagm_hits.size(); ++ih)
	{
		if (tagm_hits[ih]->row > 0) continue;
		DVector3 mom(0.0, 0.0, tagm_hits[ih]->E);
		DBeamPhoton *gamma = new DBeamPhoton;
		gamma->setPID(Gamma);
		gamma->setMomentum(mom);
		gamma->setPosition(pos);
		gamma->setTime(tagm_hits[ih]->t);
	    gamma->dSystem = SYS_TAGM;
		gamma->dCounter = tagm_hits[ih]->column;
		gamma->AddAssociatedObject(tagm_hits[ih]);
		_data.push_back(gamma);
	}

	vector<const DTAGHHit*> tagh_hits;
	locEventLoop->Get(tagh_hits, "TRUTH");
	for (unsigned int ih=0; ih < tagh_hits.size(); ++ih)
	{
		DVector3 mom(0.0, 0.0, tagh_hits[ih]->E);
		DBeamPhoton *gamma = new DBeamPhoton;
		gamma->setPID(Gamma);
		gamma->setMomentum(mom);
		gamma->setPosition(pos);
		gamma->setTime(tagh_hits[ih]->t);
	    gamma->dSystem = SYS_TAGH;
		gamma->dCounter = tagh_hits[ih]->counter_id;
		gamma->AddAssociatedObject(tagh_hits[ih]);
		_data.push_back(gamma);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBeamPhoton_factory_TRUTH::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBeamPhoton_factory_TRUTH::fini(void)
{
	return NOERROR;
}

