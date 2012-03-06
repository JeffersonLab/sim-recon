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


bool Compare_NeutralParticleHypotheses_FOM(const DNeutralParticleHypothesis *locTrack1, const DNeutralParticleHypothesis *locTrack2){
	return (locTrack1->dFOM > locTrack2->dFOM);
};

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
	unsigned int loc_i, loc_j;
	const DNeutralParticleHypothesis* locNeutralParticleHypothesis;
	DNeutralParticle *locNeutralParticle;
	bool locIDMatchFlag;
	vector<const DNeutralShower*> locNeutralShowers;
	vector<const DNeutralShower*> locNeutralShowers_Stored;

	vector<const DNeutralParticleHypothesis*> locNeutralParticleHypotheses;
	locEventLoop->Get(locNeutralParticleHypotheses);

	for(loc_i = 0; loc_i < locNeutralParticleHypotheses.size(); loc_i++){
		locNeutralParticleHypothesis = locNeutralParticleHypotheses[loc_i];
		locNeutralParticleHypothesis->GetT(locNeutralShowers);
		locIDMatchFlag = false;
		for (loc_j = 0; loc_j < _data.size(); loc_j++){
			_data[loc_j]->dNeutralParticleHypotheses[0]->GetT(locNeutralShowers_Stored);
			if(locNeutralShowers[0]->id == locNeutralShowers_Stored[0]->id){
				_data[loc_j]->dNeutralParticleHypotheses.push_back(locNeutralParticleHypothesis);
				locIDMatchFlag = true;
				break;
			}
		}
		if(locIDMatchFlag == true)
			continue;
		locNeutralParticle = new DNeutralParticle();
		locNeutralParticle->dNeutralParticleHypotheses.push_back(locNeutralParticleHypothesis);
		locNeutralParticle->AddAssociatedObject(locNeutralShowers[0]);
		_data.push_back(locNeutralParticle);
	}

	for(loc_i = 0; loc_i < _data.size(); loc_i++)
      sort(_data[loc_i]->dNeutralParticleHypotheses.begin(), _data[loc_i]->dNeutralParticleHypotheses.end(), Compare_NeutralParticleHypotheses_FOM);

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


