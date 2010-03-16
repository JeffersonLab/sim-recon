// $Id$
//
//    File: DPhoton_factory_THROWN.cc
// Created: Tue Mar  9 22:34:29 EST 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <TRACKING/DMCThrown.h>

#include "DPhoton_factory_THROWN.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPhoton_factory_THROWN::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPhoton_factory_THROWN::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPhoton_factory_THROWN::evnt(JEventLoop *loop, int eventnumber)
{

	vector<const DMCThrown*> throwns;
	loop->Get(throwns);
	
	for(unsigned int i=0; i<throwns.size(); i++){
		const DMCThrown *thrown = throwns[i];
		if(thrown->charge()!=0.0)continue;
		
		DPhoton *photon = new DPhoton;

		// Copy DKinematicData part
		DKinematicData *kd_photon = photon;
		const DKinematicData *kd_thrown = thrown;
		*kd_photon = *kd_thrown;
		
		photon->setTag(DPhoton::kFcal); // just set everything to FCAL for now
		
		_data.push_back(photon);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPhoton_factory_THROWN::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPhoton_factory_THROWN::fini(void)
{
	return NOERROR;
}

