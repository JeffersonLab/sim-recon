// $Id$
//
//    File: DPhysicsEvent_factory.cc
// Created: Wed Aug  4 10:37:55 EDT 2010
// Creator: davidl (on Darwin eleanor.jlab.org 10.4.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DPhysicsEvent_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPhysicsEvent_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPhysicsEvent_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPhysicsEvent_factory::evnt(JEventLoop *loop, int eventnumber)
{
	//Only one DPhysicsEvent is generated.  
	DPhysicsEvent *locPhysicsEvent = new DPhysicsEvent;
	loop->Get(locPhysicsEvent->particle_sets); //DParticleSet

	// "Publish" the DPhysicsEvent object 
	_data.push_back(locPhysicsEvent);

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPhysicsEvent_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPhysicsEvent_factory::fini(void)
{
	return NOERROR;
}
