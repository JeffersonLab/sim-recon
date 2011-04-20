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
  vector<const DParticleSet *>particle_sets;
  loop->Get(particle_sets);

  // DPhysicsEvent for each particle set
  for(unsigned int i=0; i<particle_sets.size(); i++){
    DPhysicsEvent *pe = new DPhysicsEvent;
    
    pe->particle_sets.push_back(particle_sets[i]);
    
		
    // "Publish" the DPhysicsEvent object 
    _data.push_back(pe);
  }
	
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
