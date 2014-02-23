// $Id$
//
//    File: DParticleID_factory.h
// Created: Mon Feb 28 13:47:49 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_factory_
#define _DParticleID_factory_

#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include "DParticleID.h"

class DParticleID_factory:public jana::JFactory<DParticleID>{
 public:
  DParticleID_factory(){};
  ~DParticleID_factory(){};
  
  
 private:
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber){
    // This is a trivial class that simply implements a default
    // factory. It is here so that the default can be changed 
    // easily by simply changing the tag here or on the command 
    // line.
    vector<const DParticleID*> pid_algorithms;
    eventLoop->Get(pid_algorithms,"PID1");
    for(unsigned int i=0; i< pid_algorithms.size(); i++){
      _data.push_back(const_cast<DParticleID*>(pid_algorithms[i]));
    }
    SetFactoryFlag(NOT_OBJECT_OWNER);
    
    return NOERROR;
      
  };	///< Called every event.
};

#endif // _DParticleID_factory_

