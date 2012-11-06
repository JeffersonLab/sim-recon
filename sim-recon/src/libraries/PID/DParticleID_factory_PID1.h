// $Id$
//
//    File: DParticleID_factory_PID1.h
// Created: Mon Feb 28 14:12:16 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_factory_PID1_
#define _DParticleID_factory_PID1_

#include <JANA/JFactory.h>
#include "DParticleID_PID1.h"

class DParticleID_factory_PID1:public jana::JFactory<DParticleID>{
 public:
  DParticleID_factory_PID1(){};
  ~DParticleID_factory_PID1(){};
  const char* Tag(void){return "PID1";}
  
 private:
  jerror_t evnt(jana::JEventLoop *loop, int eventnumber){
    // Create single DParticleID object and mark the factory as
    // persistent so it doesn't get deleted every event.
    DParticleID *pid_algorithm = new DParticleID_PID1(loop);
    SetFactoryFlag(PERSISTANT);
    ClearFactoryFlag(WRITE_TO_OUTPUT);
    _data.push_back(pid_algorithm);

    return NOERROR;
  }

};

#endif // _DParticleID_factory_PID1_

