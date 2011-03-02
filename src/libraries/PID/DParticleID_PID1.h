// $Id$
//
//    File: DParticleID_PID1.h
// Created: Mon Feb 28 15:25:35 EST 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleID_PID1_
#define _DParticleID_PID1_

#include <JANA/jerror.h>
#include "DParticleID.h"

class DParticleID_PID1:public DParticleID{
 public:
  DParticleID_PID1(JEventLoop *loop); // require JEventLoop in constructor;
  ~DParticleID_PID1();

  jerror_t GetdEdxChiSq(const DTrackTimeBased *track,double &dEdx,
			unsigned int &num,double &chi2);
 protected:
  
  
 private:
  int DEBUG_LEVEL;
  // Prohibit default constructor
  DParticleID_PID1();	

};

#endif // _DParticleID_PID1_

