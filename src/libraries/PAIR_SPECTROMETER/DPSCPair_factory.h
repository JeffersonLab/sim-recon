// $Id$
//
//    File: DPSCPair_factory.h
// Created: Tue Mar 24 21:35:49 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//

#ifndef _DPSCPair_factory_
#define _DPSCPair_factory_

#include <JANA/JFactory.h>
#include "DPSCPair.h"

class DPSCPair_factory:public jana::JFactory<DPSCPair>{
 public:
  DPSCPair_factory(){};
  ~DPSCPair_factory(){};

  double DELTA_T_PAIR_MAX;

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPSCPair_factory_

