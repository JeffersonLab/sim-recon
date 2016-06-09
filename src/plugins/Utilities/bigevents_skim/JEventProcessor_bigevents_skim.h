// $Id$
//
//    File: JEventProcessor_bigevents_skim.h
// Created: Thu May 12 08:01:59 EDT 2016
// Creator: zihlmann (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_bigevents_skim_
#define _JEventProcessor_bigevents_skim_

#include <JANA/JEventProcessor.h>
#include "evio_writer/DEventWriterEVIO.h"
#include <TRIGGER/DL1Trigger.h>

#include <CDC/DCDCDigiHit.h>
using namespace std;
using namespace jana;


class JEventProcessor_bigevents_skim:public jana::JEventProcessor{
 public:
  JEventProcessor_bigevents_skim();
  ~JEventProcessor_bigevents_skim();
  const char* className(void){return "JEventProcessor_bigevents_skim";}
  

 private:
  jerror_t init(void);///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);
  jerror_t erun(void);///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_bigevents_skim_

