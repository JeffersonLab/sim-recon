// $Id$
//
//    File: JEventProcessor_PSPair_online.h
// Created: Fri Mar 20 16:32:04 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_PSPair_online_
#define _JEventProcessor_PSPair_online_

#include <JANA/JEventProcessor.h>

class JEventProcessor_PSPair_online:public jana::JEventProcessor{
 public:
  JEventProcessor_PSPair_online();
  ~JEventProcessor_PSPair_online();
  const char* className(void){return "JEventProcessor_PSPair_online";}

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_PSPair_online_

