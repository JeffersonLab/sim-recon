// $Id$
//
//    File: JEventProcessor_TOF_online.h
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TOF_online_
#define _JEventProcessor_TOF_online_

#include <JANA/JEventProcessor.h>


class JEventProcessor_TOF_online:public jana::JEventProcessor{
 public:
  JEventProcessor_TOF_online();
  ~JEventProcessor_TOF_online();
  const char* className(void){return "JEventProcessor_TOF_online";}


 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TOF_online_

