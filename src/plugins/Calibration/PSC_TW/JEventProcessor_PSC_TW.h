// $Id$
//
//    File: JEventProcessor_PSC_TW.h
// Created: Fri Aug 21 10:42:28 EDT 2015
// Creator: aebarnes (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_PSC_TW_
#define _JEventProcessor_PSC_TW_

#include <JANA/JEventProcessor.h>

#include <PAIR_SPECTROMETER/DPSCHit.h>
#include <PAIR_SPECTROMETER/DPSCPair.h>

class JEventProcessor_PSC_TW:public jana::JEventProcessor{
   public:
      JEventProcessor_PSC_TW();
      ~JEventProcessor_PSC_TW();
      const char* className(void){return "JEventProcessor_PSC_TW";}

   private:
      jerror_t init(void);						///< Called once at program start.
      jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
      jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
      jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
      jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_PSC_TW_

