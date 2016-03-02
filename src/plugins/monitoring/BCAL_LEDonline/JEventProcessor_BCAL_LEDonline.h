// $Id$
//
//    File: JEventProcessor_BCAL_LEDonline.h
//

#ifndef _JEventProcessor_BCAL_LEDonline_
#define _JEventProcessor_BCAL_LEDonline_

#include <JANA/JEventProcessor.h>


class JEventProcessor_BCAL_LEDonline:public jana::JEventProcessor{
 public:
  JEventProcessor_BCAL_LEDonline();
  ~JEventProcessor_BCAL_LEDonline();
  const char* className(void){return "JEventProcessor_BCAL_LEDonline";}

  int NOtrig, GTPtrig, FPtrig, FPGTPtrig, trigUS, trigDS, trigCosmic;

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_BCAL_LEDonline_

