// $Id$
//
//    File: JEventProcessor_merge_rawevents.h
//

#ifndef _JEventProcessor_merge_rawevents_
#define _JEventProcessor_merge_rawevents_

#include <JANA/JEventProcessor.h>
#include "evio_writer/DEventWriterEVIO.h"

#include <vector>

class JEventProcessor_merge_rawevents:public jana::JEventProcessor{
 public:

  JEventProcessor_merge_rawevents();
  ~JEventProcessor_merge_rawevents();
  const char* className(void){return "JEventProcessor_merge_rawevents";}

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  DEventWriterEVIO* dEventWriterEVIO;
};

#endif // _JEventProcessor_merge_rawevents_

