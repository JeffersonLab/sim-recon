#ifndef _DStartTime_factory_
#define _DStartTime_factory_

#include <JANA/JObject.h>
#include <JANA/JFactory.h>

#include "DStartTime.h"

class DStartTime_factory: public jana::JFactory<DStartTime> {
 public:
  DStartTime_factory() {}
  ~DStartTime_factory() {}

 private:
  jerror_t init();						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *loop, int runnumber);	        ///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *loop, int eventnumber);	///< Called every event.
  jerror_t erun();						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini();						///< Called after last event of last event source has been processed.

};

#endif //_DStartTime_factory_
