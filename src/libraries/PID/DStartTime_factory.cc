#include "DStartTime_factory.h"

jerror_t DStartTime_factory::init() {
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DStartTime_factory::brun(jana::JEventLoop *loop, int runnumber) {
  //  return RESOURCE_UNAVAILABLE;
  
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DStartTime_factory::evnt(jana::JEventLoop *loop, int eventnumber) {
  _data.push_back(new DStartTime(0.0));
  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DStartTime_factory::erun(void) {
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DStartTime_factory::fini(void) {
  return NOERROR;
}


  
