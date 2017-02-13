// $Id$
//
//    File: JEventProcessor_BCAL_LED.h
//

#ifndef _JEventProcessor_BCAL_LED_
#define _JEventProcessor_BCAL_LED_

#include <JANA/JEventProcessor.h>


class JEventProcessor_BCAL_LED:public jana::JEventProcessor{
 public:
  JEventProcessor_BCAL_LED();
  ~JEventProcessor_BCAL_LED();
  const char* className(void){return "JEventProcessor_BCAL_LED";}

//  int NOtrig, GTPtrig, FPtrig, FPGTPtrig, trigUS, trigDS, trigCosmic;
//  int low_down_1_counter, low_down_2_counter, low_down_3_counter, low_down_4_counter, low_up_1_counter, low_up_2_counter, low_up_3_counter, low_up_4_counter, high_down_1_counter, high_down_2_counter, high_down_3_counter, high_down_4_counter, high_up_1_counter, high_up_2_counter, high_up_3_counter, high_up_4_counter;
//  int unidentified, ledcounter;
  int adccount1, adccount2, adccount3, nbins;
  double maxnumberofevents;

  int local_eventnum=0;
  

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
};

#endif // _JEventProcessor_BCAL_LED_
