// $Id$
//
//    File: JEventProcessor_pi0bcalskim.h
// Created: Mon feb 6 15:46:00 EST 2015
// Creator: wmcginle(on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_pi0bcalskim_
#define _JEventProcessor_pi0bcalskim_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include "evio_writer/DEventWriterEVIO.h"


#include <vector>

using namespace jana;
using namespace std;

class JEventProcessor_pi0bcalskim:public jana::JEventProcessor{
 public:

  JEventProcessor_pi0bcalskim();
  ~JEventProcessor_pi0bcalskim();
  const char* className(void){return "JEventProcessor_pi0bcalskim";}


 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);	
  					///< Called after last event of last event source has been processed.
  double MIN_SH1_E;
  double MIN_SH2_E;

  const DEventWriterEVIO* dEventWriterEVIO;		

 
  int WRITE_EVIO;
  int num_epics_events;

};
#endif // _JEventProcessor_pi0bcalskim_

