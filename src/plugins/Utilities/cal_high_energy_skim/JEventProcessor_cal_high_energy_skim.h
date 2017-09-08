// $Id$
//
//    File: JEventProcessor_cal_high_energy_skim.h
//

#ifndef _JEventProcessor_cal_high_energy_skim_
#define _JEventProcessor_cal_high_energy_skim_

#include <JANA/JEventProcessor.h>
#include <JANA/JApplication.h>
#include "evio_writer/DEventWriterEVIO.h"

#include <TH1F.h>

#include <vector>

using namespace jana;
using namespace std;

class JEventProcessor_cal_high_energy_skim:public jana::JEventProcessor{
 public:

  JEventProcessor_cal_high_energy_skim();
  ~JEventProcessor_cal_high_energy_skim();
  const char* className(void){return "JEventProcessor_cal_high_energy_skim";}


 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);	
  					///< Called after last event of last event source has been processed.
  double MIN_BCAL_E;
  double MIN_FCAL_E;
 
  //int WRITE_EVIO;
  int MAKE_DIAGNOSTICS;
  int num_epics_events;

  TH1F *h_BCAL_shen;
  TH1F *h_FCAL_shen;

};
#endif // _JEventProcessor_cal_high_energy_skim_

