// $Id$
//
//    File: JEventProcessor_pi0fcalskim.h
// Created: Mon Dec  1 14:57:11 EST 2014
// Creator: shepherd (on Linux ifarm1101 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_pi0fcalskim_
#define _JEventProcessor_pi0fcalskim_

#include <JANA/JEventProcessor.h>
#include "evio_writer/DEventWriterEVIO.h"

class JEventProcessor_pi0fcalskim:public jana::JEventProcessor{
 public:
  JEventProcessor_pi0fcalskim();
  ~JEventProcessor_pi0fcalskim();
  const char* className(void){return "JEventProcessor_pi0fcalskim";}

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  const DEventWriterEVIO* dEventWriterEVIO;		

  double MIN_MASS;
  double MAX_MASS;
  double MIN_E;
  double MIN_R;
  double MAX_DT;
  double MAX_ETOT;
  int    MIN_BLOCKS;

};

#endif // _JEventProcessor_pi0fcalskim_

