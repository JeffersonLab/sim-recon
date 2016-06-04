// $Id$
//
//    File: JEventProcessor_TRIG_online.h
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TRIG_online_
#define _JEventProcessor_TRIG_online_

#include <map>

#include "TH1I.h"
#include "TH2I.h"

#include <JANA/JEventProcessor.h>

class JEventProcessor_TRIG_online:public jana::JEventProcessor{
 public:
  JEventProcessor_TRIG_online();
  ~JEventProcessor_TRIG_online();
  const char* className(void){return "JEventProcessor_TRIG_online";}


 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  bool timing;
  vector<uint32_t> dTrigBits;

  TH1I* h1trig_trgbits; 

  TH1I* h1trig_tot;
  TH2I* h2trig_fcalVSbcal;
  TH1I* h1trig_fcal_time;
  TH1I* h1trig_bcal_time;
  TH2I* h2trig_tfcalVStbcal;
  TH2I* h2trig_tfcalVSfcal;
  TH2I* h2trig_tbcalVSbcal;
 
  map<uint32_t, TH2I*> h2trigbits_fcalVSbcal;
  map<uint32_t, TH1I*> h1trigbits_fcal_time;
  map<uint32_t, TH1I*> h1trigbits_bcal_time;
  map<uint32_t, TH2I*> h2trigbits_tfcalVStbcal;
  map<uint32_t, TH2I*> h2trigbits_tfcalVSfcal;
  map<uint32_t, TH2I*> h2trigbits_tbcalVSbcal;
};

#endif // _JEventProcessor_TRIG_online_

