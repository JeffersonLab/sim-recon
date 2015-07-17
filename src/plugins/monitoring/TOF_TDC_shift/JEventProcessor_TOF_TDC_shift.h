// $Id$
//
//    File: JEventProcessor_TOF_TDC_shift.h
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TOF_TDC_shift_
#define _JEventProcessor_TOF_TDC_shift_

#include <JANA/JEventProcessor.h>

#include <TDirectory.h>
#include <TH2I.h>
#include "TTree.h"

using namespace std;

class JEventProcessor_TOF_TDC_shift:public jana::JEventProcessor{
 public:
  JEventProcessor_TOF_TDC_shift();
  ~JEventProcessor_TOF_TDC_shift();
  const char* className(void){return "JEventProcessor_TOF_TDC_shift";}

  static const Int_t NMAX = 200;

  /*
  // Number of hits
  Int_t nTOF;
  Int_t nTOF_TDC;
  // ADC times
  Float_t tof_adc_time[NMAX];
  // TDC times
  Float_t tof_tdc_time[NMAX];

  // DCODAROCInfo
  Int_t nDCODAROCInfo;
  static const Int_t NROCINFO = 96;
  Int_t rocid[NROCINFO];
  ULong64_t rocTime[NROCINFO];
  */

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  TH2I *hrocTimeRemainder_AdcTdcTimeDiff;

  ofstream OUTPUT;
};

#endif // _JEventProcessor_TOF_TDC_shift_

