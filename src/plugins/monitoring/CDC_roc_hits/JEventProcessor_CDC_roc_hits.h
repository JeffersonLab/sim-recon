// $Id$
//
//    File: JEventProcessor_CDC_roc_hits.h
// Created: Wed Oct 22 2014
// Creator: Naomi Jarvis
//

#ifndef _JEventProcessor_CDC_roc_hits_
#define _JEventProcessor_CDC_roc_hits_

#include <JANA/JEventProcessor.h>

#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125PulseIntegral.h"
#include "DAQ/Df125PulsePedestal.h"
#include "DAQ/Df125CDCPulse.h"
#include "DAQ/Df125Config.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TFile.h>


// static root hist pointers for rootspy

  static TH1I *cdc_nevents;

  static TH2D *cdc_hits_roc25;
  static TH2D *cdc_hits_roc26;
  static TH2D *cdc_hits_roc27;
  static TH2D *cdc_hits_roc28;

  static TH2D *cdc_amp_roc25;
  static TH2D *cdc_amp_roc26;
  static TH2D *cdc_amp_roc27;
  static TH2D *cdc_amp_roc28;



class JEventProcessor_CDC_roc_hits:public jana::JEventProcessor{
 public:
  JEventProcessor_CDC_roc_hits();
  ~JEventProcessor_CDC_roc_hits();
  const char* className(void){return "JEventProcessor_CDC_roc_hits";}


 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.


  int TSTART;  // start of window on le_time for time_selected and _t_ histos (units of sample/10)
  int TSTOP;   // end of window on le_time for time_selected and _t_ histos (units of sample/10)

  TH1D *cdc_sumamp_roc25 = NULL;
  TH1D *cdc_sumamp_roc26 = NULL;
  TH1D *cdc_sumamp_roc27 = NULL;
  TH1D *cdc_sumamp_roc28 = NULL;

  TH2D *cdc_netamp_roc25 = NULL;
  TH2D *cdc_netamp_roc26 = NULL;
  TH2D *cdc_netamp_roc27 = NULL;
  TH2D *cdc_netamp_roc28 = NULL;

  TH2D *cdc_time_roc25 = NULL;
  TH2D *cdc_time_roc26 = NULL;
  TH2D *cdc_time_roc27 = NULL;
  TH2D *cdc_time_roc28 = NULL;

  TH2D *cdc_ped_roc25 = NULL;
  TH2D *cdc_ped_roc26 = NULL;
  TH2D *cdc_ped_roc27 = NULL;
  TH2D *cdc_ped_roc28 = NULL;

  TH1D *cdc_time_selected = NULL;

  TH1D *cdc_sumamp_t_roc25 = NULL;
  TH1D *cdc_sumamp_t_roc26 = NULL;
  TH1D *cdc_sumamp_t_roc27 = NULL;
  TH1D *cdc_sumamp_t_roc28 = NULL;

  TH2D *cdc_amp_t_roc25 = NULL;
  TH2D *cdc_amp_t_roc26 = NULL;
  TH2D *cdc_amp_t_roc27 = NULL;
  TH2D *cdc_amp_t_roc28 = NULL;


};

#endif // _JEventProcessor_CDC_roc_hits_

