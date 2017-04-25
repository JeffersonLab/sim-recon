// $Id$
//
//    File: JEventProcessor_BCAL_LED.h
//

#ifndef _JEventProcessor_BCAL_LED_
#define _JEventProcessor_BCAL_LED_

#include <JANA/JEventProcessor.h>

#include <stdint.h>
#include <vector>
#include "TTree.h"
#include "JEventProcessor_BCAL_LED.h"
#include <JANA/JApplication.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace jana;

#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALUnifiedHit.h"
#include "DAQ/Df250PulseIntegral.h"
#include "DAQ/Df250WindowRawData.h"
#include "TRIGGER/DL1Trigger.h"

#include <TDirectory.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TProfile2D.h>
#include <TStyle.h>

class JEventProcessor_BCAL_LED:public jana::JEventProcessor{
 public:
  JEventProcessor_BCAL_LED();
  ~JEventProcessor_BCAL_LED();
  const char* className(void){return "JEventProcessor_BCAL_LED";}

//  int NOtrig, GTPtrig, FPtrig, FPGTPtrig, trigUS, trigDS, trigCosmic;
//  int low_down_1_counter, low_down_2_counter, low_down_3_counter, low_down_4_counter, low_up_1_counter, low_up_2_counter, low_up_3_counter, low_up_4_counter, high_down_1_counter, high_down_2_counter, high_down_3_counter, high_down_4_counter, high_up_1_counter, high_up_2_counter, high_up_3_counter, high_up_4_counter;
//  int unidentified, ledcounter;
  int adccount, nbins; //adccount2, adccount3, nbins;
  double maxnumberofevents;
  

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  //DO NOT MAKE THESE STATIC GLOBAL EVER AGAIN!!

  // root hist pointers

  //all channels
  TProfile *bcal_peak_vevent = NULL;

  //2 sides
  TProfile *up_peak_vevent = NULL;
  TProfile *down_peak_vevent = NULL;
  //4 columns
  TProfile *column1_peak_vevent = NULL;
  TProfile *column2_peak_vevent = NULL;
  TProfile *column3_peak_vevent = NULL;
  TProfile *column4_peak_vevent = NULL;

  TProfile *column1_up_peak_vevent = NULL;
  TProfile *column2_up_peak_vevent = NULL;
  TProfile *column3_up_peak_vevent = NULL;
  TProfile *column4_up_peak_vevent = NULL;
  TProfile *column1_down_peak_vevent = NULL;
  TProfile *column2_down_peak_vevent = NULL;
  TProfile *column3_down_peak_vevent = NULL;
  TProfile *column4_down_peak_vevent = NULL;


  TProfile *column1_up_peak_vevent1 = NULL;
  TProfile *column1_down_peak_vevent1 = NULL;
  TProfile *column1_up_peak_vevent2 = NULL;
  TProfile *column1_down_peak_vevent2 = NULL;
  TProfile *column1_up_peak_vevent3 = NULL;
  TProfile *column1_down_peak_vevent3 = NULL;
  TProfile *column1_up_peak_vevent4 = NULL;
  TProfile *column1_down_peak_vevent4 = NULL;

  TProfile *column2_up_peak_vevent1 = NULL;
  TProfile *column2_down_peak_vevent1 = NULL;
  TProfile *column2_up_peak_vevent2 = NULL;
  TProfile *column2_down_peak_vevent2 = NULL;
  TProfile *column2_up_peak_vevent3 = NULL;
  TProfile *column2_down_peak_vevent3 = NULL;
  TProfile *column2_up_peak_vevent4 = NULL;
  TProfile *column2_down_peak_vevent4 = NULL;

  TProfile *column3_up_peak_vevent1 = NULL;
  TProfile *column3_down_peak_vevent1 = NULL;
  TProfile *column3_up_peak_vevent2 = NULL;
  TProfile *column3_down_peak_vevent2 = NULL;
  TProfile *column3_up_peak_vevent3 = NULL;
  TProfile *column3_down_peak_vevent3 = NULL;
  TProfile *column3_up_peak_vevent4 = NULL;
  TProfile *column3_down_peak_vevent4 = NULL;

  TProfile *column4_up_peak_vevent1 = NULL;
  TProfile *column4_down_peak_vevent1 = NULL;
  TProfile *column4_up_peak_vevent2 = NULL;
  TProfile *column4_down_peak_vevent2 = NULL;
  TProfile *column4_up_peak_vevent3 = NULL;
  TProfile *column4_down_peak_vevent3 = NULL;
  TProfile *column4_up_peak_vevent4 = NULL;
  TProfile *column4_down_peak_vevent4 = NULL;



  TProfile *low_up_1 = NULL;
  TProfile *low_up_2 = NULL;
  TProfile *low_up_3 = NULL;
  TProfile *low_up_4 = NULL;
  TProfile *low_down_1 = NULL;
  TProfile *low_down_2 = NULL;
  TProfile *low_down_3 = NULL;
  TProfile *low_down_4 = NULL;

  TProfile *high_up_1 = NULL;
  TProfile *high_up_2 = NULL;
  TProfile *high_up_3 = NULL;
  TProfile *high_up_4 = NULL;
  TProfile *high_down_1 = NULL;
  TProfile *high_down_2 = NULL;
  TProfile *high_down_3 = NULL;
  TProfile *high_down_4 = NULL;

};

#endif // _JEventProcessor_BCAL_LED_
