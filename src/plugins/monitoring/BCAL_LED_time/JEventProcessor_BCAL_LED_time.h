// $Id$
//
//    File: JEventProcessor_BCAL_LED_time.h
//

#ifndef _JEventProcessor_BCAL_LED_time_
#define _JEventProcessor_BCAL_LED_time_

#include <JANA/JEventProcessor.h>

#include <stdint.h>
#include <vector>
#include "TTree.h"
#include "JEventProcessor_BCAL_LED_time.h"
#include <JANA/JApplication.h>
#include <iostream>
#include <fstream>

using namespace std;
using namespace jana;

#include "BCAL/DBCALDigiHit.h"
#include "BCAL/DBCALHit.h"
#include "BCAL/DBCALPoint.h"
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

#include "ANALYSIS/DTreeInterface.h"


class JEventProcessor_BCAL_LED_time:public jana::JEventProcessor{
 public:
  JEventProcessor_BCAL_LED_time();
  ~JEventProcessor_BCAL_LED_time();
  const char* className(void){return "JEventProcessor_BCAL_LED_time";}

//  int NOtrig, GTPtrig, FPtrig, FPGTPtrig, trigUS, trigDS, trigCosmic;
//  int low_down_1_counter, low_down_2_counter, low_down_3_counter, low_down_4_counter, low_up_1_counter, low_up_2_counter, low_up_3_counter, low_up_4_counter, high_down_1_counter, high_down_2_counter, high_down_3_counter, high_down_4_counter, high_up_1_counter, high_up_2_counter, high_up_3_counter, high_up_4_counter;
//  int unidentified, ledcounter;
  int adccount1, adccount2, adccount3, nbins;
  double maxnumberofevents;
  

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.



// root hist pointers

//all channels
TProfile *bcal_time_vevent = NULL;

//2 sides
TProfile *up_time_vevent = NULL;
TProfile *down_time_vevent = NULL;
//4 columns
TProfile *column1_time_vevent = NULL;
TProfile *column2_time_vevent = NULL;
TProfile *column3_time_vevent = NULL;
TProfile *column4_time_vevent = NULL;

TProfile *column1_up_peak_vevent = NULL;
TProfile *column2_up_peak_vevent = NULL;
TProfile *column3_up_peak_vevent = NULL;
TProfile *column4_up_peak_vevent = NULL;
TProfile *column1_down_peak_vevent = NULL;
TProfile *column2_down_peak_vevent = NULL;
TProfile *column3_down_peak_vevent = NULL;
TProfile *column4_down_peak_vevent = NULL;


TProfile *column1_up_time_vevent1 = NULL;
TProfile *column1_down_time_vevent1 = NULL;
TProfile *column1_up_time_vevent2 = NULL;
TProfile *column1_down_time_vevent2 = NULL;
TProfile *column1_up_time_vevent3 = NULL;
TProfile *column1_down_time_vevent3 = NULL;
TProfile *column1_up_time_vevent4 = NULL;
TProfile *column1_down_time_vevent4 = NULL;

TProfile *column2_up_time_vevent1 = NULL;
TProfile *column2_down_time_vevent1 = NULL;
TProfile *column2_up_time_vevent2 = NULL;
TProfile *column2_down_time_vevent2 = NULL;
TProfile *column2_up_time_vevent3 = NULL;
TProfile *column2_down_time_vevent3 = NULL;
TProfile *column2_up_time_vevent4 = NULL;
TProfile *column2_down_time_vevent4 = NULL;

TProfile *column3_up_time_vevent1 = NULL;
TProfile *column3_down_time_vevent1 = NULL;
TProfile *column3_up_time_vevent2 = NULL;
TProfile *column3_down_time_vevent2 = NULL;
TProfile *column3_up_time_vevent3 = NULL;
TProfile *column3_down_time_vevent3 = NULL;
TProfile *column3_up_time_vevent4 = NULL;
TProfile *column3_down_time_vevent4 = NULL;

TProfile *column4_up_time_vevent1 = NULL;
TProfile *column4_down_time_vevent1 = NULL;
TProfile *column4_up_time_vevent2 = NULL;
TProfile *column4_down_time_vevent2 = NULL;
TProfile *column4_up_time_vevent3 = NULL;
TProfile *column4_down_time_vevent3 = NULL;
TProfile *column4_up_time_vevent4 = NULL;
TProfile *column4_down_time_vevent4 = NULL;
    
TProfile *low_up_1 = NULL;
TProfile *low_up_2 = NULL;
TProfile *low_up_3 = NULL;
TProfile *low_up_4 = NULL;
TProfile *low_down_1 = NULL;
TProfile *low_down_2 = NULL;
TProfile *low_down_3 = NULL;
TProfile *low_down_4 = NULL;

TProfile *low_up = NULL;
TProfile *low_down = NULL;

TProfile *high_up_1 = NULL;
TProfile *high_up_2 = NULL;
TProfile *high_up_3 = NULL;
TProfile *high_up_4 = NULL;
TProfile *high_down_1 = NULL;
TProfile *high_down_2 = NULL;
TProfile *high_down_3 = NULL;
TProfile *high_down_4 = NULL;

TProfile *high_up = NULL;
TProfile *high_down = NULL;

// Histograms added by Elton for z distributions

TProfile* h2_ledboth_Tall_vs_event = NULL;
TProfile* h2_ledboth_sector_vs_event = NULL;
TH1I* h1_ledall_layer = NULL;
TH1I* h1_led0_layer = NULL;

TH1I* h1_ledup_z_all = NULL;;
TH2I* h2_ledup_z_vs_cellid = NULL;
TH1I* h1_ledup_layer = NULL;
TH1I* h1_ledup_sector = NULL;
TH1I* h1_ledup_sector_config = NULL;
TH1I* h1_ledup_Tdiff_all = NULL;
TH1I* h1_ledup_Tup_all = NULL;
TH1I* h1_ledup_Tdown_all = NULL;
TH2I* h2_ledup_Tup_vs_z = NULL;
TH2I* h2_ledup_Tdown_vs_z = NULL;
TProfile* h2_ledup_Tup_vs_event = NULL;
TProfile* h2_ledup_Tdown_vs_event = NULL;
TProfile* h2_ledup_Tall_vs_event = NULL;
TProfile* h2_ledup_sector_vs_event = NULL;

TH1I* h1_leddown_z_all = NULL;
TH2I* h2_leddown_z_vs_cellid = NULL;
TH1I* h1_leddown_layer = NULL;
TH1I* h1_leddown_sector = NULL;
TH1I* h1_leddown_sector_config = NULL;
TH1I* h1_leddown_Tdiff_all = NULL;
TH1I* h1_leddown_Tup_all = NULL;
TH1I* h1_leddown_Tdown_all = NULL;
TH2I* h2_leddown_Tup_vs_z = NULL;
TH2I* h2_leddown_Tdown_vs_z = NULL;
TProfile* h2_leddown_Tup_vs_event = NULL;
TProfile* h2_leddown_Tdown_vs_event = NULL;
TProfile* h2_leddown_Tall_vs_event = NULL;
TProfile* h2_leddown_sector_vs_event = NULL;

};

#endif // _JEventProcessor_BCAL_LED_time_
