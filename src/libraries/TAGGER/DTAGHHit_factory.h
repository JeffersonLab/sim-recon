// $Id$
//
//    File: DTAGHHit_factory.h
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluey.phys.uconn.edu)
//

#ifndef _DTAGHHit_factory_
#define _DTAGHHit_factory_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include "DTAGHHit.h"
#include "DTAGHGeometry.h"

class DTAGHHit_factory: public jana::JFactory<DTAGHHit> {
   public:
      DTAGHHit_factory() {};
      ~DTAGHHit_factory() {};

      static const int k_counter_dead = 0;
      static const int k_counter_good = 1;
      static const int k_counter_bad = 2;
      static const int k_counter_noisy = 3;

      // fadc pulse integration window width (samples)
      int fadc_pulse_window_size;

      // overall scale factors
      double fadc_a_scale;  // pixels per fADC pulse integral count
      double fadc_t_scale;  // ns per fADC time count
      double tdc_t_scale;   // ns per F1TDC count
      double t_min;

      // calibration constants stored in row, column format
      double fadc_gains[TAGH_MAX_COUNTER+1];
      double fadc_pedestals[TAGH_MAX_COUNTER+1];
      double fadc_time_offsets[TAGH_MAX_COUNTER+1];
      double tdc_time_offsets[TAGH_MAX_COUNTER+1];
      double counter_quality[TAGH_MAX_COUNTER+1];

      bool load_ccdb_constants(std::string table_name,
                               std::string column_name,
                               double table[TAGH_MAX_COUNTER+1]);
   private:
      jerror_t init(void);                                          ///< Called once at program start
      jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);    ///< Called everytime a new run number is detected
      jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);  ///< Called every event
      jerror_t erun(void);                                          ///< Called everytime run number changes, if brun has been called
      jerror_t fini(void);                                          ///< Called after last event of last event source has been processed
};

#endif // _DTAGHHit_factory_
