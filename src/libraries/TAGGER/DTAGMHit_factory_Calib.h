// $Id$
//
//    File: DTAGMHit_factory_Calib.h
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluey.phys.uconn.edu)
//

#ifndef _DTAGMHit_factory_Calib_
#define _DTAGMHit_factory_Calib_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include "TTAB/DTTabUtilities.h"

#include "DTAGMHit.h"
#include "DTAGMGeometry.h"

class DTAGMHit_factory_Calib: public jana::JFactory<DTAGMHit> {
   public:
      DTAGMHit_factory_Calib() {};
      ~DTAGMHit_factory_Calib() {};
      const char* Tag(void){return "Calib";}

      static const int k_fiber_dead = 0;
      static const int k_fiber_good = 1;
      static const int k_fiber_bad = 2;
      static const int k_fiber_noisy = 3;

      // config. parameter
      double DELTA_T_ADC_TDC_MAX; 
      //int USE_ADC, PEAK_CUT;
      int USE_ADC, CUT_FACTOR;

      // overall scale factors
      double fadc_a_scale;  // pixels per fADC pulse integral count
      double fadc_t_scale;  // ns per fADC time count
      double t_base;
      double t_tdc_base;

      // calibration constants stored in row, column format
      double fadc_gains[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double fadc_pedestals[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double fadc_time_offsets[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double tdc_time_offsets[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double fiber_quality[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double tw_c0[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double tw_c1[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double tw_c2[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double tw_c3[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double ref[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
      double int_cuts[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];

      bool load_ccdb_constants(std::string table_name,
                               std::string column_name,
                               double table[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1]);
   private:
      jerror_t init(void);                                          ///< Called once at program start
      jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);    ///< Called everytime a new run number is detected
      jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);  ///< Called every event
      jerror_t erun(void);                                          ///< Called everytime run number changes, if brun has been called
      jerror_t fini(void);                                          ///< Called after last event of last event source has been processed

      bool CHECK_FADC_ERRORS;
};

#endif // _DTAGMHit_factory_Calib_
