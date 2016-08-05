// $Id$
//
//    File: DTAGHHit_factory_Calib.h
// Created: Wed Aug  3 12:55:19 EDT 2016
// Creator: nsparks (on Linux cua2.jlab.org 3.10.0-327.22.2.el7.x86_64 x86_64)
//

#ifndef _DTAGHHit_factory_Calib_
#define _DTAGHHit_factory_Calib_

#include <JANA/JFactory.h>
#include "DTAGHHit.h"
#include "DTAGHGeometry.h"

class DTAGHHit_factory_Calib:public jana::JFactory<DTAGHHit>{
    public:
        DTAGHHit_factory_Calib(){};
        ~DTAGHHit_factory_Calib(){};
        const char* Tag(void){return "Calib";}

        static const int k_counter_dead = 0;
        static const int k_counter_good = 1;
        static const int k_counter_bad = 2;
        static const int k_counter_noisy = 3;

        // config. parameters
        double DELTA_T_ADC_TDC_MAX;
        double ADC_THRESHOLD;

        // overall scale factors
        double fadc_a_scale;  // pixels per fADC pulse integral count
        double fadc_t_scale;  // ns per fADC time count
        double t_base;
        double t_tdc_base;

        // calibration constants stored in row, column format
        double fadc_gains[TAGH_MAX_COUNTER+1];
        double fadc_pedestals[TAGH_MAX_COUNTER+1];
        double fadc_time_offsets[TAGH_MAX_COUNTER+1];
        double tdc_time_offsets[TAGH_MAX_COUNTER+1];
        double counter_quality[TAGH_MAX_COUNTER+1];
        double tdc_twalk_c0[TAGH_MAX_COUNTER+1];
        double tdc_twalk_c1[TAGH_MAX_COUNTER+1];
        double tdc_twalk_c2[TAGH_MAX_COUNTER+1];
        double tdc_twalk_c3[TAGH_MAX_COUNTER+1];

        bool load_ccdb_constants(std::string table_name,
            std::string column_name,
            double table[TAGH_MAX_COUNTER+1]);

    private:
        jerror_t init(void);						///< Called once at program start.
        jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
        jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
        jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
        jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DTAGHHit_factory_Calib_
