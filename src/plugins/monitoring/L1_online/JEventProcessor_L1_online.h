// $Id$
//
//    File: JEventProcessor_L1_online.h
// Created: Fri Mar 20 16:32:04 EDT 2015
//

#ifndef _JEventProcessor_L1_online_
#define _JEventProcessor_L1_online_

#include <JANA/JEventProcessor.h>

class JEventProcessor_L1_online:public jana::JEventProcessor{
public:
    JEventProcessor_L1_online();
    ~JEventProcessor_L1_online();
    const char* className(void){return "JEventProcessor_L1_online";}

private:
    jerror_t init(void); ///< Called once at program start.
    jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber); ///< Called everytime a new run number is detected.
    jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber); ///< Called every event.
    jerror_t erun(void); ///< Called everytime run number changes, provided brun has been called.
    jerror_t fini(void); ///< Called after last event of last event source has been processed.

    unsigned int trig_bit[33], trig_bit_fp[33];

    int fcal_row_mask_min, fcal_row_mask_max, fcal_col_mask_min, fcal_col_mask_max;

    int run_number;

    int fcal_cell_thr;

    int bcal_cell_thr;


};

#endif // _JEventProcessor_L1_online_

