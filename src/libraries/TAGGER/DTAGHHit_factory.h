// $Id$
//
//    File: DTAGHHit_factory.h
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluey.phys.uconn.edu)
//

#ifndef _DTAGHHit_factory_
#define _DTAGHHit_factory_

#include <JANA/JFactory.h>
#include "DTAGHHit.h"

class DTAGHHit_factory: public jana::JFactory<DTAGHHit> {
    public:
        DTAGHHit_factory() {};
        ~DTAGHHit_factory() {};

        // config. parameters
        bool MERGE_DOUBLES;
        double DELTA_T_DOUBLES_MAX;
        int ID_DOUBLES_MAX;
        bool USE_SIDEBAND_DOUBLES;

        bool IsDoubleHit(double tdiff);
        void EraseHit(vector<DTAGHHit*> &v, DTAGHHit* hit);
        void MergeDoubles(map<int, vector<DTAGHHit*> > hitsById, map<int, vector<DTAGHHit*> > &doublesById);

    private:
        jerror_t init(void);                                          ///< Called once at program start
        jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);    ///< Called everytime a new run number is detected
        jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);  ///< Called every event
        jerror_t erun(void);                                          ///< Called everytime run number changes, if brun has been called
        jerror_t fini(void);                                          ///< Called after last event of last event source has been processed

        double dBeamBunchPeriod;
        void Reset_Data(void);
        vector<DTAGHHit*> dCreatedTAGHHits;
};

#endif // _DTAGHHit_factory_
