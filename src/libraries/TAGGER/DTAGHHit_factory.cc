// $Id$
//
//    File: DTAGHHit_factory.cc
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluex.phys.uconn.edu)
//
// nsparks moved original factory to DTAGHHit_factory_Calib.cc on Aug 3 2016.

#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include "DTAGHHit_factory.h"
using namespace jana;

inline bool DTAGHHit_SortByID(const DTAGHHit* h1, const DTAGHHit* h2)
{
    if (h1->counter_id == h2->counter_id) return h1->t < h2->t;
    return h1->counter_id < h2->counter_id;
}

//------------------
// init
//------------------
jerror_t DTAGHHit_factory::init(void)
{
    // Set default config. parameters
    MERGE_DOUBLES = true; // Merge double hits?
    gPARMS->SetDefaultParameter("TAGHHit:MERGE_DOUBLES", MERGE_DOUBLES,
    "Merge double hits?");
    DELTA_T_DOUBLES_MAX = 1.0; // ns
    gPARMS->SetDefaultParameter("TAGHHit:DELTA_T_DOUBLES_MAX", DELTA_T_DOUBLES_MAX,
    "Maximum time difference in ns between hits in adjacent counters"
    " for them to be merged into a single hit");
    DELTA_ID_DOUBLES_MAX = 1; // counters
    gPARMS->SetDefaultParameter("TAGHHit:DELTA_ID_DOUBLES_MAX", DELTA_ID_DOUBLES_MAX,
    "Maximum counter id difference of merged hits");
    ID_DOUBLES_MAX = 274;
    gPARMS->SetDefaultParameter("TAGHHit:ID_DOUBLES_MAX", ID_DOUBLES_MAX,
    "Maximum counter id of a double hit");
    USE_SIDEBAND_DOUBLES = false;
    gPARMS->SetDefaultParameter("TAGHHit:USE_SIDEBAND_DOUBLES", USE_SIDEBAND_DOUBLES,
    "Use sideband to estimate accidental coincidences between neighbors?");

    // Setting this flag makes it so that JANA does not delete the objects in _data.
    SetFactoryFlag(NOT_OBJECT_OWNER);

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTAGHHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
    //RF Period
    vector<double> locBeamPeriodVector;
    eventLoop->GetCalib("PHOTON_BEAM/RF/beam_period", locBeamPeriodVector);
    dBeamBunchPeriod = locBeamPeriodVector[0];
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTAGHHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // A scattered (brems.) electron can hit multiple TAGH counters, due
    // to multiple scattering and small geometric overlap between certain counters.
    // These time-coincident hits should be merged to avoid double counting.
    // This factory outputs hits after merging double hits (doubles).

    // Clear hit vector for next event
    _data.clear();

    // Get (calibrated) TAGH hits
    vector<const DTAGHHit*> hits;
    loop->Get(hits, "Calib");

    // Sort TAGH hits by counter id
    sort(hits.begin(),hits.end(),DTAGHHit_SortByID);

    for (size_t i = 0; i < hits.size(); i++) {
        DTAGHHit* hit1 = const_cast<DTAGHHit*>(hits[i]);

        if (!hit1->has_fADC) continue;
        if (!MERGE_DOUBLES || hit1->counter_id > ID_DOUBLES_MAX) {
            _data.push_back(hit1);
            continue;
        }
        if (hit1->is_double) continue;

        for (size_t j = i+1; j < hits.size(); j++) {
            DTAGHHit* hit2 = const_cast<DTAGHHit*>(hits[j]);

            if (!hit2->has_fADC) continue;
            size_t d = abs(hit2->counter_id-hit1->counter_id);
            if (d == 0) continue;
            if (d > DELTA_ID_DOUBLES_MAX) break;

            if (IsDoubleHit(hit2->t-hit1->t)) {
                hit2->is_double = true;
                hit1->is_double = true;
            }
        }
        _data.push_back(hit1);
    }

    return NOERROR;
}

bool DTAGHHit_factory::IsDoubleHit(double tdiff)
{
    if (!USE_SIDEBAND_DOUBLES) {
        return fabs(tdiff) < DELTA_T_DOUBLES_MAX;
    } else {
        return (tdiff > -DELTA_T_DOUBLES_MAX - dBeamBunchPeriod)
        && (tdiff < DELTA_T_DOUBLES_MAX - dBeamBunchPeriod);
    }
}

//------------------
// erun
//------------------
jerror_t DTAGHHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTAGHHit_factory::fini(void)
{
    return NOERROR;
}
