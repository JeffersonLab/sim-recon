// $Id$
//
//    File: DTAGHHit_factory.cc
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluex.phys.uconn.edu)
//
// nsparks moved original factory to DTAGHHit_factory_Calib.cc on Aug 3 2016.

#include <iostream>
#include <iomanip>
using namespace std;

#include "DTAGHHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DTAGHHit_factory::init(void)
{
    // Set default config. parameters
    MERGE_DOUBLES = true; // Merge double hits?
    gPARMS->SetDefaultParameter("TAGHHit:MERGE_DOUBLES", MERGE_DOUBLES,
    "Merge double hits?");
    DELTA_T_DOUBLES_MAX = 0.8; // ns
    gPARMS->SetDefaultParameter("TAGHHit:DELTA_T_DOUBLES_MAX", DELTA_T_DOUBLES_MAX,
    "Maximum time difference in ns between hits in adjacent counters"
    " for them to be merged into a single hit");
    ID_DOUBLES_MAX = 274; // 192 is last counter with an overlap in energy-boundary table
    gPARMS->SetDefaultParameter("TAGHHit:ID_DOUBLES_MAX", ID_DOUBLES_MAX,
    "Maximum counter id of a double hit");
    USE_SIDEBAND_DOUBLES = false;
    gPARMS->SetDefaultParameter("TAGHHit:USE_SIDEBAND_DOUBLES", USE_SIDEBAND_DOUBLES,
    "Use sideband to estimate accidental coincidences between neighbors?");

    //Setting this flag makes it so that JANA does not delete the objects in _data. This factory will manage this memory.
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

    // Free memory allocated for DTAGHHit pointers in previous event
    Reset_Data();

    // Get (calibrated) TAGH hits
    vector<const DTAGHHit*> hits;
    loop->Get(hits, "Calib");

    // Sort TAGH hits by counter id by putting them in a map
    map<int, vector<DTAGHHit*> > hitsById;
    for (auto&& hit : hits) {
        if (!hit->has_fADC) continue; // Skip hits that have no ADC info.
        hitsById[hit->counter_id].push_back(const_cast<DTAGHHit*>(hit));
    }

    // Merge double hits
    map<int, vector<DTAGHHit*> > doublesById;
    if (MERGE_DOUBLES && hits.size() > 1)
        MergeDoubles(hitsById, doublesById);

    // Add double hits to _data
    for (auto&& p : doublesById) {
        for (auto&& h : p.second) {
            _data.push_back(h);
        }
    }

    // Add single-counter hits to _data
    for (auto&& p : hitsById) {
        for (auto&& h : p.second) {
            if (!h->is_double) _data.push_back(h);
        }
    }

    return NOERROR;
}

bool DTAGHHit_factory::IsDoubleHit(double tdiff) {
    if (!USE_SIDEBAND_DOUBLES) {
        return fabs(tdiff) < DELTA_T_DOUBLES_MAX;
    } else {
        return (tdiff > -DELTA_T_DOUBLES_MAX - dBeamBunchPeriod)
        && (tdiff < DELTA_T_DOUBLES_MAX - dBeamBunchPeriod);
    }
}

void DTAGHHit_factory::MergeDoubles(map<int, vector<DTAGHHit*> > hitsById, map<int, vector<DTAGHHit*> > &doublesById) {
    int prev_id = -1; bool has_doubles = false;
    vector<DTAGHHit*> prev_hits;
    for (auto&& p : hitsById) {
        int id = p.first;
        if (id > ID_DOUBLES_MAX) continue;
        if (id - prev_id == 1) {
            for (auto&& h1 : prev_hits) {
                for (auto&& h2 : p.second) {
                    if (IsDoubleHit(h1->t-h2->t)) {
                        has_doubles = true;
                        if (h1->is_double && h2->is_double) {
                            EraseHit(doublesById[h1->counter_id], h1);
                            EraseHit(doublesById[h2->counter_id], h2);
                        }
                        h1->is_double = true; h2->is_double = true;
                        DTAGHHit *h = new DTAGHHit;
                        dCreatedTAGHHits.push_back(h);
                        *h = *h1;
                        h->t = 0.5*(h1->t + h2->t); h->E = 0.5*(h1->E + h2->E);
                        doublesById[h->counter_id].push_back(h);
                    }
                }
            }
        }
        prev_id = id;
        prev_hits = p.second;
    } // Merge any new double hits
    if (has_doubles) MergeDoubles(doublesById, doublesById);
}

void DTAGHHit_factory::EraseHit(vector<DTAGHHit*> &v, DTAGHHit* hit) {
    int index = -1; bool flag = false;
    for (auto&& i : v) {
        index++;
        if ((i->t == hit->t) && (i->E == hit->E)) {
            flag = true; break;
        }
    }
    if (index >= 0 && flag) v.erase(v.begin() + index);
}

void DTAGHHit_factory::Reset_Data(void)
{
    //delete objects that this factory created (since the NOT_OBJECT_OWNER flag is set)
    for (auto&& hit : dCreatedTAGHHits)
        delete hit;
    _data.clear();
    dCreatedTAGHHits.clear();
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
