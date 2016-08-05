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
    MERGE_DOUBLES = true; // Merge in-time hits of adjacent counters?
    gPARMS->SetDefaultParameter("TAGHHit:MERGE_DOUBLES", MERGE_DOUBLES,
    "Merge in-time hits of adjacent counters?");
    DELTA_T_DOUBLES_MAX = 0.5; // ns
    gPARMS->SetDefaultParameter("TAGHHit:DELTA_T_DOUBLES_MAX", DELTA_T_DOUBLES_MAX,
    "Maximum time difference in ns between hits in adjacent counters"
    " for them to be merged into a single hit");
    COUNTER_ID_DOUBLES_MAX = 192;
    gPARMS->SetDefaultParameter("TAGHHit:COUNTER_ID_DOUBLES_MAX", COUNTER_ID_DOUBLES_MAX,
    "Maximum counter id of a double hit");
    USE_SIDEBAND_DOUBLES = false;
    gPARMS->SetDefaultParameter("TAGHHit:USE_SIDEBAND_DOUBLES", USE_SIDEBAND_DOUBLES,
    "Use sideband to estimate accidental coincidences between neighbors");

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
    // A scattered (brems.) electron can hit two adjacent TAGH counters, due
    // to multiple scattering and small geometric overlap between some counters.
    // These neighboring, coincident hits should be merged to avoid double counting.
    // This factory outputs hits after merging double hits (doubles).

    // Get (calibrated) TAGH hits
    vector<const DTAGHHit*> hits;
    loop->Get(hits, "Calib");
    // Copy data, skip hits with no ADC info.
    for (const auto& hit : hits) {
        if (!hit->has_fADC) continue;
        DTAGHHit *h = new DTAGHHit;
        *h = *hit;
        _data.push_back(h);
    }
    // Merge any double hits
    if (MERGE_DOUBLES && _data.size() > 1)
        MergeDoubles();

    return NOERROR;
}

// Is the hit combo a double hit?
bool DTAGHHit_factory::IsDoubleHit(int counter_id_diff, double tdiff) {
    if (!USE_SIDEBAND_DOUBLES) {
        return (abs(counter_id_diff) == 1) && (fabs(tdiff) < DELTA_T_DOUBLES_MAX);
    } else {
        return (abs(counter_id_diff) == 1) && (tdiff > -DELTA_T_DOUBLES_MAX - dBeamBunchPeriod)
        && (tdiff < DELTA_T_DOUBLES_MAX - dBeamBunchPeriod);
    }
}

// Find indices of adjacent, in-time hits
pair<vector<size_t>, vector<size_t> > DTAGHHit_factory::FindDoubles() {
    pair<vector<size_t>, vector<size_t> > indices;
    for (size_t i = 0; i < _data.size()-1; i++) {
        const DTAGHHit* hit1 = _data[i];
        if (hit1->counter_id > COUNTER_ID_DOUBLES_MAX) continue;
        for (size_t j = i+1; j < _data.size(); j++) {
            const DTAGHHit* hit2 = _data[j];
            int counter_id_diff = hit1->counter_id-hit2->counter_id;
            double tdiff = hit1->t-hit2->t;
            if (IsDoubleHit(counter_id_diff,tdiff)) {
                if (counter_id_diff == -1) {
                    indices.first.push_back(i); indices.second.push_back(j);
                } else {
                    indices.first.push_back(j); indices.second.push_back(i);
                }
            }
        }
    }
    return indices;
}

// Is index in list?
bool DTAGHHit_factory::In(vector<size_t> &indices, size_t index) {
    for (const auto& val : indices) {
        if (index == val) return true;
    }
    return false;
}

// Merge adjacent, in-time hits
void DTAGHHit_factory::MergeDoubles() {
    pair<vector<size_t>, vector<size_t> > indices = FindDoubles();
    if (indices.first.size() == 0) return; // Nothing to merge
    vector<DTAGHHit*> new_data;
    for (size_t i = 0; i < _data.size(); i++) {
        if (!In(indices.first,i) && !In(indices.second,i)) new_data.push_back(_data[i]);
    }
    vector<DTAGHHit*> merged_data;
    for (size_t i = 0; i < indices.first.size(); i++) {
        DTAGHHit* hit1 = _data[indices.first[i]];
        const DTAGHHit* hit2 = _data[indices.second[i]];
        hit1->t = 0.5*(hit1->t + hit2->t); hit1->E = 0.5*(hit1->E + hit2->E);
        hit1->is_double = true;
        merged_data.push_back(hit1);
    }
    for (const auto& d : merged_data) new_data.push_back(d);
    _data = new_data;
    if (_data.size() > 1) MergeDoubles(); // Merge any new doubles
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
