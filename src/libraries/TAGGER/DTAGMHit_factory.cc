// $Id$
//
//    File: DTAGMHit_factory.cc
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluex.phys.uconn.edu)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTAGMDigiHit.h"
#include "DTAGMTDCDigiHit.h"
#include "DTAGMGeometry.h"
#include "DTAGMHit_factory.h"
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250Config.h>

using namespace jana;

const int DTAGMHit_factory::k_fiber_good;
const int DTAGMHit_factory::k_fiber_bad;
const int DTAGMHit_factory::k_fiber_noisy;

//------------------
// init
//------------------
jerror_t DTAGMHit_factory::init(void)
{
    DELTA_T_ADC_TDC_MAX = 10.0; // ns
    gPARMS->SetDefaultParameter("TAGMHit:DELTA_T_ADC_TDC_MAX", DELTA_T_ADC_TDC_MAX,
                "Maximum difference in ns between a (calibrated) fADC time and"
                " F1TDC time for them to be matched in a single hit");
    // initialize calibration constants
    fadc_a_scale = 0;
    fadc_t_scale = 0;
    t_base = 0;
    t_tdc_base=0;

    // calibration constants stored in row, column format
    for (int row = 0; row <= TAGM_MAX_ROW; ++row) {
        for (int col = 0; col <= TAGM_MAX_COLUMN; ++col) {
            fadc_gains[row][col] = 0;
            fadc_pedestals[row][col] = 0;
            fadc_time_offsets[row][col] = 0;
            tdc_time_offsets[row][col] = 0;
            fiber_quality[row][col] = 0;
        }
    }

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTAGMHit_factory::brun(jana::JEventLoop *eventLoop, int32_t runnumber)
{
    // Only print messages for one thread whenever run number change
    static pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
    static set<int> runs_announced;
    pthread_mutex_lock(&print_mutex);
    bool print_messages = false;
    if(runs_announced.find(runnumber) == runs_announced.end()){
        print_messages = true;
        runs_announced.insert(runnumber);
    }
    pthread_mutex_unlock(&print_mutex);

    /// set the base conversion scales
    fadc_a_scale    = 1.1;        // pixels per count
    fadc_t_scale    = 0.0625;     // ns per count
    t_base           = 0.;      // ns

    pthread_mutex_unlock(&print_mutex);
    if(print_messages) jout << "In DTAGMHit_factory, loading constants..." << std::endl;

    // load base time offset
    map<string,double> base_time_offset;
    if (eventLoop->GetCalib("/PHOTON_BEAM/microscope/base_time_offset",base_time_offset))
        jout << "Error loading /PHOTON_BEAM/microscope/base_time_offset !" << endl;
    if (base_time_offset.find("TAGM_BASE_TIME_OFFSET") != base_time_offset.end())
        t_base = base_time_offset["TAGM_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get TAGM_BASE_TIME_OFFSET from /PHOTON_BEAM/microscope/base_time_offset !" << endl;
    if (base_time_offset.find("TAGM_TDC_BASE_TIME_OFFSET") != base_time_offset.end())
        t_tdc_base = base_time_offset["TAGM_TDC_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get TAGM_TDC_BASE_TIME_OFFSET from /PHOTON_BEAM/microscope/base_time_offset !" << endl;

    if (load_ccdb_constants("fadc_gains", "gain", fadc_gains) &&
    load_ccdb_constants("fadc_pedestals", "pedestal", fadc_pedestals) &&
    load_ccdb_constants("fadc_time_offsets", "offset", fadc_time_offsets) &&
    load_ccdb_constants("tdc_time_offsets", "offset", tdc_time_offsets) &&
    load_ccdb_constants("fiber_quality", "code", fiber_quality) &&
    load_ccdb_constants("tdc_timewalk_corrections", "c0", tw_c0) &&
    load_ccdb_constants("tdc_timewalk_corrections", "c1", tw_c1) &&
    load_ccdb_constants("tdc_timewalk_corrections", "c2", tw_c2) &&
    load_ccdb_constants("tdc_timewalk_corrections", "threshold", thresh) &&
    load_ccdb_constants("tdc_timewalk_corrections", "pp_0", P_0))
    {
        return NOERROR;
    }
    return UNRECOVERABLE_ERROR;
}

//------------------
// evnt
//------------------
jerror_t DTAGMHit_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    /// Generate DTAGMHit object for each DTAGMDigiHit object.
    /// This is where the first set of calibration constants
    /// is applied to convert from digitzed units into natural
    /// units.
    ///
    /// Note that this code does NOT get called for simulated
    /// data in HDDM format. The HDDM event source will copy
    /// the precalibrated values directly into the _data vector.

    // extract the TAGM geometry
    vector<const DTAGMGeometry*> tagmGeomVect;
    eventLoop->Get( tagmGeomVect );
    if (tagmGeomVect.size() < 1)
        return OBJECT_NOT_AVAILABLE;
    const DTAGMGeometry& tagmGeom = *(tagmGeomVect[0]);

    const DTTabUtilities* locTTabUtilities = NULL;
    loop->GetSingle(locTTabUtilities);

    // First loop over all TAGMDigiHits and make DTAGMHits out of them
    vector<const DTAGMDigiHit*> digihits;
    loop->Get(digihits);
    for (unsigned int i=0; i < digihits.size(); i++) {
        const DTAGMDigiHit *digihit = digihits[i];

        // The following condition signals an error state in the flash algorithm
        // Do not make hits out of these
        const Df250PulsePedestal* PPobj = NULL;
        digihit->GetSingle(PPobj);
        if (PPobj != NULL){
            if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
        }
        else continue;

        // Get pedestal, prefer associated event pedestal if it exists,
        // otherwise, use the average pedestal from CCDB
        double pedestal = fadc_pedestals[digihit->row][digihit->column];

        const Df250PulseIntegral* PIobj = NULL;
        digihit->GetSingle(PIobj);
        if (PIobj != NULL) {
            // the measured pedestal must be scaled by the ratio of the number
            // of samples used to calculate the pedestal and the actual pulse
            // Changed to match D. Lawrence Dec 4 2014 changes
            double single_sample_ped = (double)PIobj->pedestal;
            double nsamples_integral = (double)PIobj->nsamples_integral;
            double nsamples_pedestal = (double)PIobj->nsamples_pedestal;
            pedestal          = single_sample_ped * nsamples_integral/nsamples_pedestal;
        }

        // throw away hits from bad or noisy fibers
        int quality = fiber_quality[digihit->row][digihit->column];
        if (quality == k_fiber_bad || quality == k_fiber_noisy)
            continue;

        // Skip events where fADC algorithm fails
        //if (digihit->pulse_time == 0) continue; // Should already be caught above

        DTAGMHit *hit = new DTAGMHit;
        int row = digihit->row;
        int column = digihit->column;
        hit->row = row;
        hit->column = column;
        double Elow = tagmGeom.getElow(column);
        double Ehigh = tagmGeom.getEhigh(column);
        hit->E = (Elow + Ehigh)/2;
        hit->t = 0;

        // Apply calibration constants
        double P = PPobj->pulse_peak - PPobj->pedestal;
        double A = digihit->pulse_integral;
        double T = digihit->pulse_time;
        A -= pedestal;
        hit->integral = A;
        hit->pulse_peak = P;
        hit->npix_fadc = A * fadc_a_scale * fadc_gains[row][column];
        hit->time_fadc = T * fadc_t_scale - fadc_time_offsets[row][column] + t_base;
        hit->t=hit->time_fadc;
        hit->has_TDC=false;
        hit->has_fADC=true;

        hit->AddAssociatedObject(digihit);
        _data.push_back(hit);
    }

    // Next, loop over TDC hits, matching them to the existing fADC hits
    // where possible and updating their time information. If no match is
    // found, then create a new hit with just the TDC info.
    vector<const DTAGMTDCDigiHit*> tdcdigihits;
    loop->Get(tdcdigihits);
    for (unsigned int i=0; i < tdcdigihits.size(); i++) {
        const DTAGMTDCDigiHit *digihit = tdcdigihits[i];

        // Apply calibration constants here
        int row = digihit->row;
        int column = digihit->column;
        double T = locTTabUtilities->Convert_DigiTimeToNs_F1TDC(digihit) - tdc_time_offsets[row][column] + t_tdc_base;

        // Look for existing hits to see if there is a match
        // or create new one if there is no match
        DTAGMHit *hit = 0;
        for (unsigned int j=0; j < _data.size(); ++j) {
            if (_data[j]->row == row && _data[j]->column == column &&
            fabs(T - _data[j]->time_fadc) < DELTA_T_ADC_TDC_MAX)
          {
            hit = _data[j];
          }
        }
        if (hit == 0) {
            hit = new DTAGMHit;
        hit->row = row;
        hit->column = column;
        double Elow = tagmGeom.getElow(column);
        double Ehigh = tagmGeom.getEhigh(column);
        hit->E = (Elow + Ehigh)/2;
        hit->time_fadc = 0;
        hit->npix_fadc = 0;
        hit->integral = 0;
        hit->pulse_peak = 0;
        hit->has_fADC=false;
        _data.push_back(hit);
        }
        hit->time_tdc=T;
        hit->has_TDC=true;

        // apply time-walk corrections
        double P = hit->pulse_peak;
        // c0 is not used, should it be??
        //double c0 = tw_c0[row][column];
        double c1 = tw_c1[row][column];
        double c2 = tw_c2[row][column];
        double TH = thresh[row][column];
        double pp_0 = P_0[row][column];
        if (P > 0) {
           T -= c1*(pow(P/TH,c2)-pow(pp_0/TH,c2));
        }
        
        hit->t = T;
        
        hit->AddAssociatedObject(digihit);
    }

    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTAGMHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTAGMHit_factory::fini(void)
{
    return NOERROR;
}

//---------------------
// load_ccdb_constants
//---------------------
bool DTAGMHit_factory::load_ccdb_constants(
        std::string table_name,
        std::string column_name,
        double result[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1])
{
    std::vector< std::map<std::string, double> > table;
    std::string ccdb_key = "/PHOTON_BEAM/microscope/" + table_name;
    if (eventLoop->GetCalib(ccdb_key, table))
    {
        jout << "Error loading " << ccdb_key << " from ccdb!" << std::endl;
        return false;
    }
    for (unsigned int i=0; i < table.size(); ++i) {
        int row = (table[i])["row"];
        int col = (table[i])["column"];
        result[row][col] = (table[i])[column_name];
    }
    return true;
}
