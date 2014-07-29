// $Id$
//
//    File: DTOFHit_factory.cc
// Created: Wed Aug  7 09:30:17 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include "DTOFHit_factory.h"
using namespace jana;

#define TOF_MAX_CHANNELS 176
#define TOF_NUM_PLANES   2 
#define TOF_NUM_BARS    44

static int USE_MC_CALIB = 0;

//------------------
// init
//------------------
jerror_t DTOFHit_factory::init(void)
{
	DELTA_T_ADC_TDC_MAX = 4.0; // ns
	gPARMS->SetDefaultParameter("TOF:DELTA_T_ADC_TDC_MAX", DELTA_T_ADC_TDC_MAX, "Maximum difference in ns between a (calibrated) fADC time and F1TDC time for them to be matched in a single hit");

        // should we use calibrations for simulated data? - this is a temporary workaround
        gPARMS->SetDefaultParameter("DIGI:USEMC",USE_MC_CALIB);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTOFHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
    /*
        // read in geometry information
        vector<const DTOFGeometry*> tofGeomVect;
	eventLoop->Get( tofGeomVect );
	if(tofGeomVect.size()<1)  return OBJECT_NOT_AVAILABLE;
	const DTOFGeometry& tofGeom = *(tofGeomVect[0]);
    */
	/// Set basic conversion constants
	a_scale    = 0.2/5.2E5;
	t_scale    = 0.0625;   // 62.5 ps/count
	tdc_scale    = 0.060;    // 60 ps/count

        /// Read in calibration constants
        vector<double> raw_adc_pedestals;
        vector<double> raw_adc_gains;
        vector<double> raw_adc_offsets;
        vector<double> raw_tdc_offsets;
        vector<double> raw_tdc_scales;

        jout << "In DTOFHit_factory, loading constants..." << endl;

        if(eventLoop->GetCalib("TOF/pedestals", raw_adc_pedestals))
	    jout << "Error loading /TOF/pedestals !" << endl;
        if(eventLoop->GetCalib("TOF/gains", raw_adc_gains))
	    jout << "Error loading /TOF/gains !" << endl;
	if(eventLoop->GetCalib("TOF/timing_scales", raw_tdc_scales))
	    jout << "Error loading /TOF/timing_scales !" << endl;
	if(USE_MC_CALIB>0) {
	    if(eventLoop->GetCalib("TOF/adc_timing_offsets::mc", raw_adc_offsets))
		jout << "Error loading /TOF/adc_timing_offsets !" << endl;
	    if(eventLoop->GetCalib("TOF/timing_offsets::mc", raw_tdc_offsets))
		jout << "Error loading /TOF/timing_offsets !" << endl;
	} else {
	    if(eventLoop->GetCalib("TOF/adc_timing_offsets", raw_adc_offsets))
		jout << "Error loading /TOF/adc_timing_offsets !" << endl;
	    if(eventLoop->GetCalib("TOF/timing_offsets", raw_tdc_offsets))
		jout << "Error loading /TOF/timing_offsets !" << endl;
	}

        FillCalibTable(adc_pedestals, raw_adc_pedestals);
        FillCalibTable(adc_gains, raw_adc_gains);
        FillCalibTable(adc_time_offsets, raw_adc_offsets);
        FillCalibTable(tdc_time_offsets, raw_tdc_offsets);
        FillCalibTable(tdc_scales, raw_tdc_scales);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTOFHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DTOFHit object for each DTOFDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	// First, make hits out of all fADC250 hits
	vector<const DTOFDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DTOFDigiHit *digihit = digihits[i];

		DTOFHit *hit = new DTOFHit;
		hit->plane = digihit->plane;
		hit->bar   = digihit->bar;
		hit->end   = digihit->end;
		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		// TOF constants temporarily disabled
		hit->dE = a_scale * (A - GetConstant(adc_pedestals, digihit));
		hit->t = t_scale * (T - GetConstant(adc_time_offsets, digihit));
		hit->sigma_t = 4.0;    // ns (what is the fADC time resolution?)
		hit->has_fADC = true;
		hit->has_TDC  = false; // will get set to true below if appropriate
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	// Second, loop over TDC hits, matching them to the
	// existing fADC hits where possible and updating
	// their time information. If no match is found, then
	// create a new hit with just the TDC info.
	vector<const DTOFTDCDigiHit*> tdcdigihits;
	loop->Get(tdcdigihits);
	for(unsigned int i=0; i<tdcdigihits.size(); i++){
		const DTOFTDCDigiHit *digihit = tdcdigihits[i];

		// Apply calibration constants here
		double T = (double)digihit->time;
		// TOF constants temporarily disabled
		T = GetConstant(tdc_scales , digihit)
		    * (T - GetConstant(tdc_time_offsets, digihit));
		T = t_scale * T;

		// add in timewalk corrections here
		
		// Look for existing hits to see if there is a match
		// or create new one if there is no match
		DTOFHit *hit = FindMatch(digihit->plane, hit->bar, hit->end, T);  
		if(!hit){
			hit = new DTOFHit;
			hit->plane = digihit->plane;
			hit->bar   = digihit->bar;
			hit->end   = digihit->end;
			hit->dE = 0.0;
			hit->has_fADC = false;

			_data.push_back(hit);
		}
		
		hit->t = T;
		hit->sigma_t = 0.160;    // ns (what is the TOF TDC time resolution?)
		hit->has_TDC = true;
		
		hit->AddAssociatedObject(digihit);
	}

	return NOERROR;
}

//------------------
// FindMatch
//------------------
DTOFHit* DTOFHit_factory::FindMatch(int plane, int bar, int end, double T)
{
	// Loop over existing hits (from fADC) and look for a match
	// in both the sector and the time.
	for(unsigned int i=0; i<_data.size(); i++){
		DTOFHit *hit = _data[i];
		
		if(!hit->has_fADC) continue; // only match to fADC hits, not bachelor TDC hits
		if(hit->plane != plane) continue;
		if(hit->bar != bar) continue;
		if(hit->end != end) continue;
		
		double delta_T = fabs(hit->t - T);
		if(delta_T > DELTA_T_ADC_TDC_MAX) continue;
		
		return hit;
	}
	
	return NULL;
}

//------------------
// erun
//------------------
jerror_t DTOFHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTOFHit_factory::fini(void)
{
	return NOERROR;
}


//------------------
// GetConstant
//------------------
double DTOFHit_factory::GetConstant( tof_digi_constants_t &the_table, 
				     const DTOFDigiHit *the_digihit) 
{
    // we have two ends, indexed as 0/1
    if(the_digihit->end == 0) {
	return the_table[the_digihit->plane][the_digihit->bar].first;
    } else {
	return the_table[the_digihit->plane][the_digihit->bar].second;
    }

}

//------------------
// GetConstant
//------------------
double DTOFHit_factory::GetConstant( tof_digi_constants_t &the_table, 
				     const DTOFTDCDigiHit *the_digihit) 
{
    // we have two ends, indexed as 0/1
    if(the_digihit->end == 0) {
	return the_table[the_digihit->plane][the_digihit->bar].first;
    } else {
	return the_table[the_digihit->plane][the_digihit->bar].second;
    }

}


//------------------
// FillCalibTable
//------------------
void DTOFHit_factory::FillCalibTable(tof_digi_constants_t &table, vector<double> &raw_table)
{
    int channel = 0;

    table.clear();

    for(int plane=0; plane<TOF_NUM_PLANES; plane++) {
	table.push_back( vector< pair<double,double> >(TOF_NUM_BARS) );
	for(int bar=0; bar<TOF_NUM_BARS; bar++) {
	    if( (channel > TOF_MAX_CHANNELS) || (channel+1 > TOF_MAX_CHANNELS) )   // sanity check
		return;

	    table[plane][bar] = pair<double,double>(raw_table[channel],raw_table[channel+1]);
	    
	    channel += 2;
	    
	}
    }
}
