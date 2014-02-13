// $Id$
//
//    File: DTOFHit_factory.cc
// Created: Wed Aug  7 09:30:17 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include "DTOFHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DTOFHit_factory::init(void)
{
	DELTA_T_ADC_TDC_MAX = 4.0; // ns
	gPARMS->SetDefaultParameter("TOF:DELTA_T_ADC_TDC_MAX", DELTA_T_ADC_TDC_MAX, "Maximum difference in ns between a (calibrated) fADC time and F1TDC time for them to be matched in a single hit");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTOFHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// Read in calibration constants (Needs to be done!)
	a_scale    = 0.2/5.2E5;
	a_pedestal = 0.0;
	t_scale    = 4.0;    // 4 ns/count
	t_offset   = 0;

	tdc_scale    = 0.060;    // 60 ps/count
	tdc_offset   = 0;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTOFHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DTOFHit object for each DSCDigiHit object.
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
		hit->dE = a_scale * (A - a_pedestal);
		hit->t = t_scale * (T - t_offset);
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
		T = tdc_scale * (T - tdc_offset);
		
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
		hit->sigma_t = 0.160;    // ns (what is the SC TDC time resolution?)
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

