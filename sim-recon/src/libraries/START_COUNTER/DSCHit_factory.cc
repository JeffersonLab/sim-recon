// $Id$
//
//    File: DSCHit_factory.cc
// Created: Tue Aug  6 12:53:32 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
#include <cmath>
using namespace std;

#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>
#include "DSCHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DSCHit_factory::init(void)
{
	DELTA_T_ADC_TDC_MAX = 4.0; // ns
	gPARMS->SetDefaultParameter("SC:DELTA_T_ADC_TDC_MAX", DELTA_T_ADC_TDC_MAX, "Maximum difference in ns between a (calibrated) fADC time and F1TDC time for them to be matched in a single hit");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DSCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// Read in calibration constants (Needs to be done!)
	a_scale    = 2.0E-2/5.2E-5; 
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
jerror_t DSCHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DSCHit object for each DSCDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	// First, make hits out of all fADC250 hits
	vector<const DSCDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DSCDigiHit *digihit = digihits[i];

		DSCHit *hit = new DSCHit;
		hit->sector = digihit->sector;
		
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
	vector<const DSCTDCDigiHit*> tdcdigihits;
	loop->Get(tdcdigihits);
	for(unsigned int i=0; i<tdcdigihits.size(); i++){
		const DSCTDCDigiHit *digihit = tdcdigihits[i];

		// Apply calibration constants here
		double T = (double)digihit->time;
		T = tdc_scale * (T - tdc_offset);
		
		// Look for existing hits to see if there is a match
		// or create new one if there is no match
		DSCHit *hit = FindMatch(digihit->sector, T);
		if(!hit){
			hit = new DSCHit;
			hit->sector = digihit->sector;
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
DSCHit* DSCHit_factory::FindMatch(int sector, double T)
{
	// Loop over existing hits (from fADC) and look for a match
	// in both the sector and the time.
	for(unsigned int i=0; i<_data.size(); i++){
		DSCHit *hit = _data[i];
		
		if(!hit->has_fADC) continue; // only match to fADC hits, not bachelor TDC hits
		if(hit->sector != sector) continue;
		
		double delta_T = fabs(hit->t - T);
		if(delta_T > DELTA_T_ADC_TDC_MAX) continue;
		
		return hit;
	}
	
	return NULL;
}

//------------------
// erun
//------------------
jerror_t DSCHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DSCHit_factory::fini(void)
{
	return NOERROR;
}

