// $Id$
//
//    File: DBCALHit_factory.cc
// Created: Tue Aug  6 09:26:13 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <BCAL/DBCALDigiHit.h>
#include "DBCALHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DBCALHit_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBCALHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// Read in calibration constants (Needs to be done!)
	a_scale    = 1.0; // 100 keV/count (?)
	a_pedestal = 0.0;
	t_scale    = 4.0;    // 4 ns/count
	t_offset   = 0;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DBCALHit object for each DBCALDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	vector<const DBCALDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DBCALDigiHit *digihit = digihits[i];

		DBCALHit *hit = new DBCALHit;
		hit->module = digihit->module;
		hit->layer  = digihit->layer;
		hit->sector = digihit->sector;
		hit->end    = digihit->end;
		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		hit->E = a_scale * (A - a_pedestal);
		hit->t = t_scale * (T - t_offset);
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBCALHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBCALHit_factory::fini(void)
{
	return NOERROR;
}

