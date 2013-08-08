// $Id$
//
//    File: DCDCHit_factory.cc
// Created: Tue Aug  6 11:29:56 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <CDC/DCDCDigiHit.h>
#include "DCDCHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DCDCHit_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// Read in calibration constants (Needs to be done!)
	a_scale    = 1.0E-4; // 100 keV/count (?)
	a_pedestal = 0.0;
	t_scale    = 8.0;    // 8 ns/count
	t_offset   = 0;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCDCHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DCDCHit object for each DCDCDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	vector<const DCDCDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DCDCDigiHit *digihit = digihits[i];

		DCDCHit *hit = new DCDCHit;
		hit->ring  = digihit->ring;
		hit->straw = digihit->straw;
		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		hit->q = a_scale * (A - a_pedestal);
		hit->t = t_scale * (T - t_offset);
		hit->d = 0.0;
		hit->itrack = -1;
		hit->ptype = 0;
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DCDCHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DCDCHit_factory::fini(void)
{
	return NOERROR;
}

