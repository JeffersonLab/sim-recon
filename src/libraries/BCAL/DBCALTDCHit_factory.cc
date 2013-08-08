// $Id$
//
//    File: DBCALTDCHit_factory.cc
// Created: Tue Aug  6 11:04:11 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <JANA/JEventLoop.h>
#include <BCAL/DBCALTDCDigiHit.h>
#include <BCAL/DBCALTDCHit_factory.h>
using namespace jana;

//------------------
// init
//------------------
jerror_t DBCALTDCHit_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBCALTDCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// Read in calibration constants (Needs to be done!)
	t_scale    = 0.060;    // 60 ps/count
	t_offset   = 0;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALTDCHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DBCALTDCHit object for each DBCALTDCDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	vector<const DBCALTDCDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DBCALTDCDigiHit *digihit = digihits[i];

		DBCALTDCHit *hit = new DBCALTDCHit;
		hit->module = digihit->module;
		hit->layer  = digihit->layer;
		hit->sector = digihit->sector;
		hit->end    = digihit->end;
		
		// Apply calibration constants here
		double T = (double)digihit->time;
		hit->t = t_scale * (T - t_offset);
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBCALTDCHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBCALTDCHit_factory::fini(void)
{
	return NOERROR;
}

