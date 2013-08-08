// $Id$
//
//    File: DFCALHit_factory.cc
// Created: Tue Aug  6 12:23:43 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <FCAL/DFCALDigiHit.h>
#include <FCAL/DFCALGeometry.h>
#include "DFCALHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DFCALHit_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DFCALHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// Read in calibration constants (Needs to be done!)
	a_scale    = 1.0E-4; // 100 keV/count (?)
	a_pedestal = 0.0;
	t_scale    = 4.0;    // 4 ns/count
	t_offset   = 0;
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DFCALHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DFCALHit object for each DFCALDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	// extract the FCAL Geometry (for positionOnFace())
	vector<const DFCALGeometry*> fcalGeomVect;
	eventLoop->Get( fcalGeomVect );
	if(fcalGeomVect.size()<1)return OBJECT_NOT_AVAILABLE;
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);
	
	vector<const DFCALDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DFCALDigiHit *digihit = digihits[i];
		
		DFCALHit *hit = new DFCALHit;
		hit->row    = digihit->row;
		hit->column = digihit->column;
		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		hit->E = a_scale * (A - a_pedestal);
		hit->t = t_scale * (T - t_offset);

		// Get position of blocks on front face. (This should really come from
		// hdgeant directly so the poisitions can be shifted in mcsmear.)
		DVector2 pos = fcalGeom.positionOnFace(hit->row, hit->column);
		
		hit->x = pos.X();
		hit->y = pos.Y();
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DFCALHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DFCALHit_factory::fini(void)
{
	return NOERROR;
}

