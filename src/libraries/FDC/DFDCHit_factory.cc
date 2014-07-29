// $Id$
//
//    File: DFDCHit_factory.cc
// Created: Wed Aug  7 11:55:02 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <FDC/DFDCGeometry.h>
#include <FDC/DFDCCathodeDigiHit.h>
#include <FDC/DFDCWireDigiHit.h>
#include "DFDCHit_factory.h"
using namespace jana;


static int USE_MC_CALIB = 0;

//------------------
// init
//------------------
jerror_t DFDCHit_factory::init(void)
{
  	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DFDCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// set the base conversion scales
	a_scale      = 2.4E4/1.3E5;  // cathodes
	t_scale    = 8.0/10.0;    // 8 ns/count and integer time is in 1/10th of sample
	tdc_scale    = 0.115;        // 115 ps/count

	// reset constants
	a_gains.clear();
	a_pedestals.clear();
	timing_offsets.clear();

	// now load them all
	jout << "In DFDCHit_factory, loading constants..." << endl;

	// each FDC package has the same set of constants
	LoadPackageCalibTables(eventLoop,"/FDC/package1");
	LoadPackageCalibTables(eventLoop,"/FDC/package2");
	LoadPackageCalibTables(eventLoop,"/FDC/package3");
	LoadPackageCalibTables(eventLoop,"/FDC/package4");

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DFDCHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DFDCHit object for each DFDCCathodeDigiHit and
	/// each DFDCWireDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	// Make hits out of all DFDCCathodeDigiHit hits
	vector<const DFDCCathodeDigiHit*> cathodedigihits;
	loop->Get(cathodedigihits);
	for(unsigned int i=0; i<cathodedigihits.size(); i++){
		const DFDCCathodeDigiHit *digihit = cathodedigihits[i];
		
		// The translation table has:
		// ---------------------------------------------------
		// package : 1-4
		// chamber : 1-6
		// view    : 1(="U") or 3(="D")
		// strip   : 1-192
		//
		//
		// The FDCHit class has 6 indexes which are derived
		// from these and contain some redundancy. They are:
		// ---------------------------------------------------
		// layer   : 1(V), 2(X), or 3(U)
		// module  : 1 through 8, 1 module = 3 detection layers
		// element : wire or strip number
		// plane   : for cathodes only: u(3) or v(1) plane, u@+45,v@-45 
		// gPlane  : 1 through 72
		// gLayer  : 1 through 24
		
		DFDCHit *hit = new DFDCHit;
		hit->layer   = digihit->view;
		hit->gLayer  = digihit->chamber + 6*(digihit->package - 1);
		hit->gPlane  = hit->layer + 3*(hit->gLayer - 1);
		hit->module  = 1 + (hit->gLayer-1)/3;
		hit->element = digihit->strip;
		hit->plane   = digihit->view; // "plane" is apparently the same as "layer"
		hit->r       = DFDCGeometry::getWireR(hit);
		hit->d       = 0.0; // MC data only
		hit->type    = digihit->strip_type; // n.b. DEventSourceHDDM hardwires this as "1" for cathodes!
		hit->itrack  = -1; // MC data only
		hit->ptype   = 0;// MC data only

		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		hit->q = a_scale * A;
		hit->t = t_scale * T;

		hit->q = a_scale * a_gains[hit->gPlane-1][hit->element-1] 
		    * (A - a_pedestals[hit->gPlane-1][hit->element-1] );
		hit->t = t_scale * (T - timing_offsets[hit->gPlane-1][hit->element-1]);
		//hit->t = t_scale * (T - adc_timing_offsets[hit->gPlane][hit->element]);
		
		//cerr << "FDC hitL  plane = " << hit->gPlane << "  element = " << hit->element << endl;
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	// Make hits out of all DFDCWireDigiHit hits
	vector<const DFDCWireDigiHit*> wiredigihits;
	loop->Get(wiredigihits);
	for(unsigned int i=0; i<wiredigihits.size(); i++){
		const DFDCWireDigiHit *digihit = wiredigihits[i];
		
		// The translation table has:
		// ---------------------------------------------------
		// package : 1-4
		// chamber : 1-6
		// wire    : 1-96
		//
		//
		// The FDCHit class has 6 indexes which are derived
		// from these and contain some redundancy. They are:
		// ---------------------------------------------------
		// layer   : 1(V), 2(X), or 3(U)
		// module  : 1 through 8, 1 module = 3 detection layers
		// element : wire or strip number
		// plane   : for cathodes only: u(3) or v(1) plane, u@+45,v@-45 
		// gPlane  : 1 through 72
		// gLayer  : 1 through 24
		
		DFDCHit *hit = new DFDCHit;
		hit->layer   = 2;                  // wire is always in "X" layer
		hit->gLayer  = digihit->chamber + 6*(digihit->package - 1);
		hit->gPlane  = hit->layer + 3*(hit->gLayer - 1);
		hit->module  = 1 + (hit->gLayer-1)/3;
		hit->element = digihit->wire;
		hit->plane   = hit->layer;         // "plane" is apparently the same as "layer"
		hit->r       = DFDCGeometry::getWireR(hit);
		hit->d       = 0.0;                // MC data only
		hit->type    = DFDCHit::AnodeWire; // n.b. DEventSourceHDDM hardwires this as "0" for anodes!
		hit->itrack  = -1;                 // MC data only
		hit->ptype   = 0;                  // MC data only

		// Apply calibration constants here
		double T = (double)digihit->time;
		T = tdc_scale * (T - timing_offsets[hit->gPlane-1][hit->element-1]);
		//T = tdc_scale * T;
		hit->q = 0.0; // no charge measured for wires in FDC
		hit->t = T;
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DFDCHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DFDCHit_factory::fini(void)
{
	return NOERROR;
}


//------------------
// LoadPackageCalibTables
//------------------
void DFDCHit_factory::LoadPackageCalibTables(jana::JEventLoop *eventLoop, string ccdb_prefix)
{
    vector< vector<double> >  new_gains, new_pedestals, new_strip_t0s, new_wire_t0s;
    
    if(eventLoop->GetCalib(ccdb_prefix+"/strip_gains", new_gains))
	cout << "Error loading "+ccdb_prefix+"/strip_gains !" << endl;
    if(eventLoop->GetCalib(ccdb_prefix+"/strip_pedestals", new_pedestals))
	cout << "Error loading "+ccdb_prefix+"/strip_pedestals !" << endl;
    if(eventLoop->GetCalib(ccdb_prefix+"/strip_timing_offsets", new_strip_t0s))
	cout << "Error loading "+ccdb_prefix+"/strip_timing_offsets!" << endl;
    if(eventLoop->GetCalib(ccdb_prefix+"/wire_timing_offsets", new_wire_t0s))
	cout << "Error loading "+ccdb_prefix+"/wire_timing_offsets!" << endl;

    for(int nchamber=0; nchamber<6; nchamber++) {

	// load ADC gains (only for cathode strips)
	a_gains.push_back( new_gains[2*nchamber] );
	a_gains.push_back( vector<double>() );
	a_gains.push_back( new_gains[2*nchamber+1] );

	// load ADC pedestals (only for cathode strips)
	a_pedestals.push_back( new_pedestals[2*nchamber] );
	a_pedestals.push_back( vector<double>() );
	a_pedestals.push_back( new_pedestals[2*nchamber+1] );

	// load t0's for strips and wires
	timing_offsets.push_back( new_strip_t0s[2*nchamber] );
	timing_offsets.push_back( new_wire_t0s[nchamber] );
	timing_offsets.push_back( new_strip_t0s[2*nchamber+1] );
	
    }
}

