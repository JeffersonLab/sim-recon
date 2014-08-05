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


//------------------
// init
//------------------
jerror_t DFDCHit_factory::init(void)
{
        a_scale = 0.;
	t_scale = 0.;
	tdc_scale = 0.;

  	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DFDCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	/// set the base conversion scales
	a_scale      = 2.4E4/1.3E5;  // cathodes
	t_scale      = 8.0/10.0;     // 8 ns/count and integer time is in 1/10th of sample
	tdc_scale    = 0.115;        // 115 ps/count

	// reset constants tables
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

	// Verify that the right number of layers were loaded
	char str[256];
	if(a_gains.size() != FDC_NUM_PLANES) {
		sprintf(str, "Bad # of planes for FDC gains from CCDB! CCDB=%zu , should be %d", 
			a_gains.size(), FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	if(a_pedestals.size() != FDC_NUM_PLANES) {
		sprintf(str, "Bad # of planes for FDC pedestals from CCDB! CCDB=%zu , should be %d", 
			a_pedestals.size(), FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	if(timing_offsets.size() != FDC_NUM_PLANES) {
		sprintf(str, "Bad # of planes for FDC timing offsets from CCDB! CCDB=%zu , should be %d", 
			timing_offsets.size(), FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}

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
	char str[256];

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

		// Make sure gPlane and stripare in valid range
		if( (hit->gPlane < 1) || (hit->gPlane > FDC_NUM_PLANES)) {
			sprintf(str, "DFDCDigiHit plane out of range! gPlane=%d (should be 1-%d)", hit->gPlane, FDC_NUM_PLANES);
			throw JException(str);
		}
		if( (digihit->strip < 1) || (digihit->strip > STRIPS_PER_PLANE)) {
			sprintf(str, "DFDCDigiHit straw out of range! strip=%d for plane=%d (should be 1-%d)", digihit->strip, hit->gPlane, STRIPS_PER_PLANE);
			throw JException(str);
		}

		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;

		hit->q = a_scale * a_gains[hit->gPlane-1][hit->element-1] 
		    * (A - a_pedestals[hit->gPlane-1][hit->element-1] );
		hit->t = t_scale * (T - timing_offsets[hit->gPlane-1][hit->element-1]);
		
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

		// Make sure gPlane and wire are in valid range
		if( (hit->gPlane < 1) || (hit->gPlane > FDC_NUM_PLANES)) {
			sprintf(str, "DFDCDigiHit plane out of range! gPlane=%d (should be 1-%d)", hit->gPlane, FDC_NUM_PLANES);
			throw JException(str);
		}
		if( (digihit->wire < 1) || (digihit->wire > WIRES_PER_PLANE)) {
			sprintf(str, "DFDCDigiHit straw out of range! wire=%d for plane=%d (should be 1-%d)", digihit->wire, hit->gPlane, WIRES_PER_PLANE);
			throw JException(str);
		}

		// Apply calibration constants here
		double T = (double)digihit->time;
		T = tdc_scale * (T - timing_offsets[hit->gPlane-1][hit->element-1]);
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
    char str[256];
    
    if(eventLoop->GetCalib(ccdb_prefix+"/strip_gains", new_gains))
	cout << "Error loading "+ccdb_prefix+"/strip_gains !" << endl;
    if(eventLoop->GetCalib(ccdb_prefix+"/strip_pedestals", new_pedestals))
	cout << "Error loading "+ccdb_prefix+"/strip_pedestals !" << endl;
    if(eventLoop->GetCalib(ccdb_prefix+"/strip_timing_offsets", new_strip_t0s))
	cout << "Error loading "+ccdb_prefix+"/strip_timing_offsets!" << endl;
    if(eventLoop->GetCalib(ccdb_prefix+"/wire_timing_offsets", new_wire_t0s))
	cout << "Error loading "+ccdb_prefix+"/wire_timing_offsets!" << endl;

    for(int nchamber=0; nchamber<6; nchamber++) {

	// check the size of table rows
	if(new_gains[2*nchamber].size() != STRIPS_PER_PLANE) {
	    sprintf(str, "Bad # of strips for FDC gain from CCDB! CCDB=%zu , should be %d", new_gains[2*nchamber].size(), STRIPS_PER_PLANE);
	    cerr << str << endl;
	    throw JException(str);
	}
	if(new_gains[2*nchamber+1].size() != STRIPS_PER_PLANE) {
	    sprintf(str, "Bad # of strips for FDC gain from CCDB! CCDB=%zu , should be %d", new_gains[2*nchamber+1].size(), STRIPS_PER_PLANE);
	    cerr << str << endl;
	    throw JException(str);
	}
	if(new_pedestals[2*nchamber].size() != STRIPS_PER_PLANE) {
	    sprintf(str, "Bad # of strips for FDC pedestals from CCDB! CCDB=%zu , should be %d", new_pedestals[2*nchamber].size(), STRIPS_PER_PLANE);
	    cerr << str << endl;
	    throw JException(str);
	}
	if(new_pedestals[2*nchamber+1].size() != STRIPS_PER_PLANE) {
	    sprintf(str, "Bad # of strips for FDC pedestals from CCDB! CCDB=%zu , should be %d", new_pedestals[2*nchamber+1].size(), STRIPS_PER_PLANE);
	    cerr << str << endl;
	    throw JException(str);
	}
	if(new_strip_t0s[2*nchamber].size() != STRIPS_PER_PLANE) {
	    sprintf(str, "Bad # of strips for FDC timing offsets from CCDB! CCDB=%zu , should be %d", new_strip_t0s[2*nchamber].size(), STRIPS_PER_PLANE);
	    cerr << str << endl;
	    throw JException(str);
	}
	if(new_strip_t0s[2*nchamber+1].size() != STRIPS_PER_PLANE) {
	    sprintf(str, "Bad # of strips for FDC timing offsets from CCDB! CCDB=%zu , should be %d", new_strip_t0s[2*nchamber+1].size(), STRIPS_PER_PLANE);
	    cerr << str << endl;
	    throw JException(str);
	}
	if(new_wire_t0s[nchamber].size() != WIRES_PER_PLANE) {
	    sprintf(str, "Bad # of wires for FDC timing offsets from CCDB! CCDB=%zu , should be %d", new_wire_t0s[2*nchamber].size(), WIRES_PER_PLANE);
	    cerr << str << endl;
	    throw JException(str);
	}


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

//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DFDCHit_factory::GetConstant(const fdc_digi_constants_t &the_table,
					  const int in_gPlane, const int in_element) const {
	
	char str[256];
	
	if( (in_gPlane <= 0) || (static_cast<unsigned int>(in_gPlane) > FDC_NUM_PLANES)) {
	        sprintf(str, "Bad gPlane # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", in_gPlane, FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	// strip and wire planes have different numbers of elements
	if( (in_element <= 0) || (static_cast<unsigned int>(in_element) > the_table[in_gPlane].size())) {
	        sprintf(str, "Bad element # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %zu", in_element, the_table[in_gPlane].size());
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[in_gPlane-1][in_element-1];
}

const double DFDCHit_factory::GetConstant(const fdc_digi_constants_t &the_table,
				    const DFDCCathodeDigiHit *in_digihit) const {

	char str[256];

	int gLayer = in_digihit->chamber + 6*(in_digihit->package - 1);
	int gPlane = in_digihit->view + 3*(gLayer - 1);
	
	if( (gPlane <= 0) || (static_cast<unsigned int>(gPlane) > FDC_NUM_PLANES)) {
		sprintf(str, "Bad gPlane # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", gPlane, FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	// strip and wire planes have different numbers of elements
	if( (in_digihit->strip <= 0) || (static_cast<unsigned int>(in_digihit->strip) > STRIPS_PER_PLANE)) {
	    sprintf(str, "Bad strip # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", in_digihit->strip, STRIPS_PER_PLANE);
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[gPlane-1][in_digihit->strip-1];
}

const double DFDCHit_factory::GetConstant(const fdc_digi_constants_t &the_table,
					  const DFDCWireDigiHit *in_digihit) const {

	char str[256];

	int gLayer = in_digihit->chamber + 6*(in_digihit->package - 1);
	int gPlane = 2 + 3*(gLayer - 1);
	
	if( (gPlane <= 0) || (static_cast<unsigned int>(gPlane) > FDC_NUM_PLANES)) {
		sprintf(str, "Bad gPlane # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", gPlane, FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	// strip and wire planes have different numbers of elements
	if( (in_digihit->wire <= 0) || (static_cast<unsigned int>(in_digihit->wire) > WIRES_PER_PLANE)) {
	    sprintf(str, "Bad wire # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", in_digihit->wire, WIRES_PER_PLANE);
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[gPlane-1][in_digihit->wire-1];
}

const double DFDCHit_factory::GetConstant(const fdc_digi_constants_t &the_table,
					  const DFDCHit *in_hit) const {

	char str[256];
	
	if( (in_hit->gPlane <= 0) || (static_cast<unsigned int>(in_hit->gPlane) > FDC_NUM_PLANES)) {
		sprintf(str, "Bad gPlane # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", in_hit->gPlane, FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	// strip and wire planes have different numbers of elements
	if( (in_hit->element <= 0) || (static_cast<unsigned int>(in_hit->element) > the_table[in_hit->gPlane].size())) {
	        sprintf(str, "Bad element # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %zu", in_hit->element, the_table[in_hit->gPlane].size());
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[in_hit->gPlane-1][in_hit->element-1];
}
/*
const double DFDCHit_factory::GetConstant(const fdc_digi_constants_t &the_table,
					  const DTranslationTable *ttab,
					  const int in_rocid, const int in_slot, const int in_channel) const {

	char str[256];
	
	DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
	DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);
	
	if( channel_info.det_sys == DTranslationTable::FDC_CATHODES ) {  
	    // FDC Cathodes
	    int gLayer = channel_info.fdc_cathodes.chamber + 6*(channel_info.fdc_cathodes.package - 1);
	    int gPlane = channel_info.fdc_cathodes.view + 3*(gLayer - 1);

	    if( (gPlane <= 0) || (gPlane > FDC_NUM_PLANES)) {
		sprintf(str, "Bad gPlane # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", gPlane, FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	    }
	    // strip and wire planes have different numbers of elements
	    if( (channel_info.fdc_cathodes.strip <= 0) 
		|| (channel_info.fdc_cathodes.strip > STRIPS_PER_PLANE)) {
		sprintf(str, "Bad strip # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", channel_info.fdc_cathodes.strip, STRIPS_PER_PLANE);
		cerr << str << endl;
		throw JException(str);
	    }

	    return the_table[gPlane-1][channel_info.fdc_cathodes.strip-1];
	} else if( channel_info.det_sys == DTranslationTable::FDC_WIRES ) {  
	    // FDC Wirees
	    int gLayer = channel_info.fdc_wires.chamber + 6*(channel_info.fdc_wires.package - 1);
	    int gPlane = 2 + 3*(gLayer - 1);  // wire planes are always layer 2

	    if( (gPlane <= 0) || (gPlane > FDC_NUM_PLANES)) {
		sprintf(str, "Bad gPlane # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", gPlane, FDC_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	    }
	    // strip and wire planes have different numbers of elements
	    if( (channel_info.fdc_wires.wire <= 0) 
		|| (channel_info.fdc_wires.wire > WIRES_PER_PLANE)) {
		sprintf(str, "Bad strip # requested in DFDCHit_factory::GetConstant()! requested=%d , should be %ud", channel_info.fdc_wires.wire, WIRES_PER_PLANE);
		cerr << str << endl;
		throw JException(str);
	    }

	    return the_table[gPlane-1][channel_info.fdc_wires.wire-1];
	} else {
	    sprintf(str, "Got bad detector type in DFDCHit_factory::GetConstant()! requested=%d", channel_info.module_type);
	    cerr << str << endl;
	    throw JException(str);

	    return -1.;  // should never reach here!
	}	
}
*/
