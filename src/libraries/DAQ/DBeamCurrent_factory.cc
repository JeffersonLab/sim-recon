// $Id$
//
//    File: DBeamCurrent_factory.cc
// Created: Tue Feb 21 04:25:04 EST 2017
// Creator: davidl (on Linux gluon48.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
using namespace std;

#include <DAQ/DCODAEventInfo.h>
#include "DBeamCurrent_factory.h"
using namespace jana;



//------------------
// init
//------------------
jerror_t DBeamCurrent_factory::init(void)
{
	BEAM_ON_MIN_nA  = 10.0;  // nA
	BEAM_TRIP_MIN_T = 3.0;   // seconds
	
	gPARMS->SetDefaultParameter("BEAM_ON_MIN_nA", BEAM_ON_MIN_nA, "Minimum current in nA to consider the beam \"on\" by DBeamCurrent");
	gPARMS->SetDefaultParameter("BEAM_TRIP_MIN_T", BEAM_TRIP_MIN_T, "Minimum amount of time in seconds that event is away from beam trips to be considered fiducial");

	ticks_per_sec      = 250.011E6; // 250MHz clock (may be overwritten with calib constant in brun)
	rcdb_start_time    = 0;       // unix time of when 250MHz clock was reset. (overwritten below)
	rcdb_250MHz_offset_tics = 0;  // offset between 250MHz clock zero and RCDB recorded start time of event (overwritten below)

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBeamCurrent_factory::brun(jana::JEventLoop *loop, int32_t runnumber)
{
	// Clear maps in case we are called more than once
	boundaries.clear();
	trip.clear();
	recover.clear();

	// Data is stored in CCDB as a large ASCII string with
	// pairs of values:
	//
	//  t  Ibeam
	//
	//  where
	//      t = time in seconds since start of run
	//  Ibeam = beam current in nA
	//
	
	// n.b. we have to do this by first getting the JCalibration
	// object and then using it directly. If we use the GetCalib()
	// method in JEventLoop, it will try parsing the string and 
	// return more than a 1 element map.
	map<string,string> mstr;
	loop->GetJCalibration()->GetCalib("/ELECTRON_BEAM/current_map_epics", mstr);
	if(mstr.empty()) return NOERROR;
	string &electron_beam_current = mstr.begin()->second;
	
	map<string,string> mcalib;
	loop->GetJCalibration()->GetCalib("/ELECTRON_BEAM/timestamp_to_unix", mcalib);
	if(mcalib.size() == 3){
		//ticks_per_sec           = atof(mcalib["tics_per_sec"].c_str());
		rcdb_250MHz_offset_tics = atoi(mcalib["rcdb_250MHz_offset_tics"].c_str());
		rcdb_start_time         = atoi(mcalib["rcdb_start_time"].c_str());
	}

	
	// Parse text to create Boundary objects and maps of trip/recovery points
	istringstream ss(electron_beam_current);
	double last_Ibeam = 0.0;
	for(string line; getline(ss, line, '\n'); ){
		double t     = atof(line.c_str());
		double Ibeam = atof(line.substr(1+line.find(" ")).c_str());
		Boundary b = {t, Ibeam, 0.0, 0.0};
		boundaries.push_back(b);
		
		// Record trip/recovery points
		if( fabs(Ibeam - last_Ibeam) > BEAM_ON_MIN_nA ){
			if(Ibeam < BEAM_ON_MIN_nA){
				trip.push_back(t);
			}else if(last_Ibeam < BEAM_ON_MIN_nA){
				recover.push_back(t);
			}
		}
		last_Ibeam = Ibeam;
	}
	
	// Loop through all boundaries and update the time to
	// all trip points.
	double t_max = IntegratedTime();
	for(auto &b : boundaries){
	
		// Prev and next values should correspond appropriately
		// to trip or recovery boundaries as to whether this 
		// boundary starts a "ON" or "OFF" period. i.e.
		//
		// If beam current is "OFF" then prev should be the
		// previous trip time and next the next recovery time.
		//
		// If beam current is "ON" then prev should be the
		// previous recovery time and next the next trip time.

		double t_prev = 0.0;
		double t_next = t_max;
		if(b.Ibeam < BEAM_ON_MIN_nA){
			// beam "OFF"
			for( double t : recover ) if( t>b.t ) {t_next = t; break;}
			for( double t : trip    ) if( t>b.t ) {break;} else {t_prev = t;}
		}else{
			// beam "ON"
			for( double t : trip    ) if( t>b.t ) {t_next = t; break;}
			for( double t : recover ) if( t>b.t ) {break;} else {t_prev = t;}
		}

		b.t_trip_prev = b.t - t_prev;
		b.t_trip_next = t_next - b.t;
	}
	
	jout << "Electron beam current trip map loaded with " << boundaries.size() << " boundaries (" << trip.size() << " trips over " << t_max << " sec)" << endl;

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBeamCurrent_factory::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	if(boundaries.empty()) return NOERROR;

	// Get time of this event relative to start of run in seconds
	// n.b. don't use loop->GetSingle() here. It results in infinite
	// recursion.
	vector<const DCODAEventInfo*> codainfos;
	loop->Get(codainfos);
	if(codainfos.empty()) return NOERROR;
	const DCODAEventInfo *codainfo = codainfos[0];

	// Get tme relative to RCDB recorded start time of event
	// (all times in trip map are recorded relative to this as well)	
	double t = (codainfo->avg_timestamp - (double)rcdb_250MHz_offset_tics)/ticks_per_sec;

	// Find closest entry given current time
	auto it = boundaries.begin();
	while(it!=boundaries.end()){
		if(it->t > t) break;
		it++;
	}
	
	if(it != boundaries.begin() ){
		it--;
		Boundary &b = *it;
		double t_rel = t - b.t;	// time relative to previous boundary

		DBeamCurrent *bc = new DBeamCurrent;
		bc->Ibeam  = b.Ibeam;
		bc->t_prev = b.t_trip_prev + t_rel;
		bc->t_next = b.t_trip_next - t_rel;
		bc->is_fiducial = false;
		if(b.Ibeam >= BEAM_ON_MIN_nA){
			if(bc->t_prev>=BEAM_TRIP_MIN_T){
				if(bc->t_next>=BEAM_TRIP_MIN_T) bc->is_fiducial=true;
			}
		}
		
		_data.push_back(bc);
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBeamCurrent_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBeamCurrent_factory::fini(void)
{
	return NOERROR;
}

//------------------
// IntegratedFiducialTime
//------------------
double DBeamCurrent_factory::IntegratedFiducialTime(double t_start, double t_end)
{
	/// Loop over all boundaries to find the total fiducial time
	/// between the two given times. Both t_start and t_end should
	/// be in seconds and relative to the start of the run (NOT
	/// the start of this particular file!) The value returned is
	/// in seconds. The fraction of the time the beam was on for
	/// the entire run can be obtained with:
	///
	///    IntegratedFiducialTime()/IntegratedTime()
	///
	
	
	// Loop over start of "recover" regions
	double t_fiducial = 0.0; // total fiducial time so far
	double t = t_start;      // point we have already integrated to
	for(double t1 : recover){
	
		// check if next recovery region starts after where we are
		if(t1 > t)t = t1;
	
		// find start of next "tripped" region
		double t2 = t_end;
		for(double tt : trip)if(tt>t || tt>t_end){t2 = tt; break;}
		
		// At this point t is between t1 and t2 (possibly at t1)
		// which is a beam "ON" region. Calculate time in this
		// region excluding BEAM_TRIP_MIN_T seconds from edges.
		
		// move t to start of valid integration range
		if( (t-t1)<BEAM_TRIP_MIN_T ) t = t1 + BEAM_TRIP_MIN_T;
		
		// make t3 the end of integration range being careful
		// not to include time past t_end or within BEAM_TRIP_MIN_T
		// of next trip.
		double t3 = t2 - BEAM_TRIP_MIN_T;
		if(t3>t_end) t3 = t_end;
		
		// Calculate how much fiducial time in this region is
		// within our integration window and if it is positive,
		// then add it.
		double delta_t = t3-t;
		if(delta_t > 0.0) t_fiducial += delta_t;
		
		// Move integration point
		t = t2;
		
		// break early if we have included whole window
		if(t >= t_end) break;
	}

	return t_fiducial;
}

//------------------
// IntegratedTime
//------------------
double DBeamCurrent_factory::IntegratedTime(void)
{
	/// Return total integrated time for this run.
	/// WARNING: This is NOT the time of only the
	/// current file. See IntegratedFiducialTime
	/// for more details.
	if(boundaries.empty()) return 0.0;
	return boundaries[boundaries.size()-1].t;
}


