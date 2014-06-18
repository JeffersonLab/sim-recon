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
#include "DCDCWire.h"
using namespace jana;

#define CDC_MAX_CHANNELS  3522

static int USE_MC_CALIB = 0;

//------------------
// init
//------------------
jerror_t DCDCHit_factory::init(void)
{
        // should we use calibrations for simulated data? - this is a temporary workaround
        gPARMS->SetDefaultParameter("DIGI:USEMC",USE_MC_CALIB);

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
        // calculate the number of straws in each ring
        vector<int> Nstraws;
	CalcNstraws(eventLoop, runnumber, Nstraws);

	/// set the base conversion scales
	a_scale    = 1.0E6/1.3E5; 
	t_scale    = 8.0;    // 8 ns/count

	/// Read in calibration constants
        vector<double> raw_gains;
        vector<double> raw_pedestals;
        vector<double> raw_time_offsets;

	jout << "In DFCALHit_factory, loading constants..." << endl;

        if(eventLoop->GetCalib("/CDC/wire_gains", raw_gains))
	    jout << "Error loading /CDC/wire_gains !" << endl;
        if(eventLoop->GetCalib("/CDC/pedestals", raw_pedestals))
	    jout << "Error loading /CDC/pedestals !" << endl;
	if(USE_MC_CALIB>0) {
	    if(eventLoop->GetCalib("/CDC/timing_offsets::mc", raw_time_offsets))
		jout << "Error loading /CDC/timing_offsets !" << endl;
	} else {
	    if(eventLoop->GetCalib("/CDC/timing_offsets", raw_time_offsets))
		jout << "Error loading /CDC/timing_offsets !" << endl;
	}

	// fill the tables
        FillCalibTable(gains, raw_gains, Nstraws);
        FillCalibTable(pedestals, raw_pedestals, Nstraws);
        FillCalibTable(time_offsets, raw_time_offsets, Nstraws);

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

		// note that ring/straw counting starts at 1
		hit->q = a_scale * gains[hit->ring-1][hit->straw-1] * (A - pedestals[hit->ring-1][hit->straw-1]);
		hit->t = t_scale * (T - time_offsets[hit->ring-1][hit->straw-1]);
		hit->d = 0.0;
		hit->itrack = -1;
		hit->ptype = 0;

		//cerr << "CDC Hit:  ring = " << hit->ring << "  straw = " << hit->straw << endl;
		
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

//------------------
// CalcNstraws
//------------------
void DCDCHit_factory::CalcNstraws(jana::JEventLoop *eventLoop, int runnumber, vector<int> &Nstraws)
{
    DGeometry *dgeom;
    vector<vector<DCDCWire *> >cdcwires;

    // Get pointer to DGeometry object
    DApplication* dapp=dynamic_cast<DApplication*>(eventLoop->GetJApplication());
    dgeom  = dapp->GetDGeometry(runnumber);
  
    // Get the CDC wire table from the XML
    dgeom->GetCDCWires(cdcwires);
  
    // Fill array with the number of straws for each layer
    for (unsigned int i=0;i<cdcwires.size();i++){
	Nstraws.push_back( cdcwires[i].size() );
    }

    // clear up all of the wire information
    for (unsigned int i=0;i<cdcwires.size();i++) {
	for (unsigned int j=0;j<cdcwires[i].size();j++) {
	    delete cdcwires[i][j];
	}
    }    
    cdcwires.clear();
}


//------------------
// FillCalibTable
//------------------
void DCDCHit_factory::FillCalibTable(vector< vector<double> > &table, vector<double> &raw_table, 
				     vector<int> &Nstraws)
{
    int ring = 0;
    int straw = 0;

    // reset table before filling it
    table.clear();
    table.resize( Nstraws.size() );

    for(int channel=0; channel<static_cast<int>(raw_table.size()); channel++,straw++) {
        // make sure that we don't try to load info for channels that don't exist
        if(channel == CDC_MAX_CHANNELS) break;

	// if we've hit the end of the ring, move on to the next
	if(straw == Nstraws[ring]) {
	    ring++;
	    straw = 0;
	}

        table[ring].push_back( raw_table[channel] );
    }
/*
    cerr << "LOADED TABLE " << endl;
    cerr << " number of rings = " << table.size() << endl;
    for(int i=0; i<table.size(); i++)
	cerr << " ring #" << i+1 << " has " << table[i].size() << " straws" << endl;
*/
}
