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
#include <DAQ/Df125PulseIntegral.h>
using namespace jana;

#define CDC_MAX_CHANNELS  3522

static double DIGI_THRESHOLD = -1000000.0;

//------------------
// init
//------------------
jerror_t DCDCHit_factory::init(void)
{
        gPARMS->SetDefaultParameter("CDC:DIGI_THRESHOLD",DIGI_THRESHOLD, "Do not convert CDC digitized hits into DCDCHit objects that would have q less than this");

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCDCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
        // calculate the number of straws in each ring
 	CalcNstraws(eventLoop, runnumber, Nstraws);
	Nrings = Nstraws.size();

	/// set the base conversion scales
	//a_scale    = 1.0E6/1.3E5; 
	a_scale    = 4.0E3/1.0E2; 
	t_scale    = 8.0/10.0;    // 8 ns/count and integer time is in 1/10th of sample

	/// Read in calibration constants
        vector<double> raw_gains;
        vector<double> raw_pedestals;
        vector<double> raw_time_offsets;

	jout << "In DFCALHit_factory, loading constants..." << endl;

        if(eventLoop->GetCalib("/CDC/wire_gains", raw_gains))
	    jout << "Error loading /CDC/wire_gains !" << endl;
        if(eventLoop->GetCalib("/CDC/pedestals", raw_pedestals))
	    jout << "Error loading /CDC/pedestals !" << endl;
	if(eventLoop->GetCalib("/CDC/timing_offsets", raw_time_offsets))
	    jout << "Error loading /CDC/timing_offsets !" << endl;

	// fill the tables
        FillCalibTable(gains, raw_gains, Nstraws);
        FillCalibTable(pedestals, raw_pedestals, Nstraws);
        FillCalibTable(time_offsets, raw_time_offsets, Nstraws);
	
	// Verify that the right number of rings was read for each set of constants
	char str[256];
	if(gains.size() != Nrings){
		sprintf(str, "Bad # of rings for CDC gain from CCDB! CCDB=%zu , should be %d", gains.size(), Nrings);
		cerr << str << endl;
		throw JException(str);
	}
	if(pedestals.size() != Nrings){
		sprintf(str, "Bad # of rings for CDC pedestal from CCDB! CCDB=%zu , should be %d", pedestals.size(), Nrings);
		cerr << str << endl;
		throw JException(str);
	}
	if(time_offsets.size() != Nrings){
		sprintf(str, "Bad # of rings for CDC time_offset from CCDB! CCDB=%zu , should be %d", time_offsets.size(), Nrings);
		cerr << str << endl;
		throw JException(str);
	}
	
	// Verify the right number of straws was read for each ring for each set of constants
	for(unsigned int i=0; i<Nrings; i++){
		if(gains[i].size() != Nstraws[i]){
			sprintf(str, "Bad # of straws for CDC gain from CCDB! CCDB=%zu , should be %d for ring %d", gains[i].size(), Nstraws[i], i+1);
			cerr << str << endl;
			throw JException(str);
		}
		if(pedestals[i].size() != Nstraws[i]){
			sprintf(str, "Bad # of straws for CDC pedestal from CCDB! CCDB=%zu , should be %d for ring %d", pedestals[i].size(), Nstraws[i], i+1);
			cerr << str << endl;
			throw JException(str);
		}
		if(time_offsets[i].size() != Nstraws[i]){
			sprintf(str, "Bad # of straws for CDC time_offset from CCDB! CCDB=%zu , should be %d for ring %d", time_offsets[i].size(), Nstraws[i], i+1);
			cerr << str << endl;
			throw JException(str);
		}
	}

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
	char str[256];
	for(unsigned int i=0; i<digihits.size(); i++){
		const DCDCDigiHit *digihit = digihits[i];
		const int &ring  = digihit->ring;
		const int &straw = digihit->straw;
		
		// Make sure ring and straw are in valid range
		if(ring > (int)Nrings){
			sprintf(str, "DCDCDigiHit ring out of range! ring=%d (should be 1-%d)", ring, Nrings);
			throw JException(str);
		}
		if(straw > (int)Nstraws[ring-1]){
			sprintf(str, "DCDCDigiHit straw out of range! straw=%d for ring=%d (should be 1-%d)", straw, ring, Nstraws[ring-1]);
			throw JException(str);
		}
		
		// Get pedestal. Preference is given to pedestal measured
		// for event. Otherwise, use statistical one from CCDB
		double pedestal = pedestals[ring-1][straw-1];
		vector<const Df125PulseIntegral*> PIvect;
		digihit->Get(PIvect);
		if(!PIvect.empty()) pedestal = (double)PIvect[0]->pedestal;
 		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		
		double q = a_scale * gains[ring-1][straw-1] * (A - pedestal);
		double t = t_scale * (T - time_offsets[ring-1][straw-1]);

		if( q < DIGI_THRESHOLD) continue;

		DCDCHit *hit = new DCDCHit;
		hit->ring  = ring;
		hit->straw = straw;

		// note that ring/straw counting starts at 1
		hit->q = q;
		hit->t = t;
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
void DCDCHit_factory::CalcNstraws(jana::JEventLoop *eventLoop, int runnumber, vector<unsigned int> &Nstraws)
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
				     vector<unsigned int> &Nstraws)
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
	if(straw == (int)Nstraws[ring]) {
	    ring++;
	    straw = 0;
	}

        table[ring].push_back( raw_table[channel] );
    }

}
