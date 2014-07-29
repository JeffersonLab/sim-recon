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
#include "DBCALGeometry.h"
#include "DBCALHit_factory.h"
#include <DAQ/Df250PulseIntegral.h>
using namespace jana;

#define BCAL_NUM_MODULES   48
#define BCAL_NUM_LAYERS     4
#define BCAL_NUM_SECTORS    4

#define kMaxChannels     1536

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
	/// set the base conversion scales
	a_scale    = 0.1;   // to get units of MeV
	//  Crude calibration 
	//    A minimally ionising particle deposits and integral of 230 ADC counts per cell, 
	//    which corresponds to approximately 22 MeV.  Thus, the factor is 0.1 to get MeV
	//a_pedestal = 10000;  // default pedestal of 100 ADC units over 100 samples 
	t_scale    = 0.0625;   // 62.5 ps/count

	/// Read in calibration constants
	vector<double> raw_gains;
	vector<double> raw_pedestals;
	vector<double> raw_time_offsets;

	jout << "In DBCALHit_factory, loading constants..." << endl;

	if(eventLoop->GetCalib("/BCAL/ADC_gains", raw_gains))
	    jout << "Error loading /BCAL/ADC_gains !" << endl;
	if(eventLoop->GetCalib("/BCAL/ADC_pedestals", raw_pedestals))
	    jout << "Error loading /BCAL/ADC_pedestals !" << endl;
	if(eventLoop->GetCalib("/BCAL/ADC_timing_offsets", raw_time_offsets))
	    jout << "Error loading /BCAL/ADC_timing_offsets !" << endl;

	FillCalibTable(gains, raw_gains);
	FillCalibTable(pedestals, raw_pedestals);
	FillCalibTable(time_offsets, raw_time_offsets);


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
		
		// Get pedestal.  Prefer associated event pedestal if it exist.
		// Otherwise, use the average pedestal from CCDB
		double pedestal = GetConstant(pedestals,digihit);
		vector<const Df250PulseIntegral*> PIvect;
		digihit->Get(PIvect);
		if(!PIvect.empty()){
		    const Df250PulseIntegral *PIobj = PIvect[0];
		    pedestal = PIobj->pedestal;
		}

		DBCALHit *hit = new DBCALHit;
		hit->module = digihit->module;
		hit->layer  = digihit->layer;
		hit->sector = digihit->sector;
		hit->end    = digihit->end;
		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;

		double gain = GetConstant(gains,digihit);

		hit->E = a_scale * gain * (A - pedestal);
		hit->t = t_scale * (T - GetConstant(time_offsets,digihit));
		
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

//------------------
// GetConstant
//------------------
double DBCALHit_factory::GetConstant( map<int,cell_calib_t> &the_table, 
				     const DBCALDigiHit *the_digihit) 
{
    int the_cell = DBCALGeometry::cellId(the_digihit->module,
					 the_digihit->layer,
					 the_digihit->sector);

    if(the_digihit->end == DBCALGeometry::kUpstream) {
	// handle the upstream end
	return the_table[the_cell].first;
    } else {
	// handle the downstream end
	return the_table[the_cell].second;
    }

}


//------------------
// FillCalibTable
//------------------
void DBCALHit_factory::FillCalibTable( map<int,cell_calib_t> &table, 
				       const vector<double> &raw_table) 
{
    
    int channel = 0;
    
    table.clear();

    for(int module=1; module<=BCAL_NUM_MODULES; module++) {
	for(int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
	    for(int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
		if( (channel > kMaxChannels) || (channel+1 > kMaxChannels) )   // sanity check
		    return;
		
		int cell_id = DBCALGeometry::cellId(module,layer,sector);
		
		table[cell_id] = cell_calib_t(raw_table[channel],raw_table[channel+1]);

		channel += 2;
		
	    }
	}
    }

}
