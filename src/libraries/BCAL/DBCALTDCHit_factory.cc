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

#define BCAL_NUM_MODULES   48
#define BCAL_NUM_TDC_LAYERS     3
#define BCAL_NUM_SECTORS    4

#define kMaxChannels     1152

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
	/// set the base conversion scale
	t_scale    = 0.060;    // 60 ps/count
	//t_offset   = 0;

	/// Read in calibration constants
	vector<double> raw_time_offsets;

	jout << "In DBCALTDCHit_factory, loading constants..." << endl;

	if(eventLoop->GetCalib("/BCAL/TDC_offsets", raw_time_offsets))
	    jout << "Error loading /BCAL/TDC_offsets !" << endl;

	FillCalibTable(time_offsets, raw_time_offsets);


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
		hit->t = t_scale * (T - GetConstant(time_offsets,digihit));
	
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


//------------------
// GetConstant
//------------------
double DBCALTDCHit_factory::GetConstant( map<int,cell_calib_t> &the_table, 
					 const DBCALTDCDigiHit *the_digihit) 
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
void DBCALTDCHit_factory::FillCalibTable( map<int,cell_calib_t> &table, 
					  const vector<double> &raw_table) 
{
    
    int channel = 0;

    table.clear();

    for(int module=1; module<=BCAL_NUM_MODULES; module++) {
	for(int layer=1; layer<=BCAL_NUM_TDC_LAYERS; layer++) {
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
