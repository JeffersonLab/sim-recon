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
#include <DAQ/Df250PulseIntegral.h>
using namespace jana;


#define FCAL_MAX_CHANNELS   2800

//------------------
// init
//------------------
jerror_t DFCALHit_factory::init(void)
{
        // initialize calibration tables
        vector< vector<double > > new_gains(kBlocksTall, vector<double>(kBlocksWide));
        vector< vector<double > > new_pedestals(kBlocksTall, vector<double>(kBlocksWide));
	vector< vector<double > > new_t0s(kBlocksTall, vector<double>(kBlocksWide));
	vector< vector<double > > new_qualities(kBlocksTall, vector<double>(kBlocksWide));

	gains = new_gains;
	pedestals = new_pedestals;
	time_offsets = new_t0s;
	block_qualities = new_qualities;
    
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DFCALHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	// extract the FCAL Geometry
	vector<const DFCALGeometry*> fcalGeomVect;
	eventLoop->Get( fcalGeomVect );
	if(fcalGeomVect.size()<1)return OBJECT_NOT_AVAILABLE;
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);
	
	/// set the base conversion scales
	a_scale    = 4.0E1/2.5E5; 
	t_scale    = 0.0625;   // 62.5 ps/count

	/// Read in calibration constants
	vector< double > raw_gains;
	vector< double > raw_pedestals;
	vector< double > raw_time_offsets;
	vector< double > raw_block_qualities;    // we should change this to an int?

	jout << "In DFCALHit_factory, loading constants..." << endl;

	if(eventLoop->GetCalib("/FCAL/gains", raw_gains))
	    jout << "Error loading /FCAL/gains !" << endl;
	if(eventLoop->GetCalib("/FCAL/pedestals", raw_pedestals))
	    jout << "Error loading /FCAL/pedestals !" << endl;
	if(eventLoop->GetCalib("/FCAL/timing_offsets", raw_time_offsets))
	    jout << "Error loading /FCAL/timing_offsets !" << endl;
	if(eventLoop->GetCalib("/FCAL/block_quality", raw_block_qualities))
	    jout << "Error loading /FCAL/block_quality !" << endl;

	FillCalibTable(gains, raw_gains, fcalGeom);
	FillCalibTable(pedestals, raw_pedestals, fcalGeom);
	FillCalibTable(time_offsets, raw_time_offsets, fcalGeom);
	FillCalibTable(block_qualities, raw_block_qualities, fcalGeom);

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

		// Get pedestal.  Prefer associated event pedestal if it exist.
		// Otherwise, use the average pedestal from CCDB
		double pedestal = pedestals[digihit->row][digihit->column];
		vector<const Df250PulseIntegral*> PIvect;
		digihit->Get(PIvect);
		if(!PIvect.empty()){
		    const Df250PulseIntegral *PIobj = PIvect[0];
		    pedestal = PIobj->pedestal;
		}
		
		DFCALHit *hit = new DFCALHit;
		hit->row    = digihit->row;
		hit->column = digihit->column;

		// throw away hits from bad or noisy channels
		fcal_quality_state quality = static_cast<fcal_quality_state>(block_qualities[hit->row][hit->column]);
		if( (quality==BAD) || (quality==NOISY) ) 
		    return NOERROR;

		// Apply calibration constants
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		hit->E = a_scale * gains[hit->row][hit->column] * (A - pedestal);
		hit->t = t_scale * (T - time_offsets[hit->row][hit->column]);

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

//------------------
// FillCalibTable
//------------------
void DFCALHit_factory::FillCalibTable( fcal_digi_constants_t &table, 
				       const vector<double> &raw_table, 
				       const DFCALGeometry &fcalGeom)
{
    for(int channel=0; channel<static_cast<int>(raw_table.size()); channel++) {
	// make sure that we don't try to load info for channels that don't exist
	if(channel == FCAL_MAX_CHANNELS) break;
	
	int row = fcalGeom.row(channel);
	int col = fcalGeom.column(channel);

	table[row][col] = raw_table[channel];
    }
    
}
