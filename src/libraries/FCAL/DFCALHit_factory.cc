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


//#define FCAL_MAX_CHANNELS   2800

//------------------
// init
//------------------
jerror_t DFCALHit_factory::init(void)
{
        // initialize calibration tables
	vector< vector<double > > new_gains(kBlocksTall, 
					    vector<double>(kBlocksWide));
	vector< vector<double > > new_pedestals(kBlocksTall, 
						vector<double>(kBlocksWide));
	vector< vector<double > > new_t0s(kBlocksTall, 
					  vector<double>(kBlocksWide));
	vector< vector<double > > new_qualities(kBlocksTall, 
						vector<double>(kBlocksWide));

	gains = new_gains;
	pedestals = new_pedestals;
	time_offsets = new_t0s;
	block_qualities = new_qualities;

	// default values
	a_scale = 0.;
	t_scale = 0.;
    
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
	char str[256];

	// extract the FCAL Geometry (for positionOnFace())
	vector<const DFCALGeometry*> fcalGeomVect;
	eventLoop->Get( fcalGeomVect );
	if(fcalGeomVect.size()<1)return OBJECT_NOT_AVAILABLE;
	const DFCALGeometry& fcalGeom = *(fcalGeomVect[0]);
	
	vector<const DFCALDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DFCALDigiHit *digihit = digihits[i];

		// Check to see if the hit corresponds to a valid channel
		if(fcalGeom.isBlockActive(digihit->row,digihit->column) == false) {
		    sprintf(str, "DFCALHit corresponds to inactive channel!  row=%d, col=%d", 
			    digihit->row, digihit->column);
		    throw JException(str);
		}

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
    char str[256];

    // sanity check that we have the right geometry
    // (deprecate this?) 
    if(fcalGeom.numActiveBlocks() != FCAL_MAX_CHANNELS) {
	sprintf(str, "FCAL geometry is wrong size! channels=%d (should be %d)", 
		fcalGeom.numActiveBlocks(), FCAL_MAX_CHANNELS);
	throw JException(str);
    }

    // check to see if the table is the right size
    if( fcalGeom.numActiveBlocks() != static_cast<int>(raw_table.size()) ) {
         sprintf(str, "FCAL constant table is wrong size! channels=%d (should be %d)", 
		 fcalGeom.numActiveBlocks(), static_cast<int>(raw_table.size()));
         throw JException(str);
    }

    for(int channel=0; channel < static_cast<int>(raw_table.size()); channel++) {
	// make sure that we don't try to load info for channels that don't exist
	if(channel == fcalGeom.numActiveBlocks()) break;
	
	int row = fcalGeom.row(channel);
	int col = fcalGeom.column(channel);

	// results from DFCALGeometry should be self consistent, but add in some
	// sanity checking just to be sure
	if(fcalGeom.isBlockActive(row,col) == false) {
	    sprintf(str, "Loading FCAL constant for inactive channel!  row=%d, col=%d", row, col);
	    throw JException(str);
	}

	table[row][col] = raw_table[channel];
    }
}

//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
					   const int in_row, const int in_column) const {
	
	char str[256];
	
	if( (in_row <= 0) || (in_row > kBlocksTall)) {
		sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", in_row, kBlocksTall);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_column <= 0) || (in_column > kBlocksWide)) {
		sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", in_column, kBlocksWide);
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[in_row][in_column];
}

const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
					   const DFCALDigiHit *in_digihit) const {

	char str[256];
	
	if( (in_digihit->row <= 0) || (in_digihit->row > kBlocksTall)) {
		sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", in_digihit->row, kBlocksTall);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->column <= 0) || (in_digihit->column > kBlocksWide)) {
		sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", in_digihit->column, kBlocksWide);
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[in_digihit->row][in_digihit->column];
}

const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
					   const DFCALHit *in_hit) const {

	char str[256];
	
	if( (in_hit->row <= 0) || (in_hit->row > kBlocksTall)) {
		sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", in_hit->row, kBlocksTall);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_hit->column <= 0) || (in_hit->column > kBlocksWide)) {
		sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", in_hit->column, kBlocksWide);
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[in_hit->row][in_hit->column];
}
/*
const double DFCALHit_factory::GetConstant(const fcal_digi_constants_t &the_table,
					   const DTranslationTable *ttab,
					   const int in_rocid, const int in_slot, const int in_channel) const {

	char str[256];
	
	DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
	DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);
	
	if( (channel_info.fcal.row <= 0) 
	    || (channel_info.fcal.row > static_cast<unsigned int>(kBlocksTall))) {
		sprintf(str, "Bad row # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", channel_info.fcal.row, kBlocksTall);
		cerr << str << endl;
		throw JException(str);
	}
	if( (channel_info.fcal.col <= 0) 
	    || (channel_info.fcal.col > static_cast<unsigned int>(kBlocksWide))) {
	    sprintf(str, "Bad column # requested in DFCALHit_factory::GetConstant()! requested=%d , should be %ud", channel_info.fcal.row, kBlocksWide);
		cerr << str << endl;
		throw JException(str);
	}
	
	return the_table[channel_info.fcal.row][channel_info.fcal.col];
}
*/
