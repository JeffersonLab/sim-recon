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
// FillCalibTable
//------------------
void DBCALTDCHit_factory::FillCalibTable( map<int,cell_calib_t> &table, 
				       const vector<double> &raw_table) 
{
    char str[256];
    int channel = 0;

    // reset the table before filling it
    table.clear();

    for(int module=1; module<=BCAL_NUM_MODULES; module++) {
	for(int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
	    for(int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
		if( (channel > BCAL_MAX_CHANNELS) || (channel+1 > BCAL_MAX_CHANNELS) ) {  // sanity check
		    sprintf(str, "Too many channels for BCAL table! channel=%d (should be %d)", 
			    channel, BCAL_MAX_CHANNELS);
		    cerr << str << endl;
		    throw JException(str);
		}
		
		int cell_id = DBCALGeometry::cellId(module,layer,sector);
		
		table[cell_id] = cell_calib_t(raw_table[channel],raw_table[channel+1]);

		channel += 2;
	    }
	}
    }

    // check to make sure that we loaded enough channels
    if(channel != BCAL_MAX_CHANNELS) { 
	sprintf(str, "Not enough channels for BCAL table! channel=%d (should be %d)", 
		channel, BCAL_MAX_CHANNELS);
	cerr << str << endl;
	throw JException(str);
    }
}

//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DBCALTDCHit_factory::GetConstant( const bcal_digi_constants_t &the_table, 
					    const int in_module, const int in_layer,
					    const int in_sector, const int in_end) const
{
    	char str[256];
	
	if( (in_module <= 0) || (in_module > BCAL_NUM_MODULES)) {
		sprintf(str, "Bad module # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_module, BCAL_NUM_MODULES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_layer <= 0) || (in_layer > BCAL_NUM_LAYERS)) {
		sprintf(str, "Bad layer # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_layer, BCAL_NUM_LAYERS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_sector <= 0) || (in_sector > BCAL_NUM_SECTORS)) {
		sprintf(str, "Bad sector # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_sector, BCAL_NUM_SECTORS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_end != DBCALGeometry::kUpstream) && (in_end != DBCALGeometry::kDownstream) ) {
		sprintf(str, "Bad end # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 0-1", in_end);
		cerr << str << endl;
		throw JException(str);
	}

	const int the_cell = DBCALGeometry::cellId( in_module,
						    in_layer,
						    in_sector);
	
	if(in_end == DBCALGeometry::kUpstream) {
	    // handle the upstream end
	    return the_table.at(the_cell).first;
	} else {
	    // handle the downstream end
	    return the_table.at(the_cell).second;
	}

}

const double DBCALTDCHit_factory::GetConstant( const bcal_digi_constants_t &the_table, 
					    const DBCALTDCHit *in_hit) const
{
    	char str[256];
	
	if( (in_hit->module <= 0) || (in_hit->module > BCAL_NUM_MODULES)) {
		sprintf(str, "Bad module # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_hit->module, BCAL_NUM_MODULES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_hit->layer <= 0) || (in_hit->layer > BCAL_NUM_LAYERS)) {
		sprintf(str, "Bad layer # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_hit->layer, BCAL_NUM_LAYERS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_hit->sector <= 0) || (in_hit->sector > BCAL_NUM_SECTORS)) {
		sprintf(str, "Bad sector # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_hit->sector, BCAL_NUM_SECTORS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_hit->end != DBCALGeometry::kUpstream) && (in_hit->end != DBCALGeometry::kDownstream) ) {
		sprintf(str, "Bad end # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 0-1", in_hit->end);
		cerr << str << endl;
		throw JException(str);
	}

	const int the_cell = DBCALGeometry::cellId(in_hit->module,
						   in_hit->layer,
						   in_hit->sector);
	
	if(in_hit->end == DBCALGeometry::kUpstream) {
	    // handle the upstream end
	    return the_table.at(the_cell).first;
	    //return the_table[the_cell].first;
	} else {
	    // handle the downstream end
	    return the_table.at(the_cell).second;
	    //return the_table[the_cell].second;
	}

}

const double DBCALTDCHit_factory::GetConstant( const bcal_digi_constants_t &the_table, 
					       const DBCALTDCDigiHit *in_digihit) const
{
    	char str[256];
	
	if( (in_digihit->module <= 0) || (in_digihit->module > static_cast<unsigned int>(BCAL_NUM_MODULES))) {
		sprintf(str, "Bad module # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->module, BCAL_NUM_MODULES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->layer <= 0) || (in_digihit->layer > static_cast<unsigned int>(BCAL_NUM_LAYERS))) {
		sprintf(str, "Bad layer # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->layer, BCAL_NUM_LAYERS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->sector <= 0) || (in_digihit->sector >static_cast<unsigned int>( BCAL_NUM_SECTORS))) {
		sprintf(str, "Bad sector # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->sector, BCAL_NUM_SECTORS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->end != DBCALGeometry::kUpstream) && (in_digihit->end != DBCALGeometry::kDownstream) ) {
	        sprintf(str, "Bad end # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 0-1", in_digihit->end);
		cerr << str << endl;
		throw JException(str);
	}

	const int the_cell = DBCALGeometry::cellId(in_digihit->module,
						   in_digihit->layer,
						   in_digihit->sector);
	
	if(in_digihit->end == DBCALGeometry::kUpstream) {
	    // handle the upstream end
	    return the_table.at(the_cell).first;
	} else {
	    // handle the downstream end
	    return the_table.at(the_cell).second;
	}

}
/*
const double DBCALTDCHit_factory::GetConstant( const bcal_digi_constants_t &the_table,
					    const DTranslationTable *ttab,
					    const int in_rocid, const int in_slot, const int in_channel) const {
    
	char str[256];
	
	DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
	DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);
	
	if( (channel_info.bcal.module <= 0) 
	    || (channel_info.bcal.module > static_cast<unsigned int>(BCAL_NUM_MODULES))) {
		sprintf(str, "Bad module # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", channel_info.bcal.module, BCAL_NUM_MODULES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (channel_info.bcal.layer <= 0) 
	    || (channel_info.bcal.layer > static_cast<unsigned int>(BCAL_NUM_LAYERS))) {
		sprintf(str, "Bad layer # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", channel_info.bcal.layer, BCAL_NUM_LAYERS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (channel_info.bcal.sector <= 0) 
	    || (channel_info.bcal.sector > static_cast<unsigned int>(BCAL_NUM_SECTORS))) {
		sprintf(str, "Bad sector # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 1-%d", channel_info.bcal.sector, BCAL_NUM_SECTORS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (channel_info.bcal.end != DBCALGeometry::kUpstream) 
	    && (channel_info.bcal.end != DBCALGeometry::kDownstream) ) {
		sprintf(str, "Bad end # requested in DBCALTDCHit_factory::GetConstant()! requested=%d , should be 0-1", channel_info.bcal.end);
		cerr << str << endl;
		throw JException(str);
	}

	int the_cell = DBCALGeometry::cellId(channel_info.bcal.module,
					     channel_info.bcal.layer,
					     channel_info.bcal.sector);
	
	if(channel_info.bcal.end == DBCALGeometry::kUpstream) {
	    // handle the upstream end
	    return the_table.at(the_cell).first;
	} else {
	    // handle the downstream end
	    return the_table.at(the_cell).second;
	}
}
*/
