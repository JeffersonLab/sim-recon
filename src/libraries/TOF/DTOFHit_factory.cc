// $Id$
//
//    File: DTOFHit_factory.cc
// Created: Wed Aug  7 09:30:17 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>
using namespace std;

#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include "DTOFHit_factory.h"
using namespace jana;


//------------------
// init
//------------------
jerror_t DTOFHit_factory::init(void)
{
	DELTA_T_ADC_TDC_MAX = 4.0; // ns
	gPARMS->SetDefaultParameter("TOF:DELTA_T_ADC_TDC_MAX", DELTA_T_ADC_TDC_MAX, "Maximum difference in ns between a (calibrated) fADC time and F1TDC time for them to be matched in a single hit");

	/// Set basic conversion constants
	a_scale    = 0.2/5.2E5;
	t_scale    = 0.0625;   // 62.5 ps/count
	tdc_scale  = 0.025;    // 25 ps/count
        t_min      = -100.;    // ns

	TOF_NUM_PLANES = 2;
	TOF_NUM_BARS = 44;

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTOFHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
    
        // read in geometry information
        vector<const DTOFGeometry*> tofGeomVect;
	eventLoop->Get( tofGeomVect );
	if(tofGeomVect.size()<1)  return OBJECT_NOT_AVAILABLE;
	const DTOFGeometry& tofGeom = *(tofGeomVect[0]);
	
        /// Read in calibration constants
        vector<double> raw_adc_pedestals;
        vector<double> raw_adc_gains;
        vector<double> raw_adc_offsets;
        vector<double> raw_tdc_offsets;
        vector<double> raw_tdc_scales;

        jout << "In DTOFHit_factory, loading constants..." << endl;

	// load scale factors
	map<string,double> scale_factors;
	if(eventLoop->GetCalib("/TOF/digi_scales", scale_factors))
	    jout << "Error loading /TOF/digi_scales !" << endl;
	if( scale_factors.find("TOF_ADC_ASCALE") != scale_factors.end() ) {
		;	//a_scale = scale_factors["TOF_ADC_ASCALE"];
	} else {
	    jerr << "Unable to get TOF_ADC_ASCALE from /TOF/digi_scales !" << endl;
	}
	if( scale_factors.find("TOF_ADC_TSCALE") != scale_factors.end() ) {
		; //t_scale = scale_factors["TOF_ADC_TSCALE"];
	} else {
	    jerr << "Unable to get TOF_ADC_TSCALE from /TOF/digi_scales !" << endl;
	}
	if( scale_factors.find("TOF_TDC_SCALE") != scale_factors.end() ) {
		; //tdc_scale = scale_factors["TOF_TDC_SCALE"];
	} else {
	    jerr << "Unable to get TOF_TDC_SCALE from /TOF/digi_scales !" << endl;
	}

        if(eventLoop->GetCalib("TOF/pedestals", raw_adc_pedestals))
	    jout << "Error loading /TOF/pedestals !" << endl;
        if(eventLoop->GetCalib("TOF/gains", raw_adc_gains))
	    jout << "Error loading /TOF/gains !" << endl;
	if(eventLoop->GetCalib("TOF/timing_scales", raw_tdc_scales))
	    jout << "Error loading /TOF/timing_scales !" << endl;
	if(eventLoop->GetCalib("TOF/adc_timing_offsets", raw_adc_offsets))
	    jout << "Error loading /TOF/adc_timing_offsets !" << endl;
	if(eventLoop->GetCalib("TOF/timing_offsets", raw_tdc_offsets))
	    jout << "Error loading /TOF/timing_offsets !" << endl;

        FillCalibTable(adc_pedestals, raw_adc_pedestals, tofGeom);
        FillCalibTable(adc_gains, raw_adc_gains, tofGeom);
        FillCalibTable(adc_time_offsets, raw_adc_offsets, tofGeom);
        FillCalibTable(tdc_time_offsets, raw_tdc_offsets, tofGeom);
        FillCalibTable(tdc_scales, raw_tdc_scales, tofGeom);
/*
	CheckCalibTable(adc_pedestals,"/TOF/pedestals");
	CheckCalibTable(adc_gains,"/TOF/gains");
	CheckCalibTable(adc_time_offsets,"/TOF/adc_timing_offsets");
	CheckCalibTable(tdc_time_offsets,"/TOF/timing_offsets");
	CheckCalibTable(tdc_scales,"/TOF/timing_scales");
*/
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTOFHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
	/// Generate DTOFHit object for each DTOFDigiHit object.
	/// This is where the first set of calibration constants
	/// is applied to convert from digitzed units into natural
	/// units.
	///
	/// Note that this code does NOT get called for simulated
	/// data in HDDM format. The HDDM event source will copy
	/// the precalibrated values directly into the _data vector.

	// First, make hits out of all fADC250 hits
	vector<const DTOFDigiHit*> digihits;
	loop->Get(digihits);
	for(unsigned int i=0; i<digihits.size(); i++){
		const DTOFDigiHit *digihit = digihits[i];

		// Get pedestal.  Prefer associated event pedestal if it exists.
		// Otherwise, use the average pedestal from CCDB
		double pedestal = GetConstant(adc_pedestals,digihit);
		vector<const Df250PulseIntegral*> PIvect;
		digihit->Get(PIvect);
		if(!PIvect.empty()){
		    const Df250PulseIntegral *PIobj = PIvect[0];
		    pedestal = PIobj->pedestal;
		}

		DTOFHit *hit = new DTOFHit;
		hit->plane = digihit->plane;
		hit->bar   = digihit->bar;
		hit->end   = digihit->end;
		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;


		hit->dE = a_scale * (A - pedestal);
		hit->t = t_scale * (T - GetConstant(adc_time_offsets, digihit)) + t_min;
		hit->sigma_t = 4.0;    // ns (what is the fADC time resolution?)
		hit->has_fADC = true;
		hit->has_TDC  = false; // will get set to true below if appropriate

/*
		cout << "TOF ADC hit =  (" << hit->plane << "," << hit->bar << "," << hit->end << ")  " 
		     << t_scale << " " << T << "  "
		     << GetConstant(adc_time_offsets, digihit) << " " 
		     << t_scale*GetConstant(adc_time_offsets, digihit) << " " << hit->t << endl;
*/
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	// Second, loop over TDC hits, matching them to the
	// existing fADC hits where possible and updating
	// their time information. If no match is found, then
	// create a new hit with just the TDC info.
	vector<const DTOFTDCDigiHit*> tdcdigihits;
	loop->Get(tdcdigihits);
	for(unsigned int i=0; i<tdcdigihits.size(); i++){
		const DTOFTDCDigiHit *digihit = tdcdigihits[i];

		// Apply calibration constants here
		double T = (double)digihit->time;

		T = tdc_scale * (T - GetConstant(tdc_time_offsets, digihit)) + t_min;
		// future: allow for seperate TDC scales for each channel
		//T = GetConstant(tdc_scales, digihit)
		//  * (T - GetConstant(tdc_time_offsets, digihit));

/*
		cout << "TOF TDC hit =  (" << digihit->plane << "," << digihit->bar << "," << digihit->end << ")  " 
		     << tdc_scale << " " << T << "  "
		     << GetConstant(tdc_time_offsets, digihit) << " " 
		     << tdc_scale*GetConstant(tdc_time_offsets, digihit) << " " << T << endl;
*/

		// add in timewalk corrections here
		
		// Look for existing hits to see if there is a match
		// or create new one if there is no match
		DTOFHit *hit = FindMatch(digihit->plane, digihit->bar, digihit->end, T);  
		if(!hit){
			hit = new DTOFHit;
			hit->plane = digihit->plane;
			hit->bar   = digihit->bar;
			hit->end   = digihit->end;
			hit->dE = 0.0;
			hit->has_fADC = false;

			_data.push_back(hit);
		}
		
		hit->t = T;
		hit->sigma_t = 0.160;    // ns (what is the TOF TDC time resolution?)
		hit->has_TDC = true;
		
		hit->AddAssociatedObject(digihit);
	}

	return NOERROR;
}

//------------------
// FindMatch
//------------------
DTOFHit* DTOFHit_factory::FindMatch(int plane, int bar, int end, double T)
{
	DTOFHit* best_match = NULL;

	// Loop over existing hits (from fADC) and look for a match
	// in both the sector and the time.
	for(unsigned int i=0; i<_data.size(); i++){
		DTOFHit *hit = _data[i];
		
		if(!hit->has_fADC) continue; // only match to fADC hits, not bachelor TDC hits
		if(hit->plane != plane) continue;
		if(hit->bar != bar) continue;
		if(hit->end != end) continue;
		
		double delta_T = fabs(hit->t - T);
		if(delta_T > DELTA_T_ADC_TDC_MAX) continue;
		
		// if there are multiple hits, pick the one that is closest in time
		if(best_match != NULL) {
			if(delta_T < fabs(best_match->t - T))
				best_match = hit;
		} else {
			best_match = hit;
		}
	}
	
	return best_match;
}

//------------------
// erun
//------------------
jerror_t DTOFHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTOFHit_factory::fini(void)
{
	return NOERROR;
}


//------------------
// FillCalibTable
//------------------
void DTOFHit_factory::FillCalibTable(tof_digi_constants_t &table, vector<double> &raw_table, 
				     const DTOFGeometry &tofGeom)
{
    char str[256];
    int channel = 0;
 
    // reset the table before filling it
    table.clear();

    for(int plane=0; plane<tofGeom.NLAYERS; plane++) {
	table.push_back( vector< pair<double,double> >(TOF_NUM_BARS) );
	for(int bar=0; bar<tofGeom.NLONGBARS; bar++) {
	    if( (channel > TOF_MAX_CHANNELS) || (channel+1 > TOF_MAX_CHANNELS) ) {  // sanity check
		sprintf(str, "Too many channels for TOF table! channel=%d (should be %d)", 
			channel, TOF_MAX_CHANNELS);
		cerr << str << endl;
		throw JException(str);
	    }
	    
	    table[plane][bar] = pair<double,double>(raw_table[channel],raw_table[channel+1]);
	    
	    channel += 2;
	}
    }

    // check to make sure that we loaded enough channels
    if(channel != TOF_MAX_CHANNELS) { 
	sprintf(str, "Not enough channels for TOF table! channel=%d (should be %d)", 
		channel, TOF_MAX_CHANNELS);
	cerr << str << endl;
	throw JException(str);
    }
}


//------------------------------------
// GetConstant
//   Allow a few different interfaces
//   NOTE: LoadGeometry() must be called before calling these functions
//
//   TOF Geometry as defined in the Translation Table:
//     plane = 0-1
//     bar   = 1-44
//     end   = 0-1
//   Note the different counting schemes used
//------------------------------------
const double DTOFHit_factory::GetConstant( const tof_digi_constants_t &the_table, 
					    const int in_plane, const int in_bar, const int in_end) const
{
    	char str[256];
	
	if( (in_plane < 0) || (in_plane > TOF_NUM_PLANES)) {
		sprintf(str, "Bad module # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_plane, TOF_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_bar <= 0) || (in_bar > TOF_NUM_BARS)) {
		sprintf(str, "Bad layer # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_bar, TOF_NUM_BARS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_end != 0) && (in_end != 1) ) {
		sprintf(str, "Bad end # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 0-1", in_end);
		cerr << str << endl;
		throw JException(str);
	}

	// we have two ends, indexed as 0/1 
	// could be north/south or up/down depending on the bar orientation
	if(in_end == 0) {
	    return the_table[in_plane][in_bar].first;
	} else {
	    return the_table[in_plane][in_bar].second;
	}
}

const double DTOFHit_factory::GetConstant( const tof_digi_constants_t &the_table, 
					    const DTOFHit *in_hit) const
{
    	char str[256];
	
	if( (in_hit->plane < 0) || (in_hit->plane > TOF_NUM_PLANES)) {
		sprintf(str, "Bad module # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_hit->plane, TOF_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_hit->bar <= 0) || (in_hit->bar > TOF_NUM_BARS)) {
		sprintf(str, "Bad layer # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_hit->bar, TOF_NUM_BARS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_hit->end != 0) && (in_hit->end != 1) ) {
		sprintf(str, "Bad end # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 0-1", in_hit->end);
		cerr << str << endl;
		throw JException(str);
	}

	// we have two ends, indexed as 0/1 
	// could be north/south or up/down depending on the bar orientation
	if(in_hit->end == 0) {
	    return the_table[in_hit->plane][in_hit->bar].first;
	} else {
	    return the_table[in_hit->plane][in_hit->bar].second;
	}
}

const double DTOFHit_factory::GetConstant( const tof_digi_constants_t &the_table, 
					    const DTOFDigiHit *in_digihit) const
{
    	char str[256];
	
	if( (in_digihit->plane < 0) || (in_digihit->plane > TOF_NUM_PLANES)) {
		sprintf(str, "Bad module # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->plane, TOF_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->bar <= 0) || (in_digihit->bar > TOF_NUM_BARS)) {
		sprintf(str, "Bad layer # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->bar, TOF_NUM_BARS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->end != 0) && (in_digihit->end != 1) ) {
		sprintf(str, "Bad end # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 0-1", in_digihit->end);
		cerr << str << endl;
		throw JException(str);
	}

	// we have two ends, indexed as 0/1 
	// could be north/south or up/down depending on the bar orientation
	if(in_digihit->end == 0) {
	    return the_table[in_digihit->plane][in_digihit->bar].first;
	} else {
	    return the_table[in_digihit->plane][in_digihit->bar].second;
	}
}

const double DTOFHit_factory::GetConstant( const tof_digi_constants_t &the_table, 
					    const DTOFTDCDigiHit *in_digihit) const
{
    	char str[256];
	
	if( (in_digihit->plane < 0) || (in_digihit->plane > TOF_NUM_PLANES)) {
		sprintf(str, "Bad module # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->plane, TOF_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->bar <= 0) || (in_digihit->bar > TOF_NUM_BARS)) {
		sprintf(str, "Bad layer # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", in_digihit->bar, TOF_NUM_BARS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (in_digihit->end != 0) && (in_digihit->end != 1) ) {
		sprintf(str, "Bad end # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 0-1", in_digihit->end);
		cerr << str << endl;
		throw JException(str);
	}

	// we have two ends, indexed as 0/1 
	// could be north/south or up/down depending on the bar orientation
	if(in_digihit->end == 0) {
	    return the_table[in_digihit->plane][in_digihit->bar].first;
	} else {
	    return the_table[in_digihit->plane][in_digihit->bar].second;
	}
}


/*
const double DTOFHit_factory::GetConstant( const tof_digi_constants_t &the_table,
					    const DTranslationTable *ttab,
					    const int in_rocid, const int in_slot, const int in_channel) const {
    
	char str[256];
	
	DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
	DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);
	
	if( (channel_info.tof.plane <= 0) 
	    || (channel_info.tof.plane > static_cast<unsigned int>(TOF_NUM_PLANES))) {
		sprintf(str, "Bad plane # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", channel_info.tof.plane, TOF_NUM_PLANES);
		cerr << str << endl;
		throw JException(str);
	}
	if( (channel_info.tof.bar <= 0) 
	    || (channel_info.tof.bar > static_cast<unsigned int>(TOF_NUM_BARS))) {
		sprintf(str, "Bad bar # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 1-%d", channel_info.tof.bar, TOF_NUM_BARS);
		cerr << str << endl;
		throw JException(str);
	}
	if( (channel_info.tof.end != 0) && (channel_info.tof.end != 1) ) {
		sprintf(str, "Bad end # requested in DTOFHit_factory::GetConstant()! requested=%d , should be 0-1", channel_info.tof.end);
		cerr << str << endl;
		throw JException(str);
	}

	int the_cell = DTOFGeometry::cellId(channel_info.tof.module,
					     channel_info.tof.layer,
					     channel_info.tof.sector);
	
	if(channel_info.tof.end == DTOFGeometry::kUpstream) {
	    // handle the upstream end
	    return the_table.at(the_cell).first;
	} else {
	    // handle the downstream end
	    return the_table.at(the_cell).second;
	}
}
*/
