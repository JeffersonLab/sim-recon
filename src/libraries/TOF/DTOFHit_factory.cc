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
#include <limits>

using namespace std;

#include <TOF/DTOFDigiHit.h>
#include <TOF/DTOFTDCDigiHit.h>
#include "DTOFHit_factory.h"
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250Config.h>
#include <DAQ/DCODAROCInfo.h>
using namespace jana;

static bool COSMIC_DATA = false;

//------------------
// init
//------------------
jerror_t DTOFHit_factory::init(void)
{
  DELTA_T_ADC_TDC_MAX = 10.0; // ns
  //	DELTA_T_ADC_TDC_MAX = 30.0; // ns, value based on the studies from cosmic events
	gPARMS->SetDefaultParameter("TOF:DELTA_T_ADC_TDC_MAX", DELTA_T_ADC_TDC_MAX, "Maximum difference in ns between a (calibrated) fADC time and F1TDC time for them to be matched in a single hit");
	
	int analyze_cosmic_data = 0;
	gPARMS->SetDefaultParameter("TOF:COSMIC_DATA", analyze_cosmic_data,
				    "Special settings for analysing cosmic data");
	if(analyze_cosmic_data > 0)
		COSMIC_DATA = true;


	/// Set basic conversion constants
	a_scale    = 0.2/5.2E5;
	t_scale    = 0.0625;   // 62.5 ps/count
	tdc_scale  = 0.025;    // 25 ps/count
        t_base     = 0.;       // ns
	t_base_tdc = 0.; // ns

	if(COSMIC_DATA)
		// Hardcoding of 110 taken from cosmics events
		tdc_adc_time_offset = 110.;
	else 
		tdc_adc_time_offset = 0.;

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

	// load base time offset
	map<string,double> base_time_offset;
	if (eventLoop->GetCalib("/TOF/base_time_offset",base_time_offset))
		jout << "Error loading /TOF/base_time_offset !" << endl;
	if (base_time_offset.find("TOF_BASE_TIME_OFFSET") != base_time_offset.end())
		t_base = base_time_offset["TOF_BASE_TIME_OFFSET"];
	else
		jerr << "Unable to get TOF_BASE_TIME_OFFSET from /TOF/base_time_offset !" << endl;	

	if (base_time_offset.find("TOF_TDC_BASE_TIME_OFFSET") != base_time_offset.end())
		t_base_tdc = base_time_offset["TOF_TDC_BASE_TIME_OFFSET"];
	else
		jerr << "Unable to get TOF_TDC_BASE_TIME_OFFSET from /TOF/base_time_offset !" << endl;
	
	// load constant tables
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
	if(eventLoop->GetCalib("TOF/timewalk_parms", timewalk_parameters))
	    jout << "Error loading /TOF/timewalk_parms !" << endl;

        FillCalibTable(adc_pedestals, raw_adc_pedestals, tofGeom);
        FillCalibTable(adc_gains, raw_adc_gains, tofGeom);
        FillCalibTable(adc_time_offsets, raw_adc_offsets, tofGeom);
        FillCalibTable(tdc_time_offsets, raw_tdc_offsets, tofGeom);
        FillCalibTable(tdc_scales, raw_tdc_scales, tofGeom);



	// load shift factors (only for fall 2014 runs)
	map<string,double> tof_tdc_shift;
	tdc_shift = -1;

	if(eventLoop->GetCalib("/TOF/tdc_shift", tof_tdc_shift))
	    jout << "Error loading /TOF/tdc_shift !" << endl;
	if( tof_tdc_shift.find("TOF_TDC_SHIFT") != tof_tdc_shift.end() ) {
	  tdc_shift = tof_tdc_shift["TOF_TDC_SHIFT"];
	  cout << "getting tdc_shift" << endl << endl;
	  cout << "tdc_shift = " << tdc_shift << endl << endl;
	} else {
	    jerr << "Unable to get TOF_TDC_SHIFT from /TOF/tdc_shift !" << endl;
	}

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
		const Df250PulseIntegral* PIobj = NULL;
		const Df250Config *configObj = NULL;
		digihit->GetSingle(PIobj);
		PIobj->GetSingle(configObj);
		if ((PIobj != NULL) && (configObj != NULL)) {
			// the measured pedestal must be scaled by the ratio of the number
			// of samples used to calculate the pedestal and the actual pulse
			pedestal = static_cast<double>(configObj->NSA_NSB) * PIobj->pedestal;
		}

		
		// Apply calibration constants here
		double A = (double)digihit->pulse_integral;
		double T = (double)digihit->pulse_time;
		T =  t_scale * T - GetConstant(adc_time_offsets, digihit) + t_base;
		pedestal=double(digihit->pedestal*digihit->nsamples_integral
				/digihit->nsamples_pedestal);
		double dA=A-pedestal;

		if (digihit->pulse_time==0) continue;
		if (dA<0) continue;

		DTOFHit *hit = new DTOFHit;
		hit->plane = digihit->plane;
		hit->bar   = digihit->bar;
		hit->end   = digihit->end;
		hit->dE=dA;  // this will be scaled to energy units later

		if(COSMIC_DATA)
		  hit->dE = (A - 55*pedestal); // value of 55 is taken from (NSB,NSA)=(10,45) in the confg file

		hit->t_TDC=numeric_limits<double>::quiet_NaN();
		hit->t_fADC=T;
		hit->t = hit->t_fADC;  // set initial time to the ADC time, in case there's no matching TDC hit

		hit->has_fADC=true;
		hit->has_TDC=false;

/*
		cout << "TOF ADC hit =  (" << hit->plane << "," << hit->bar << "," << hit->end << ")  " 
		     << t_scale << " " << T << "  "
		     << GetConstant(adc_time_offsets, digihit) << " " 
		     << t_scale*GetConstant(adc_time_offsets, digihit) << " " << hit->t << endl;
*/
		
		hit->AddAssociatedObject(digihit);
		
		_data.push_back(hit);
	}

	//Find the Trigger Time from the TI:
	// and determine which 4ns pulse the trigger is
	// within 24 ns.
	vector<const DCODAROCInfo*> TIInfo;
	loop->Get(TIInfo);
	
	unsigned long TriggerTime = 0;
	for (unsigned int k=0; k<TIInfo.size(); k++){
	  if (TIInfo[k]->rocid == 78){
	    TriggerTime = TIInfo[k]->timestamp;
	    break;
	  }
	}

	// TriggerBIT is not really a bit...
	int TriggerBIT = TriggerTime%6;  

	// Next, loop over TDC hits, matching them to the
	// existing fADC hits where possible and updating
	// their time information. If no match is found, then
	// create a new hit with just the TDC info.
	vector<const DTOFTDCDigiHit*> tdcdigihits;
	loop->Get(tdcdigihits);
	for(unsigned int i=0; i<tdcdigihits.size(); i++){
		const DTOFTDCDigiHit *digihit = tdcdigihits[i];

		// Apply calibration constants here
		double T = (double)digihit->time;

		// Get overall shift value for this run, shift for each value
		// of  TriggerBIT is given by tdc_shift - TriggerBIT
		// Shift will be positive for TDC times.

		// initalize to 0 in case there is no constant available
		int nshifts = 0;

		if(tdc_shift != -1){
		  int shift = tdc_shift - TriggerBIT;
		  if(shift < 0) shift += 6;
		  // TDC bins are 25 ps wide, so each block of 4 ns is 160 bins
		  nshifts = 160 * shift;
		}
		// cout << endl << "TriggerBIT = " << TriggerBIT << " nshifts = " << nshifts << endl << endl;
		// cout << "uncorrected: " << tdc_scale * (T - GetConstant(tdc_time_offsets, digihit)) + t_base_tdc + tdc_adc_time_offset << endl
		//      << "corrected  : " << tdc_scale * (T - GetConstant(tdc_time_offsets, digihit) + nshifts) + t_base_tdc + tdc_adc_time_offset << endl;
		  
		T = tdc_scale *T - GetConstant(tdc_time_offsets, digihit) 
		  + tdc_scale * nshifts
		  + t_base_tdc + tdc_adc_time_offset; 
		
 		// future: allow for seperate TDC scales for each channel
		//T = GetConstant(tdc_scales, digihit)
		//  * (T - GetConstant(tdc_time_offsets, digihit));

/*
		cout << "TOF TDC hit =  (" << digihit->plane << "," << digihit->bar << "," << digihit->end << ")  " 
		     << tdc_scale << " " << T << "  "
		     << GetConstant(tdc_time_offsets, digihit) << " " 
		     << tdc_scale*GetConstant(tdc_time_offsets, digihit) << " " << T << endl;
*/
		
		// Look for existing hits to see if there is a match
		// or create new one if there is no match
		DTOFHit *hit = FindMatch(digihit->plane, digihit->bar, digihit->end, T);
		//DTOFHit *hit = FindMatch(digihit->plane, hit->bar, hit->end, T);
		if(!hit){
			hit = new DTOFHit;
			hit->plane = digihit->plane;
			hit->bar   = digihit->bar;
			hit->end   = digihit->end;
			hit->dE = 0.0;
			hit->t_fADC=numeric_limits<double>::quiet_NaN();
			hit->has_fADC=false;

			_data.push_back(hit);
		}
		hit->has_TDC=true;
		hit->t_TDC=T;

		if (hit->dE>0.){
		  // time walk correction
		  // The correction is the form t=t_tdc- C1 (A^C2 - A0^C2)
		  int id=88*hit->plane+44*hit->end+hit->bar-1;
		  double A=hit->dE;
		  double C1=timewalk_parameters[id][1];
		  double C2=timewalk_parameters[id][2];
		  double A0=timewalk_parameters[id][3];
		  T-=C1*(pow(A,C2)-pow(A0,C2));
		}
		hit->t=T;

		hit->AddAssociatedObject(digihit);
	}

	// Apply calibration constants to convert pulse integrals to energy 
	// units
	for (unsigned int i=0;i<_data.size();i++){
	  _data[i]->dE*=a_scale;
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
		
		if(!isfinite(hit->t_fADC)) continue; // only match to fADC hits, not bachelor TDC hits
		if(hit->plane != plane) continue;
		if(hit->bar != bar) continue;
		if(hit->end != end) continue;
		
		//double delta_T = fabs(hit->t - T);
		double delta_T = fabs(T - hit->t);
		if(delta_T > DELTA_T_ADC_TDC_MAX) continue;

		return hit;
		
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
      int plane_index=2*TOF_NUM_BARS*plane;
      table.push_back( vector< pair<double,double> >(TOF_NUM_BARS) );
      for(int bar=0; bar<TOF_NUM_BARS; bar++) {
	table[plane][bar] 
	  = pair<double,double>(raw_table[plane_index+bar],
				raw_table[plane_index+TOF_NUM_BARS+bar]);
	channel+=2;	      
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
	    return the_table[in_hit->plane][in_hit->bar-1].first;
	} else {
	    return the_table[in_hit->plane][in_hit->bar-1].second;
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
	    return the_table[in_digihit->plane][in_digihit->bar-1].first;
	} else {
	    return the_table[in_digihit->plane][in_digihit->bar-1].second;
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
	    return the_table[in_digihit->plane][in_digihit->bar-1].first;
	} else {
	    return the_table[in_digihit->plane][in_digihit->bar-1].second;
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
