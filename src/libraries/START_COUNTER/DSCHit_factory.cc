// $Id$
//
//    File: DSCHit_factory.cc
// Created: Tue Aug  6 12:53:32 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
using namespace std;

#include <START_COUNTER/DSCDigiHit.h>
#include <START_COUNTER/DSCTDCDigiHit.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250Config.h>
#include "DSCHit_factory.h"

using namespace jana;

bool DSCHit_fadc_cmp(const DSCDigiHit *a,const DSCDigiHit *b)
{
    if (a->sector==b->sector) return (a->pulse_time<b->pulse_time);
    return (a->sector<b->sector);
}

bool DSCHit_tdc_cmp(const DSCTDCDigiHit *a,const DSCTDCDigiHit *b)
{
    if (a->sector==b->sector) return (a->time<b->time);
    return (a->sector<b->sector);
}

//------------------
// init
//------------------
jerror_t DSCHit_factory::init(void)
{
    DELTA_T_ADC_TDC_MAX = 20.0; // ns
    //DELTA_T_ADC_TDC_MAX = 50.0; // ns
    //DELTA_T_ADC_TDC_MAX = 3600.0; // ns
    gPARMS->SetDefaultParameter("SC:DELTA_T_ADC_TDC_MAX", DELTA_T_ADC_TDC_MAX,
            "Maximum difference in ns between a (calibrated) fADC time and"
            " F1TDC time for them to be matched in a single hit");

    HIT_TIME_WINDOW = 60.0; //ns
    gPARMS->SetDefaultParameter("SC:HIT_TIME_WINDOW", HIT_TIME_WINDOW,
           "Time window of trigger corrected TDC time in which a hit in"
	   " in the TDC will match to a hit in the fADC to form an ST hit");
			     
    //ADC_THRESHOLD = 200.; // adc counts (= 50 mV threshold)
    ADC_THRESHOLD = 120.; // adc counts (= 10 Mv threshold)
    gPARMS->SetDefaultParameter("SC:ADC_THRESHOLD", ADC_THRESHOLD,
            "Software pulse integral threshold");

    /// set the base conversion scales
    a_scale    = 0.0001; 
    t_scale    = 0.0625;   // 62.5 ps/count
    t_base     = 0.;       // ns
    t_tdc_base = 0.;

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DSCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
    // Only print messages for one thread whenever run number change
    static pthread_mutex_t print_mutex = PTHREAD_MUTEX_INITIALIZER;
    static set<int> runs_announced;
    pthread_mutex_lock(&print_mutex);
    bool print_messages = false;
    if(runs_announced.find(runnumber) == runs_announced.end()){
        print_messages = true;
        runs_announced.insert(runnumber);
    }
    pthread_mutex_unlock(&print_mutex);

    /// Read in calibration constants
    if(print_messages) jout << "In DSCHit_factory, loading constants..." << endl;

    // load scale factors
    map<string,double> scale_factors;
    // a_scale (SC_ADC_SCALE)
    if (eventLoop->GetCalib("/START_COUNTER/digi_scales", scale_factors))
        jout << "Error loading /START_COUNTER/digi_scales !" << endl;
    if (scale_factors.find("SC_ADC_ASCALE") != scale_factors.end())
        a_scale = scale_factors["SC_ADC_ASCALE"];
    else
        jerr << "Unable to get SC_ADC_ASCALE from /START_COUNTER/digi_scales !" 
            << endl;
    // t_scale (SC_ADC_SCALE)
    if (scale_factors.find("SC_ADC_TSCALE") != scale_factors.end())
        t_scale = scale_factors["SC_ADC_TSCALE"];
    else
        jerr << "Unable to get SC_ADC_TSCALE from /START_COUNTER/digi_scales !" 
            << endl;

    // load base time offset
    map<string,double> base_time_offset;
    // t_base (SC_BASE_TIME_OFFSET)
    if (eventLoop->GetCalib("/START_COUNTER/base_time_offset",base_time_offset))
        jout << "Error loading /START_COUNTER/base_time_offset !" << endl;
    if (base_time_offset.find("SC_BASE_TIME_OFFSET") != base_time_offset.end())
        t_base = base_time_offset["SC_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get SC_BASE_TIME_OFFSET from /START_COUNTER/base_time_offset !" << endl;
    // t_tdc_base (SC_TDC_BASE_TIME_OFFSET)
    if (base_time_offset.find("SC_TDC_BASE_TIME_OFFSET") != base_time_offset.end())
        t_tdc_base = base_time_offset["SC_TDC_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get SC_BASE_TIME_OFFSET from /START_COUNTER/base_time_offset !" << endl;

    // load constant tables
    // a_gains (gains)
    if (eventLoop->GetCalib("/START_COUNTER/gains", a_gains))
        jout << "Error loading /START_COUNTER/gains !" << endl;
    // a_pedestals (pedestals)
    if (eventLoop->GetCalib("/START_COUNTER/pedestals", a_pedestals))
        jout << "Error loading /START_COUNTER/pedestals !" << endl;
    // adc_time_offsets (adc_timing_offsets)
    if (eventLoop->GetCalib("/START_COUNTER/adc_timing_offsets", adc_time_offsets))
        jout << "Error loading /START_COUNTER/adc_timing_offsets !" << endl;
    // tdc_time_offsets (tdc_timing_offsets)
    if (eventLoop->GetCalib("/START_COUNTER/tdc_timing_offsets", tdc_time_offsets))
        jout << "Error loading /START_COUNTER/tdc_timing_offsets !" << endl;
    // timewalk_parameters (timewalk_parms)
    if(eventLoop->GetCalib("START_COUNTER/timewalk_parms", timewalk_parameters))
        jout << "Error loading /START_COUNTER/timewalk_parms !" << endl;


    /* 
    // load higher order corrections
    map<string,double> in_prop_corr;
    map<string,double> in_atten_corr;

    if(!eventLoop->GetCalib("/START_COUNTER/propagation_speed",in_prop_corr ))
    jout << "Error loading /START_COUNTER/propagation_speed !" << endl;
    if(!eventLoop->GetCalib("/START_COUNTER/attenuation_factor", in_atten_corr))
    jout << "Error loading /START_COUNTER/attenuation_factor !" << endl;

    // propogation correction:  A + Bx
    propogation_corr_factors.push_back(in_prop_corr["A"]);
    propogation_corr_factors.push_back(in_prop_corr["B"]);

    // attenuation correction:  A + Bx + Cx^2 + Dx^3 + Ex^4 + Fx^5
    attenuation_corr_factors.push_back(in_atten_corr["A"]);
    attenuation_corr_factors.push_back(in_atten_corr["B"]);
    attenuation_corr_factors.push_back(in_atten_corr["C"]);
    attenuation_corr_factors.push_back(in_atten_corr["D"]);
    attenuation_corr_factors.push_back(in_atten_corr["E"]);
    attenuation_corr_factors.push_back(in_atten_corr["F"]);
    */

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DSCHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
    /// Generate DSCHit object for each DSCDigiHit object.
    /// This is where the first set of calibration constants
    /// is applied to convert from digitzed units into natural
    /// units.
    ///
    /// Note that this code does NOT get called for simulated
    /// data in HDDM format. The HDDM event source will copy
    /// the precalibrated values directly into the _data vector.

    // First, make hits out of all fADC250 hits
    vector<const DSCDigiHit*> digihits;
    loop->Get(digihits);
    sort(digihits.begin(),digihits.end(),DSCHit_fadc_cmp);

	const DTTabUtilities* locTTabUtilities = NULL;
	loop->GetSingle(locTTabUtilities);

    char str[256];

    for (unsigned int i = 0; i < digihits.size(); i++)
    {
        const DSCDigiHit *digihit = digihits[i];

        // There is a slight difference between Mode 7 and 8 data
        // The following condition signals an error state in the flash algorithm
        // Do not make hits out of these
        const Df250PulsePedestal* PPobj = NULL;
        digihit->GetSingle(PPobj);
        if (PPobj != NULL)
	  {
            if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
	  }

	// Make sure sector is in valid range
        if( (digihit->sector <= 0) && (digihit->sector > MAX_SECTORS)) 
        {
            sprintf(str, "DSCDigiHit sector out of range! sector=%d (should be 1-%d)", 
                    digihit->sector, MAX_SECTORS);
            throw JException(str);
        }

        // Initialize pedestal to one found in CCDB, but override it
        // with one found in event if is available
        double pedestal = a_pedestals[digihit->sector-1];
        const Df250PulseIntegral *pulse_integral = NULL;
        digihit->GetSingle(pulse_integral);
        if (pulse_integral != NULL) 
	  {
            double single_sample_ped = (double)pulse_integral->pedestal;
            double nsamples_integral = (double)pulse_integral->nsamples_integral;
            double nsamples_pedestal = (double)pulse_integral->nsamples_pedestal;
            pedestal          = single_sample_ped * nsamples_integral/nsamples_pedestal;
	  }      	

        // Apply calibration constants here
        double A = (double)digihit->pulse_integral;
        double T = (double)digihit->pulse_time;
        double dA = A - pedestal;
	double ped_corr_pulse_peak = PPobj->pulse_peak - PPobj->pedestal;

        if (digihit->pulse_time == 0) continue; // Should already be caught, but I'll leave it
        //if (dA < ADC_THRESHOLD) continue; // Will comment out until this is set to something useful by default

        DSCHit *hit = new DSCHit;
        hit->sector = digihit->sector;

        // Sectors are numbered from 1-30
        hit->dE = dA; // This will be scaled to energy units later
        hit->t_fADC = t_scale * T - adc_time_offsets[hit->sector-1] + t_base;
        hit->t_TDC = numeric_limits<double>::quiet_NaN();

        hit->has_TDC = false;
        hit->has_fADC = true;

        hit->t = hit->t_fADC; // set time from fADC in case no TDC hit
	hit->pulse_height = ped_corr_pulse_peak;

        // add in higher order corrections?

        hit->AddAssociatedObject(digihit);

        _data.push_back(hit);
    }


    // Get the trigger time from the f1 TDC
    vector<const DF1TDCHit*> tdchit;
    eventLoop->Get(tdchit);

	// Next, loop over TDC hits, matching them to the
	// existing fADC hits where possible and updating
	// their time information. If no match is found, then
	// create a new hit with just the TDC info.
	vector<const DSCTDCDigiHit*> tdcdigihits;
	loop->Get(tdcdigihits);
	sort(tdcdigihits.begin(),tdcdigihits.end(),DSCHit_tdc_cmp);

	for (unsigned int i = 0; i < tdcdigihits.size(); i++)
	  {
	    const DSCTDCDigiHit *digihit = tdcdigihits[i];

	    // Make sure sector is in valid range
	    if((digihit->sector <= 0) && (digihit->sector > MAX_SECTORS))
	      {
		sprintf(str, "DSCDigiHit sector out of range! sector=%d (should be 1-%d)",
			digihit->sector, MAX_SECTORS);
		throw JException(str);
	      }

	    unsigned int id = digihit->sector - 1;
	    double T = locTTabUtilities->Convert_DigiTimeToNs_F1TDC(digihit) - tdc_time_offsets[id] + t_tdc_base;

	    // cout << "T = " << T << endl;
	    // jout << "T = " << T << endl;

	    // Look for existing hits to see if there is a match
	    //   or create new one if there is no match
	    // Require that the trigger corrected TDC time fall within
	    //   a reasonable time window so that when a hit is associated with
	    //   a hit in the TDC and not the ADC it is a "decent" TDC hit
	    if (fabs(T) < HIT_TIME_WINDOW)
	      {
		//jout << " T cut = " << T << endl;
		DSCHit *hit = FindMatch(digihit->sector, T);
		if (! hit)
		  {
		    hit = new DSCHit;
		    hit->sector = digihit->sector;
		    hit->dE = 0.0;
		    hit->t_fADC= numeric_limits<double>::quiet_NaN();
		    hit->has_fADC=false;
		    _data.push_back(hit);
		  }

		hit->has_TDC=true;
		hit->t_TDC=T;
		//jout << "t_tDC = " << hit->t_TDC << endl;


		if (hit->dE > 0.0)
		  {
		    // // Correct for time walk
		    // // The correction is the form t=t_tdc- C1 (A^C2 - A0^C2)
		    // double A  = hit->dE;
		    // double C1 = timewalk_parameters[id][1];
		    // double C2 = timewalk_parameters[id][2];
		    // double A0 = timewalk_parameters[id][3];
		    // T -= C1*(pow(A,C2) - pow(A0,C2));

		    // Correct for timewalk using pulse peak instead of pulse integral
		    double A        = hit->pulse_height;
		    double C1       = timewalk_parameters[id][0];
		    double C2       = timewalk_parameters[id][1];
		    double A_THRESH = timewalk_parameters[id][2];
		    double A0       = timewalk_parameters[id][3];
		    T -= C1*(pow(A/A_THRESH, C2) - pow(A0/A_THRESH, C2));
		  }
		hit->t=T;
		//jout << " T cut TW Corr = " << T << endl;
		hit->AddAssociatedObject(digihit);
	      } // Hit time window cut

	}


    // Apply calibration constants to convert pulse integrals to energy 
    // units
    for (unsigned int i=0;i<_data.size();i++)
      {
        _data[i]->dE*=a_scale * a_gains[_data[i]->sector-1];
      }


    return NOERROR;
}

//------------------
// FindMatch
//------------------
DSCHit* DSCHit_factory::FindMatch(int sector, double T)
{
    DSCHit *best_match=NULL;

    // Loop over existing hits (from fADC) and look for a match
    // in both the sector and the time.
    for(unsigned int i = 0; i < _data.size(); i++)
    {
        DSCHit *hit = _data[i];

        if (! isfinite(hit->t_fADC))
	  continue; // only match to fADC hits, not bachelor TDC hits

        if (hit->sector != sector)
	  continue; // require identical sectors fired

        double delta_T = fabs(hit->t - T);
        if (delta_T > DELTA_T_ADC_TDC_MAX)
	  continue;
	  
        if (best_match != NULL)
        {
            if (delta_T < fabs(best_match->t - T))
                best_match = hit;
        } else best_match = hit;
    }

    return best_match;
}

//------------------
// erun
//------------------
jerror_t DSCHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DSCHit_factory::fini(void)
{
    return NOERROR;
}


//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DSCHit_factory::GetConstant(const vector<double> &the_table,
        const int in_sector) const
{
    char str[256];

    if ( (in_sector < 0) || (in_sector >= MAX_SECTORS)) 
    {
        sprintf(str, "Bad sector # requested in DSCHit_factory::GetConstant()!"
                " requested=%d , should be %ud", in_sector, MAX_SECTORS);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_sector];
}

const double DSCHit_factory::GetConstant(const vector<double> &the_table,
        const DSCDigiHit *in_digihit) const 
{
    char str[256];

    if ( (in_digihit->sector < 0) || (in_digihit->sector >= MAX_SECTORS)) 
    {
        sprintf(str, "Bad sector # requested in DSCHit_factory::GetConstant()!"
                " requested=%d , should be %ud", 
                in_digihit->sector, MAX_SECTORS);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_digihit->sector];
}

const double DSCHit_factory::GetConstant(const vector<double> &the_table,
        const DSCHit *in_hit) const 
{

    char str[256];

    if ( (in_hit->sector < 0) || (in_hit->sector >= MAX_SECTORS)) 
    {
        sprintf(str, "Bad sector # requested in DSCHit_factory::GetConstant()!"
                " requested=%d , should be %ud",
                in_hit->sector, MAX_SECTORS);
        cerr << str << endl;
        throw JException(str);
    }

    return the_table[in_hit->sector];
}
/*
   const double DSCHit_factory::GetConstant(const vector<double> &the_table,
   const DTranslationTable *ttab,
   const int in_rocid,
   const int in_slot, 
   const int in_channel) const
   {
   char str[256];

   DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
   DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);

   if( (channel_info.sc.sector <= 0) 
   || (channel_info.sc.sector > static_cast<unsigned int>(MAX_SECTORS))) {
   sprintf(str, "Bad sector # requested in DSCHit_factory::GetConstant()!"
   " requested=%d , should be %ud",
   channel_info.sc.sector, MAX_SECTORS);
   cerr << str << endl;
   throw JException(str);
   }

   return the_table[channel_info.sc.sector];
   }
   */
