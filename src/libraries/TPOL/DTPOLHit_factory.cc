
#include <iostream>
#include <iomanip>
#include <limits>
#include <cmath>
using namespace std;

#include <TPOL/DTPOLSectorDigiHit.h>
#include <TPOL/DTPOLRingDigiHit.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250Config.h>
#include "DTPOLHit_factory.h"
#include "DAQ/DF1TDCHit.h"

using namespace jana;

// static consts need initialization
const double DTPOLHit_factory::SECTOR_DIVISION = 360. / DTPOLHit_factory::NSECTORS;
const double DTPOLHit_factory::INNER_RADIUS = 22. / 2;       // From "ACTIVE INNER DIAMETER" in catalog
const double DTPOLHit_factory::OUTER_RADIUS = 70. / 2;       // From "ACTIVE OUTER DIAMETER" in catalog
const double DTPOLHit_factory::RING_DIVISION   = (OUTER_RADIUS - INNER_RADIUS ) / DTPOLHit_factory::NRINGS;

// Order hits by sector/ring. If same sector/ring, use ADC time to order hits.
bool DTPOLSectorHit_fadc_cmp(const DTPOLSectorDigiHit *a,const DTPOLSectorDigiHit *b){
  if (a->sector==b->sector) return (a->pulse_time<b->pulse_time);
  return (a->sector<b->sector);
}

bool DTPOLRingHit_fadc_cmp(const DTPOLRingDigiHit *a,const DTPOLRingDigiHit *b){
  if (a->ring==b->ring) return (a->pulse_time<b->pulse_time);
  return (a->ring<b->ring);
}

//------------------
// init
//------------------
jerror_t DTPOLHit_factory::init(void)
{
    ADC_THRESHOLD = 120.;
    gPARMS->SetDefaultParameter("TPOL:ADC_THRESHOLD", ADC_THRESHOLD,
            "Software pulse integral threshold");

    /// set the base conversion scales
    a_scale    = 0.0001; 
    t_scale    = 0.0625;   // 62.5 ps/count

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTPOLHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
    /// Read in calibration constants
    jout << "In DTPOLHit_factory, loading constants..." << endl;

    // load scale factors
    map<string,double> scale_factors;
    // a_scale (SC_ADC_SCALE)
    if (eventLoop->GetCalib("/TPOL/digi_scales", scale_factors))
        jout << "Error loading /TPOL/digi_scales !" << endl;
    if (scale_factors.find("TPOL_ADC_ASCALE") != scale_factors.end())
        a_scale = scale_factors["TPOL_ADC_ASCALE"];
    else
        jerr << "Unable to get TPOL_ADC_ASCALE from /TPOL/digi_scales !" 
            << endl;
    // t_scale (SC_ADC_SCALE)
    if (scale_factors.find("TPOL_ADC_TSCALE") != scale_factors.end())
        t_scale = scale_factors["TPOL_ADC_TSCALE"];
    else
        jerr << "Unable to get TPOL_ADC_TSCALE from /TPOL/digi_scales !" 
            << endl;

    // load base time offset
    map<string,double> base_time_offset;
    if (eventLoop->GetCalib("/TPOL/base_time_offset",base_time_offset))
        jout << "Error loading /TPOL/base_time_offset !" << endl;
    if (base_time_offset.find("TPOL_BASE_TIME_OFFSET") != base_time_offset.end())
        t_base = base_time_offset["TPOL_BASE_TIME_OFFSET"];
    else
        jerr << "Unable to get TPOL_BASE_TIME_OFFSET from /TPOL/base_time_offset !" << endl;

    // load constant tables
    // a_gains (gains)
    if (eventLoop->GetCalib("/TPOL/gains", a_gains))
        jout << "Error loading /TPOL/gains !" << endl;
    // a_pedestals (pedestals)
    if (eventLoop->GetCalib("/TPOL/pedestals", a_pedestals))
        jout << "Error loading /TPOL/pedestals !" << endl;
    // adc_time_offsets (adc_timing_offsets)
    if (eventLoop->GetCalib("/TPOL/adc_timing_offsets", adc_time_offsets))
        jout << "Error loading /TPOL/adc_timing_offsets !" << endl;

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTPOLHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
    /// Generate DTPOLHit object for each DTPOLSectorDigiHit 
    /// and DTPOLRingDigiHit object.
    /// This is where the first set of calibration constants
    /// is applied to convert from digitzed units into natural
    /// units.
    ///
    /// Note that this code does NOT get called for simulated
    /// data in HDDM format. The HDDM event source will copy
    /// the precalibrated values directly into the _data vector.

    // First, get all fADC250 hits
    vector<const DTPOLSectorDigiHit*> sectordigihits;
    loop->Get(sectordigihits);
    sort(sectordigihits.begin(),sectordigihits.end(),DTPOLSectorHit_fadc_cmp);

    vector<const DTPOLRingDigiHit*> ringdigihits;
    loop->Get(ringdigihits);
    sort(ringdigihits.begin(),ringdigihits.end(),DTPOLRingHit_fadc_cmp);

    char str[256];

    // Loop over SECTOR hits
    for (unsigned int i = 0; i < sectordigihits.size(); i++){
      const DTPOLSectorDigiHit *sectordigihit = sectordigihits[i];

        // There is a slight difference between Mode 7 and 8 data
        // The following condition signals an error state in the flash algorithm
        // Do not make hits out of these
        const Df250PulsePedestal* PPobj = NULL;
        sectordigihit->GetSingle(PPobj);
        if (PPobj != NULL){
	  if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
	}

        // Make sure sector is in valid range
        if( sectordigihit->sector <= 0 || sectordigihit->sector > DTPOLHit_factory::NSECTORS){
	  sprintf(str, "DTPOLSectorDigiHit sector out of range! sector=%d (should be 1-%d)", 
		  sectordigihit->sector, DTPOLHit_factory::NSECTORS);
	  throw JException(str);
	}

        // Initialize pedestal to one found in CCDB, but override it
        // with one found in event if is available
        double pedestal = a_pedestals[sectordigihit->sector-1];
        const Df250PulseIntegral *pulse_integral = NULL;
        sectordigihit->GetSingle(pulse_integral);
        if (pulse_integral != NULL) {
	  double single_sample_ped = (double)pulse_integral->pedestal;
	  double nsamples_integral = (double)pulse_integral->nsamples_integral;
	  double nsamples_pedestal = (double)pulse_integral->nsamples_pedestal;
	  pedestal          = single_sample_ped * nsamples_integral/nsamples_pedestal;
	}      	

        // Apply calibration constants here
        double A = (double)sectordigihit->pulse_integral;
        double T = (double)sectordigihit->pulse_time;
        double dA = A - pedestal;

        if (sectordigihit->pulse_time == 0) continue;

        DTPOLHit *hit = new DTPOLHit;
        hit->sector = sectordigihit->sector;

        // Sectors are numbered from 1-30
        hit->dE = dA; // This will be scaled to energy units later
        hit->t_fADC = t_scale * T - adc_time_offsets[hit->sector-1] + t_base;

        hit->has_sector = true;

        hit->t = hit->t_fADC; // set time from fADC in case no TDC hit

        // add in higher order corrections?

        hit->AddAssociatedObject(sectordigihit);

        _data.push_back(hit);
    }

    /*
    // Get the trigger time from the f1 TDC
    vector<const DF1TDCHit*> tdchit;
    eventLoop->Get(tdchit);

    int tref = 0;
    for(unsigned int i=0;i<tdchit.size();i++)
    {
        if(tdchit[i]->rocid==51 && tdchit[i]->slot==17 && tdchit[i]->channel==8)
        {
            tref=tdchit[i]->time; // in clicks
            if (tref > 0) break;
            //       printf("tref %d %f\n",tdchit[i]->time,tref);
        }
    }
    //if (tref > 0)
    //{ // got reference signal
        // Next, loop over TDC hits, matching them to the
        // existing fADC hits where possible and updating
        // their time information. If no match is found, then
        // create a new hit with just the TDC info.
        vector<const DSCTDCDigiHit*> tdcdigihits;
        loop->Get(tdcdigihits);
        sort(tdcdigihits.begin(),tdcdigihits.end(),DTPOLHit_tdc_cmp);

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

            // Take care of rollover
            int tdiff = int(digihit->time) - int(tref);
            if (tdiff < 0) tdiff += rollover_count;
            else if (tdiff > rollover_count) tdiff -= rollover_count;

            // Apply calibration constants here
            double T = (double)digihit->time;

            //printf("T %d %f\n",digihit->time,0.0559*T);
            //tdc_scale=0.0559; // hard code correctd tdc conversion scale (need to put in CCDB)
            unsigned int id = digihit->sector - 1;
            T = tdc_scale * tdiff - tdc_time_offsets[id] + t_tdc_base;

            // cout << "T = " << T << endl;
            // jout << "T = " << T << endl;
            // printf("T = %f scale %f\n", T, tdc_scale);

            // Look for existing hits to see if there is a match
            //   or create new one if there is no match
	    // Require that the trigger corrected TDC time fall within 
	    //   a reasonable time window so that when a hit is associated with
	    //   a hit in the TDC and not the ADC it is a "decent" TDC hit
	    if (fabs(T) < HIT_TIME_WINDOW)
	      {
		//jout << " T cut = " << T << endl;
		DTPOLHit *hit = FindMatch(digihit->sector, T); 
		if (! hit) 
		  {
		    hit = new DTPOLHit;
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
		    // Correct for time walk
		    // The correction is the form t=t_tdc- C1 (A^C2 - A0^C2)
		    double A  = hit->dE;
		    double C1 = timewalk_parameters[id][1];
		    double C2 = timewalk_parameters[id][2];
		    double A0 = timewalk_parameters[id][3];
		    T -= C1*(pow(A,C2) - pow(A0,C2));	
		  }
		hit->t=T;
		//jout << " T cut TW Corr = " << T << endl;
		hit->AddAssociatedObject(digihit);
	      } // Hit time window cut
	      
        }
	//}


    // Apply calibration constants to convert pulse integrals to energy 
    // units
    for (unsigned int i=0;i<_data.size();i++)
      {
        _data[i]->dE*=a_scale * a_gains[_data[i]->sector-1];
      }
    */

    return NOERROR;
}

/*
//------------------
// FindMatch
//------------------
DTPOLHit* DTPOLHit_factory::FindMatch(int sector, double T)
{
    DTPOLHit *best_match=NULL;

    // Loop over existing hits (from fADC) and look for a match
    // in both the sector and the time.
    for(unsigned int i = 0; i < _data.size(); i++)
    {
        DTPOLHit *hit = _data[i];

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
*/

//------------------
// erun
//------------------
jerror_t DTPOLHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTPOLHit_factory::fini(void)
{
    return NOERROR;
}


//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DTPOLHit_factory::GetConstant(const vector<double> &the_table,
        const int in_sector) const{
    char str[256];

    if ( (in_sector < 0) || (in_sector >= DTPOLHit_factory::NSECTORS)){
      sprintf(str, "Bad sector # requested in DTPOLHit_factory::GetConstant()!"
	      " requested=%d , should be %ud", in_sector, DTPOLHit_factory::NSECTORS);
      cerr << str << endl;
      throw JException(str);
    }

    return the_table[in_sector];
}

const double DTPOLHit_factory::GetConstant(const vector<double> &the_table,
        const DSCDigiHit *in_digihit) const{
    char str[256];

    if ( (in_digihit->sector < 0) || (in_digihit->sector >= DTPOLHit_factory::NSECTORS)){
      sprintf(str, "Bad sector # requested in DTPOLHit_factory::GetConstant()!"
	      " requested=%d , should be %ud", 
	      in_digihit->sector, DTPOLHit_factory::NSECTORS);
      cerr << str << endl;
      throw JException(str);
    }

    return the_table[in_digihit->sector];
}

const double DTPOLHit_factory::GetConstant(const vector<double> &the_table,
					   const DTPOLHit *in_hit) const {
  char str[256];
  
  if ( (in_hit->sector < 0) || (in_hit->sector >= DTPOLHit_factory::NSECTORS)){
    sprintf(str, "Bad sector # requested in DTPOLHit_factory::GetConstant()!"
	    " requested=%d , should be %ud",
	    in_hit->sector, DTPOLHit_factory::NSECTORS);
    cerr << str << endl;
    throw JException(str);
  }
  
  return the_table[in_hit->sector];
}
/*
   const double DTPOLHit_factory::GetConstant(const vector<double> &the_table,
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
   sprintf(str, "Bad sector # requested in DTPOLHit_factory::GetConstant()!"
   " requested=%d , should be %ud",
   channel_info.sc.sector, MAX_SECTORS);
   cerr << str << endl;
   throw JException(str);
   }

   return the_table[channel_info.sc.sector];
   }
   */
