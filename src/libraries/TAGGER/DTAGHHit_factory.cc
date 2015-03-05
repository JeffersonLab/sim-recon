// $Id$
//
//    File: DTAGHHit_factory.cc
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluex.phys.uconn.edu)
//


#include <iostream>
#include <iomanip>
#include <limits>
using namespace std;

#include "DTAGHDigiHit.h"
#include "DTAGHTDCDigiHit.h"
#include "DTAGHGeometry.h"
#include "DTAGHHit_factory.h"
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250PulsePedestal.h>
#include "DAQ/DF1TDCHit.h"

using namespace jana;

const int DTAGHHit_factory::k_counter_dead;
const int DTAGHHit_factory::k_counter_good;
const int DTAGHHit_factory::k_counter_bad;
const int DTAGHHit_factory::k_counter_noisy;

#define DELTA_T_ADC_TDC_MATCH_NS 10.0

//------------------
// init
//------------------
jerror_t DTAGHHit_factory::init(void)
{
  ADC_THRESHOLD = 1000.0; // ADC integral counts
  gPARMS->SetDefaultParameter("TAGHHit:ADC_THRESHOLD",ADC_THRESHOLD,
			      "pedestal-subtracted pulse integral threshold");

   // initialize calibration constants
   fadc_a_scale = 0;
   fadc_t_scale = 0;
   tdc_t_scale = 0;
   t_base = 0;
   t_tdc_base=0;

   // calibration constants stored by counter index
   for (int counter = 0; counter <= TAGH_MAX_COUNTER; ++counter) {
      fadc_gains[counter] = 0;
      fadc_pedestals[counter] = 0;
      fadc_time_offsets[counter] = 0;
      tdc_time_offsets[counter] = 0;
      counter_quality[counter] = 0;
   }
    
   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTAGHHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
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

   /// set the base conversion scales
   fadc_a_scale    = 1.1;        // pixels per count
   fadc_t_scale    = 0.0625;     // ns per count
   tdc_t_scale     = 0.0600;     // ns per count
   t_base           = 0.;      // ns

   if(print_messages) jout << "In DTAGHHit_factory, loading constants..." << std::endl;

   // F1TDC tframe(ns) and rollover count
   map<string,int> tdc_parms;
   double tframe=-1.;
   rollover_count=0;
   if (eventLoop->GetCalib("/F1TDC/rollover",tdc_parms))
     jout << "Error loading /F1TDC/rollover !" <<endl;
   if (tdc_parms.find("tframe")!=tdc_parms.end())
     tframe=double(tdc_parms["tframe"]);
   else 
     jerr << "Unable to get tframe from /F1TDC/rollover !" <<endl;
   if (tdc_parms.find("count")!=tdc_parms.end())
     rollover_count=tdc_parms["count"];
   else 
     jerr << "Unable to get rollover count from /F1TDC/rollover !" <<endl;


   // load base time offset
   map<string,double> base_time_offset;
   if (eventLoop->GetCalib("/PHOTON_BEAM/hodoscope/base_time_offset",base_time_offset))
       jout << "Error loading /PHOTON_BEAM/hodoscope/base_time_offset !" << endl;
   if (base_time_offset.find("TAGH_BASE_TIME_OFFSET") != base_time_offset.end())
       t_base = base_time_offset["TAGH_BASE_TIME_OFFSET"];
   else
       jerr << "Unable to get TAGH_BASE_TIME_OFFSET from /PHOTON_BEAM/hodoscope/base_time_offset !" << endl;

    if (base_time_offset.find("TAGH_TDC_BASE_TIME_OFFSET") != base_time_offset.end())
       t_tdc_base = base_time_offset["TAGH_TDC_BASE_TIME_OFFSET"];
   else
       jerr << "Unable to get TAGH_TDC_BASE_TIME_OFFSET from /PHOTON_BEAM/hodoscope/base_time_offset !" << endl;

   // tdc_t_scale
   // By default we will use the F1TDC entries in the ccdb to find tdc_scale
   if (tframe>0. && rollover_count>0){ 
     tdc_t_scale=tframe/double(rollover_count);
     //_DBG_ << tdc_scale << endl;
   }
   else {
     jerr << "Unable to get TDC_SCALE from database !" 
	  << endl;
   }
   if(print_messages) jout << "TDC scale = " << tdc_t_scale << " ns/count" << endl;

   if (load_ccdb_constants("fadc_gains", "gain", fadc_gains) &&
       load_ccdb_constants("fadc_pedestals", "pedestal", fadc_pedestals) &&
       load_ccdb_constants("fadc_time_offsets", "offset", fadc_time_offsets) &&
       load_ccdb_constants("tdc_time_offsets", "offset", tdc_time_offsets) &&
       load_ccdb_constants("counter_quality", "code", counter_quality) )
   {
      return NOERROR;
   }
   return UNRECOVERABLE_ERROR;
}

//------------------
// evnt
//------------------
jerror_t DTAGHHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
   /// Generate DTAGHHit object for each DTAGHDigiHit object.
   /// This is where the first set of calibration constants
   /// is applied to convert from digitzed units into natural
   /// units.
   ///
   /// Note that this code does NOT get called for simulated
   /// data in HDDM format. The HDDM event source will copy
   /// the precalibrated values directly into the _data vector.

   // extract the TAGH geometry
   vector<const DTAGHGeometry*> taghGeomVect;
   eventLoop->Get( taghGeomVect );
   if (taghGeomVect.size() < 1)
      return OBJECT_NOT_AVAILABLE;
   const DTAGHGeometry& taghGeom = *(taghGeomVect[0]);

   // First loop over all TAGHDigiHits and make DTAGHHits out of them
   vector<const DTAGHDigiHit*> digihits;
   loop->Get(digihits);
   for (unsigned int i=0; i < digihits.size(); i++) {
      const DTAGHDigiHit *digihit = digihits[i];
      int counter = digihit->counter_id;

      // throw away hits from bad or noisy counters
      int quality = counter_quality[counter];
      if (quality == k_counter_bad || quality == k_counter_noisy) 
          continue;

      // Throw away hits where the fADC timing algorithm failed
      //if (digihit->pulse_time == 0) continue;
      // The following condition signals an error state in the flash algorithm
      // Do not make hits out of these
      const Df250PulsePedestal* PPobj = NULL;
      digihit->GetSingle(PPobj);
      if (PPobj != NULL){
          if (PPobj->pedestal == 0 || PPobj->pulse_peak == 0) continue;
      }

      // Get pedestal, prefer associated event pedestal if it exists,
      // otherwise, use the average pedestal from CCDB
      double pedestal = fadc_pedestals[counter];
      const Df250PulseIntegral* PIobj = NULL;
      digihit->GetSingle(PIobj);
      if (PIobj != NULL) {
	  // the measured pedestal must be scaled by the ratio of the number
	  // of samples used to calculate the integral and the pedestal          
	  // Changed to conform to D. Lawrence changes Dec. 4 2014
          double ped_sum = (double)PIobj->pedestal;
          double nsamples_integral = (double)PIobj->nsamples_integral;
          double nsamples_pedestal = (double)PIobj->nsamples_pedestal;
          pedestal          = ped_sum * nsamples_integral/nsamples_pedestal;
      }

      // Subtract pedestal from pulse integral
      double A = digihit->pulse_integral;
      A -= pedestal;
      // Throw away hits with small pedestal-subtracted integrals
      if (A < ADC_THRESHOLD) continue;

      DTAGHHit *hit = new DTAGHHit;
      hit->counter_id = counter;
      double Elow = taghGeom.getElow(counter);
      double Ehigh = taghGeom.getEhigh(counter);
      hit->E = (Elow + Ehigh)/2;

      // Apply calibration constants
      double T = digihit->pulse_time;
      hit->integral = A;
      hit->npe_fadc = A * fadc_a_scale * fadc_gains[counter];
      hit->time_fadc = T * fadc_t_scale - fadc_time_offsets[counter] + t_base;
      hit->time_tdc = numeric_limits<double>::quiet_NaN();
      hit->t = hit->time_fadc;
      hit->has_fADC = true;
      hit->has_TDC = false;

      hit->AddAssociatedObject(digihit);
      _data.push_back(hit);
   }

   // Get the trigger time from the f1 TDC
   vector<const DF1TDCHit*> tdchit;
   eventLoop->Get(tdchit);

   int tref = 0;
   for(unsigned int i=0;i<tdchit.size();i++)
   {
       if(tdchit[i]->rocid==51 && tdchit[i]->slot==17 && tdchit[i]->channel==8)
       {
           tref=tdchit[i]->time; // in clicks
           break;
           //       printf("tref %d %f\n",tdchit[i]->time,tref);
       }
   }
   if (tref > 0){
       // Next, loop over TDC hits, matching them to the existing fADC hits
       // where possible and updating their time information. If no match is
       // found, then create a new hit with just the TDC info.
       vector<const DTAGHTDCDigiHit*> tdcdigihits;
       loop->Get(tdcdigihits);
       for (unsigned int i=0; i < tdcdigihits.size(); i++) {
           const DTAGHTDCDigiHit *digihit = tdcdigihits[i];

           // Take care of rollover
           int tdiff = int(digihit->time) - int(tref);
           if (tdiff < 0) tdiff += rollover_count;
           else if (tdiff > rollover_count) tdiff -= rollover_count;

           // Apply calibration constants here
           int counter = digihit->counter_id;
           double T = tdiff* tdc_t_scale - tdc_time_offsets[counter] + t_tdc_base;

           // Look for existing hits to see if there is a match
           // or create new one if there is no match
           DTAGHHit *hit = 0;
           for (unsigned int j=0; j < _data.size(); ++j) {
               if (_data[j]->counter_id == counter &&
                       fabs(T - _data[j]->time_fadc) < DELTA_T_ADC_TDC_MATCH_NS)
               {
                   hit = _data[j];
               }
           }
           if (hit == 0) {
               hit = new DTAGHHit;
               hit->counter_id = counter;
               double Elow = taghGeom.getElow(counter);
               double Ehigh = taghGeom.getEhigh(counter);
               hit->E = (Elow + Ehigh)/2;
               hit->time_fadc = numeric_limits<double>::quiet_NaN();
               hit->integral = numeric_limits<double>::quiet_NaN();
               hit->npe_fadc = numeric_limits<double>::quiet_NaN();
               hit->has_fADC = false;
               _data.push_back(hit);
           }      
           hit->time_tdc = T;
           hit->has_TDC = true;
           hit->t = T;

           // apply time-walk corrections?

           hit->AddAssociatedObject(digihit);
       }
   }

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTAGHHit_factory::erun(void)
{
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTAGHHit_factory::fini(void)
{
    return NOERROR;
}

//---------------------
// load_ccdb_constants
//---------------------
bool DTAGHHit_factory::load_ccdb_constants(
        std::string table_name,
        std::string column_name,
        double result[TAGH_MAX_COUNTER+1])
{
    std::vector< std::map<std::string, double> > table;
    std::string ccdb_key = "/PHOTON_BEAM/hodoscope/" + table_name;
    if (eventLoop->GetCalib(ccdb_key, table))
    {
        jout << "Error loading " << ccdb_key << " from ccdb!" << std::endl;
        return false;
    }
    for (unsigned int i=0; i < table.size(); ++i) {
        int counter = (table[i])["id"];
        result[counter] = (table[i])[column_name];
    }
    return true;
}
