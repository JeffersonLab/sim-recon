// $Id$
//
//    File: DTAGHHit_factory.cc
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluex.phys.uconn.edu)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTAGHDigiHit.h"
#include "DTAGHTDCDigiHit.h"
#include "DTAGHGeometry.h"
#include "DTAGHHit_factory.h"
#include <DAQ/Df250PulseIntegral.h>
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
   // initialize calibration constants
   fadc_a_scale = 0;
   fadc_t_scale = 0;
   tdc_t_scale = 0;
   t_min = 0;

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
   /// set the base conversion scales
   fadc_a_scale    = 1.1;        // pixels per count
   fadc_t_scale    = 0.0625;     // ns per count
   tdc_t_scale     = 0.0600;     // ns per count
   t_min           = -100.;      // ns

   jout << "In DTAGHHit_factory, loading constants..." << std::endl;

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

      // Get pedestal, prefer associated event pedestal if it exists,
      // otherwise, use the average pedestal from CCDB
      double pedestal = fadc_pedestals[digihit->counter_id];
      vector<const Df250PulseIntegral*> PIvect;
      digihit->Get(PIvect);
      if (!PIvect.empty()) {
          const Df250PulseIntegral *PIobj = PIvect[0];
          pedestal = PIobj->pedestal;
      }

      DTAGHHit *hit = new DTAGHHit;
      int counter = digihit->counter_id;
      hit->counter_id = counter;
      double Elow = taghGeom.getElow(counter);
      double Ehigh = taghGeom.getEhigh(counter);
      hit->E = (Elow + Ehigh)/2;
      hit->t = 0;

      // throw away hits from bad or noisy counters
      int quality = counter_quality[counter];
      if (quality == k_counter_bad || quality == k_counter_noisy) 
          continue;

      // Apply calibration constants
      double A = digihit->pulse_integral;
      double T = digihit->pulse_time;
      A -= pedestal * digihit->nsamples_integral;
      hit->npe_fadc = A * fadc_a_scale * fadc_gains[counter];
      hit->time_fadc = T * fadc_t_scale - tdc_time_offsets[counter] + t_min;

      hit->AddAssociatedObject(digihit);
      _data.push_back(hit);
   }
      
   // Second, loop over TDC hits, matching them to the existing fADC hits
   // where possible and updating their time information. If no match is
   // found, then create a new hit with just the TDC info.
   vector<const DTAGHTDCDigiHit*> tdcdigihits;
   loop->Get(tdcdigihits);
   for (unsigned int i=0; i < tdcdigihits.size(); i++) {
      const DTAGHTDCDigiHit *digihit = tdcdigihits[i];

      // Apply calibration constants here
      int counter = digihit->counter_id;
      double T = (double)digihit->time;
      T = T * tdc_t_scale - tdc_time_offsets[counter] + t_min;

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
         hit->time_fadc = 0;
         hit->npe_fadc = 0;
         _data.push_back(hit);
      }      
      hit->t = T;

      // apply time-walk corrections?
      
      hit->AddAssociatedObject(digihit);
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
      int counter = (table[i])["counter"];
      result[counter] = (table[i])[column_name];
   }
   return true;
}
