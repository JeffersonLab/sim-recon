// $Id$
//
//    File: DTAGMHit_factory.cc
// Created: Sat Aug  2 12:23:43 EDT 2014
// Creator: jonesrt (on Linux gluex.phys.uconn.edu)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DTAGMDigiHit.h"
#include "DTAGMTDCDigiHit.h"
#include "DTAGMGeometry.h"
#include "DTAGMHit_factory.h"
#include <DAQ/Df250PulseIntegral.h>
using namespace jana;

const int DTAGMHit_factory::k_fiber_good;
const int DTAGMHit_factory::k_fiber_bad;
const int DTAGMHit_factory::k_fiber_noisy;

#define DELTA_T_ADC_TDC_MATCH_NS 10.0

//------------------
// init
//------------------
jerror_t DTAGMHit_factory::init(void)
{
   // initialize calibration constants
   fadc_a_scale = 0;
   fadc_t_scale = 0;
   tdc_t_scale = 0;
   t_min = 0;

   // calibration constants stored in row, column format
   for (int row = 0; row <= TAGM_MAX_ROW; ++row) {
      for (int col = 0; col <= TAGM_MAX_COLUMN; ++col) {
         fadc_gains[row][col] = 0;
         fadc_pedestals[row][col] = 0;
         fadc_time_offsets[row][col] = 0;
         tdc_time_offsets[row][col] = 0;
         fiber_quality[row][col] = 0;
      }
   }
    
   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTAGMHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
   /// set the base conversion scales
   fadc_a_scale    = 1.1;        // pixels per count
   fadc_t_scale    = 0.0625;     // ns per count
   tdc_t_scale     = 0.0600;     // ns per count
   t_min           = -100.;      // ns

   jout << "In DTAGMHit_factory, loading constants..." << std::endl;

   if (load_ccdb_constants("fadc_gains", "gain", fadc_gains) &&
       load_ccdb_constants("fadc_pedestals", "pedestal", fadc_pedestals) &&
       load_ccdb_constants("fadc_time_offsets", "offset", fadc_time_offsets) &&
       load_ccdb_constants("tdc_time_offsets", "offset", tdc_time_offsets) &&
       load_ccdb_constants("fiber_quality", "code", fiber_quality) )
   {
      return NOERROR;
   }
   return UNRECOVERABLE_ERROR;
}

//------------------
// evnt
//------------------
jerror_t DTAGMHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
   /// Generate DTAGMHit object for each DTAGMDigiHit object.
   /// This is where the first set of calibration constants
   /// is applied to convert from digitzed units into natural
   /// units.
   ///
   /// Note that this code does NOT get called for simulated
   /// data in HDDM format. The HDDM event source will copy
   /// the precalibrated values directly into the _data vector.

   // extract the TAGM geometry
   vector<const DTAGMGeometry*> tagmGeomVect;
   eventLoop->Get( tagmGeomVect );
   if (tagmGeomVect.size() < 1)
      return OBJECT_NOT_AVAILABLE;
   const DTAGMGeometry& tagmGeom = *(tagmGeomVect[0]);

   // First loop over all TAGMDigiHits and make DTAGMHits out of them
   vector<const DTAGMDigiHit*> digihits;
   loop->Get(digihits);
   for (unsigned int i=0; i < digihits.size(); i++) {
      const DTAGMDigiHit *digihit = digihits[i];

      // Get pedestal, prefer associated event pedestal if it exists,
      // otherwise, use the average pedestal from CCDB
      double pedestal = fadc_pedestals[digihit->row][digihit->column];
      vector<const Df250PulseIntegral*> PIvect;
      digihit->Get(PIvect);
      if (!PIvect.empty()) {
          const Df250PulseIntegral *PIobj = PIvect[0];
          pedestal = PIobj->pedestal;
      }

      DTAGMHit *hit = new DTAGMHit;
      int row = digihit->row;
      int column = digihit->column;
      hit->row = row;
      hit->column = column;
      double Elow = tagmGeom.getElow(column);
      double Ehigh = tagmGeom.getEhigh(column);
      hit->E = (Elow + Ehigh)/2;
      hit->t = 0;

      // throw away hits from bad or noisy fibers
      int quality = fiber_quality[row][column];
      if (quality == k_fiber_bad || quality == k_fiber_noisy) 
          continue;

      // Apply calibration constants
      double A = digihit->pulse_integral;
      double T = digihit->pulse_time;
      A -= pedestal * digihit->nsamples_integral;
      hit->npix_fadc = A * fadc_a_scale * fadc_gains[row][column];
      hit->time_fadc = T * fadc_t_scale - tdc_time_offsets[row][column] + t_min;

      hit->AddAssociatedObject(digihit);
      _data.push_back(hit);
   }
      
   // Second, loop over TDC hits, matching them to the existing fADC hits
   // where possible and updating their time information. If no match is
   // found, then create a new hit with just the TDC info.
   vector<const DTAGMTDCDigiHit*> tdcdigihits;
   loop->Get(tdcdigihits);
   for (unsigned int i=0; i < tdcdigihits.size(); i++) {
      const DTAGMTDCDigiHit *digihit = tdcdigihits[i];

      // Apply calibration constants here
      int row = digihit->row;
      int column = digihit->column;
      double T = (double)digihit->time;
      T = T * tdc_t_scale - tdc_time_offsets[row][column] + t_min;

      // Look for existing hits to see if there is a match
      // or create new one if there is no match
      DTAGMHit *hit = 0;
      for (unsigned int j=0; j < _data.size(); ++j) {
         if (_data[j]->row == row && _data[j]->column == column &&
             fabs(T - _data[j]->time_fadc) < DELTA_T_ADC_TDC_MATCH_NS)
         {
            hit = _data[j];
         }
      }
      if (hit == 0) {
         hit = new DTAGMHit;
         hit->row = row;
         hit->column = column;
         double Elow = tagmGeom.getElow(column);
         double Ehigh = tagmGeom.getEhigh(column);
         hit->E = (Elow + Ehigh)/2;
         hit->time_fadc = 0;
         hit->npix_fadc = 0;
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
jerror_t DTAGMHit_factory::erun(void)
{
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DTAGMHit_factory::fini(void)
{
   return NOERROR;
}

//---------------------
// load_ccdb_constants
//---------------------
bool DTAGMHit_factory::load_ccdb_constants(
                       std::string table_name,
                       std::string column_name,
                       double result[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1])
{
   std::vector< std::map<std::string, double> > table;
   std::string ccdb_key = "/PHOTON_BEAM/microscope/" + table_name;
   if (eventLoop->GetCalib(ccdb_key, table))
   {
       jout << "Error loading " << ccdb_key << " from ccdb!" << std::endl;
       return false;
   }
   for (unsigned int i=0; i < table.size(); ++i) {
      int row = (table[i])["row"];
      int col = (table[i])["col"];
      result[row][col] = (table[i])[column_name];
   }
   return true;
}
