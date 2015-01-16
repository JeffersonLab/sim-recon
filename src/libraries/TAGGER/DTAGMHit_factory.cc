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
#include <DAQ/Df250Config.h>
#include "DAQ/DF1TDCHit.h"

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
   t_base = 0;
   t_tdc_base=0;

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
   t_base           = 0.;      // ns

   jout << "In DTAGMHit_factory, loading constants..." << std::endl;

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
   if (eventLoop->GetCalib("/PHOTON_BEAM/microscope/base_time_offset",base_time_offset))
       jout << "Error loading /PHOTON_BEAM/microscope/base_time_offset !" << endl;
   if (base_time_offset.find("TAGM_BASE_TIME_OFFSET") != base_time_offset.end())
       t_base = base_time_offset["TAGM_BASE_TIME_OFFSET"];
   else
       jerr << "Unable to get TAGM_BASE_TIME_OFFSET from /PHOTON_BEAM/microscope/base_time_offset !" << endl;
   if (base_time_offset.find("TAGM_TDC_BASE_TIME_OFFSET") != base_time_offset.end())
       t_tdc_base = base_time_offset["TAGM_TDC_BASE_TIME_OFFSET"];
   else
       jerr << "Unable to get TAGM_TDC_BASE_TIME_OFFSET from /PHOTON_BEAM/microscope/base_time_offset !" << endl;

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
   jout << "TDC scale = " << tdc_t_scale << " ns/count" << endl;

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
      const Df250PulseIntegral* PIobj = NULL;
      const Df250Config *configObj = NULL;
      digihit->GetSingle(PIobj);
      PIobj->GetSingle(configObj);
      if ((PIobj != NULL) && (configObj != NULL)) {
	  // the measured pedestal must be scaled by the ratio of the number
	  // of samples used to calculate the pedestal and the actual pulse
	  pedestal = static_cast<double>(configObj->NSA_NSB) * PIobj->pedestal;
      }

      // throw away hits from bad or noisy fibers
      int quality = fiber_quality[digihit->row][digihit->column];
      if (quality == k_fiber_bad || quality == k_fiber_noisy)
            continue;

      // Skip events where fADC algorithm fails
      if (digihit->pulse_time == 0) continue;

      DTAGMHit *hit = new DTAGMHit;
      int row = digihit->row;
      int column = digihit->column;
      hit->row = row;
      hit->column = column;
      double Elow = tagmGeom.getElow(column);
      double Ehigh = tagmGeom.getEhigh(column);
      hit->E = (Elow + Ehigh)/2;
      hit->t = 0;

      // Apply calibration constants
      double A = digihit->pulse_integral;
      double T = digihit->pulse_time;
      A -= pedestal * digihit->nsamples_integral;
      hit->integral=A;
      hit->npix_fadc = A * fadc_a_scale * fadc_gains[row][column];
      hit->time_fadc = T * fadc_t_scale - tdc_time_offsets[row][column] + t_base;
      hit->t=hit->time_fadc;
      hit->has_TDC=false;
      hit->has_fADC=true;
      
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
     vector<const DTAGMTDCDigiHit*> tdcdigihits;
     loop->Get(tdcdigihits);
     for (unsigned int i=0; i < tdcdigihits.size(); i++) {
       const DTAGMTDCDigiHit *digihit = tdcdigihits[i];
       
       // Take care of rollover
       int tdiff = int(digihit->time) - int(tref);
       if (tdiff < 0) tdiff += rollover_count;
       else if (tdiff > rollover_count) tdiff -= rollover_count;
       
       // Apply calibration constants here
       int row = digihit->row;
       int column = digihit->column;
       double T = tdiff * tdc_t_scale - tdc_time_offsets[row][column] + t_tdc_base;
       
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
	 hit->has_fADC=false;
         _data.push_back(hit);
       }      
       hit->time_tdc=T;
       hit->has_TDC=true;

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
