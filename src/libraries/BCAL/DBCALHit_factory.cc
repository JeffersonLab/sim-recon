// $Id$
//
//    File: DBCALHit_factory.cc
// Created: Tue Aug  6 09:26:13 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include <BCAL/DBCALDigiHit.h>
#include "DBCALGeometry.h"
#include "DBCALHit_factory.h"
#include <DAQ/Df250PulseIntegral.h>
using namespace jana;


//------------------
// init
//------------------
jerror_t DBCALHit_factory::init(void)
{
   /// set the base conversion scales
   a_scale    = 0.1;   // to get units of MeV
   //  Crude calibration 
   //    A minimally ionising particle deposits and integral of 230 ADC counts per cell, 
   //    which corresponds to approximately 22 MeV.  Thus, the factor is 0.1 to get MeV
   //a_pedestal = 10000;  // default pedestal of 100 ADC units over 100 samples 
   t_scale    = 0.0625;   // There are 62.5 ps/count from the fADC
   t_min      = -100.;

   return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DBCALHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
   /// Read in calibration constants
   vector<double> raw_gains;
   vector<double> raw_pedestals;
   vector<double> raw_time_offsets;

   jout << "In DBCALHit_factory, loading constants..." << endl;
   
   // load scale factors
   map<string,double> scale_factors;
   if (eventLoop->GetCalib("/BCAL/digi_scales", scale_factors))
       jout << "Error loading /BCAL/digi_scales !" << endl;
   if (scale_factors.find("BCAL_ADC_ASCALE") != scale_factors.end())
       a_scale = scale_factors["BCAL_ADC_ASCALE"];
   else
       jerr << "Unable to get BCAL_ADC_ASCALE from /BCAL/digi_scales !" << endl;
   if (scale_factors.find("BCAL_ADC_TSCALE") != scale_factors.end())
       t_scale = scale_factors["BCAL_ADC_TSCALE"];
   else
       jerr << "Unable to get BCAL_ADC_TSCALE from /BCAL/digi_scales !" << endl;

   if (eventLoop->GetCalib("/BCAL/ADC_gains", raw_gains))
       jout << "Error loading /BCAL/ADC_gains !" << endl;
   if (eventLoop->GetCalib("/BCAL/ADC_pedestals", raw_pedestals))
       jout << "Error loading /BCAL/ADC_pedestals !" << endl;
   if (eventLoop->GetCalib("/BCAL/ADC_timing_offsets", raw_time_offsets))
       jout << "Error loading /BCAL/ADC_timing_offsets !" << endl;

   FillCalibTable(gains, raw_gains);
   FillCalibTable(pedestals, raw_pedestals);
   FillCalibTable(time_offsets, raw_time_offsets);

   return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DBCALHit_factory::evnt(JEventLoop *loop, int eventnumber)
{
   /// Generate DBCALHit object for each DBCALDigiHit object.
   /// This is where the first set of calibration constants
   /// is applied to convert from digitzed units into natural
   /// units.
   ///
   /// Note that this code does NOT get called for simulated
   /// data in HDDM format. The HDDM event source will copy
   /// the precalibrated values directly into the _data vector.

   vector<const DBCALDigiHit*> digihits;
   loop->Get(digihits);
   for(unsigned int i=0; i<digihits.size(); i++){
      const DBCALDigiHit *digihit = digihits[i];
      
      // Get pedestal.  Prefer associated event pedestal if it exists.
      // Otherwise, use the average pedestal from CCDB
      double pedestalpersample = GetConstant(pedestals,digihit);
      if (digihit->nsamples_pedestal > 0)
         pedestalpersample = digihit->pedestal / digihit->nsamples_pedestal;
      else
         pedestalpersample = digihit->pedestal;

      DBCALHit *hit = new DBCALHit;
      hit->module = digihit->module;
      hit->layer  = digihit->layer;
      hit->sector = digihit->sector;
      hit->end    = digihit->end;
      
      // Apply calibration constants here
      double A = (double)digihit->pulse_integral;
      double T = (double)digihit->pulse_time;

      double gain = GetConstant(gains,digihit);
      double E = 0;
      if (A > 0) {
         E = a_scale * gain * (A - (pedestalpersample * digihit->nsamples_integral));
      }

      hit->E = E/1000.;  // all calorimeter hits have E in GeV
      hit->t = t_scale * (T - GetConstant(time_offsets,digihit)) + t_min;
      
      hit->AddAssociatedObject(digihit);
      
      _data.push_back(hit);
   }

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBCALHit_factory::erun(void)
{
   return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DBCALHit_factory::fini(void)
{
   return NOERROR;
}


//------------------
// FillCalibTable
//------------------
void DBCALHit_factory::FillCalibTable( map<int,cell_calib_t> &table, 
                   const vector<double> &raw_table) 
{
    char str[256];
    int channel = 0;

    // reset the table before filling it
    table.clear();

    for (int module=1; module<=BCAL_NUM_MODULES; module++) {
       for (int layer=1; layer<=BCAL_NUM_LAYERS; layer++) {
          for (int sector=1; sector<=BCAL_NUM_SECTORS; sector++) {
             if ((channel > BCAL_MAX_CHANNELS) || (channel+1 > BCAL_MAX_CHANNELS)) {  // sanity check
                sprintf(str, "Too many channels for BCAL table!"
                             " channel=%d (should be %d)", 
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
    if (channel != BCAL_MAX_CHANNELS) { 
       sprintf(str, "Not enough channels for BCAL table!"
                    " channel=%d (should be %d)", 
               channel, BCAL_MAX_CHANNELS);
       cerr << str << endl;
       throw JException(str);
    }
}

//------------------------------------
// GetConstant
//   Allow a few different interfaces
//------------------------------------
const double DBCALHit_factory::GetConstant(const bcal_digi_constants_t &the_table, 
                                           const int in_module, 
                                           const int in_layer,
                                           const int in_sector,
                                           const int in_end) const
{
   char str[256];
   
   if ( (in_module <= 0) || (in_module > BCAL_NUM_MODULES)) {
      sprintf(str, "Bad module # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", in_module, BCAL_NUM_MODULES);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_layer <= 0) || (in_layer > BCAL_NUM_LAYERS)) {
      sprintf(str, "Bad layer # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", in_layer, BCAL_NUM_LAYERS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_sector <= 0) || (in_sector > BCAL_NUM_SECTORS)) {
      sprintf(str, "Bad sector # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", in_sector, BCAL_NUM_SECTORS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_end != DBCALGeometry::kUpstream) && (in_end != DBCALGeometry::kDownstream) ) {
      sprintf(str, "Bad end # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 0-1", in_end);
      cerr << str << endl;
      throw JException(str);
   }

   const int the_cell = DBCALGeometry::cellId( in_module, in_layer, in_sector);
   
   if (in_end == DBCALGeometry::kUpstream) {
       // handle the upstream end
       return the_table.at(the_cell).first;
   } else {
       // handle the downstream end
       return the_table.at(the_cell).second;
   }

}

const double DBCALHit_factory::GetConstant( const bcal_digi_constants_t &the_table, 
                   const DBCALHit *in_hit) const
{
   char str[256];
   
   if ( (in_hit->module <= 0) || (in_hit->module > BCAL_NUM_MODULES)) {
      sprintf(str, "Bad module # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", in_hit->module, BCAL_NUM_MODULES);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_hit->layer <= 0) || (in_hit->layer > BCAL_NUM_LAYERS)) {
      sprintf(str, "Bad layer # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", in_hit->layer, BCAL_NUM_LAYERS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_hit->sector <= 0) || (in_hit->sector > BCAL_NUM_SECTORS)) {
      sprintf(str, "Bad sector # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", in_hit->sector, BCAL_NUM_SECTORS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_hit->end != DBCALGeometry::kUpstream) && 
        (in_hit->end != DBCALGeometry::kDownstream) )
   {
      sprintf(str, "Bad end # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 0-1", in_hit->end);
      cerr << str << endl;
      throw JException(str);
   }

   const int the_cell = DBCALGeometry::cellId(in_hit->module,
                                              in_hit->layer,
                                              in_hit->sector);
   
   if (in_hit->end == DBCALGeometry::kUpstream) {
       // handle the upstream end
       return the_table.at(the_cell).first;
       //return the_table[the_cell].first;
   } else {
       // handle the downstream end
       return the_table.at(the_cell).second;
       //return the_table[the_cell].second;
   }

}

const double DBCALHit_factory::GetConstant(const bcal_digi_constants_t &the_table, 
                                           const DBCALDigiHit *in_digihit) const
{
   char str[256];
   
   if ( (in_digihit->module <= 0) || (in_digihit->module > BCAL_NUM_MODULES)) {
      sprintf(str, "Bad module # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", 
              in_digihit->module, BCAL_NUM_MODULES);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_digihit->layer <= 0) || (in_digihit->layer > BCAL_NUM_LAYERS)) {
      sprintf(str, "Bad layer # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d",
              in_digihit->layer, BCAL_NUM_LAYERS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_digihit->sector <= 0) || (in_digihit->sector > BCAL_NUM_SECTORS)) {
      sprintf(str, "Bad sector # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", 
              in_digihit->sector, BCAL_NUM_SECTORS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (in_digihit->end != DBCALGeometry::kUpstream) &&
        (in_digihit->end != DBCALGeometry::kDownstream) )
   {
      sprintf(str, "Bad end # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 0-1", in_digihit->end);
      cerr << str << endl;
      throw JException(str);
   }

   const int the_cell = DBCALGeometry::cellId(in_digihit->module,
                                              in_digihit->layer,
                                              in_digihit->sector);
   
   if (in_digihit->end == DBCALGeometry::kUpstream) {
       // handle the upstream end
       return the_table.at(the_cell).first;
   } else {
       // handle the downstream end
       return the_table.at(the_cell).second;
   }

}
/*
const double DBCALHit_factory::GetConstant(const bcal_digi_constants_t &the_table,
                                           const DTranslationTable *ttab,
                                           const int in_rocid,
                                           const int in_slot,
                                           const int in_channel) const
{
   char str[256];
   
   DTranslationTable::csc_t daq_index = { in_rocid, in_slot, in_channel };
   DTranslationTable::DChannelInfo channel_info = ttab->GetDetectorIndex(daq_index);
   
   if ( (channel_info.bcal.module <= 0) 
       || (channel_info.bcal.module > static_cast<unsigned int>(BCAL_NUM_MODULES)))
   {
      sprintf(str, "Bad module # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", channel_info.bcal.module, BCAL_NUM_MODULES);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (channel_info.bcal.layer <= 0) 
       || (channel_info.bcal.layer > static_cast<unsigned int>(BCAL_NUM_LAYERS)))
   {
      sprintf(str, "Bad layer # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", channel_info.bcal.layer, BCAL_NUM_LAYERS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (channel_info.bcal.sector <= 0) 
       || (channel_info.bcal.sector > static_cast<unsigned int>(BCAL_NUM_SECTORS)))
   {
      sprintf(str, "Bad sector # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 1-%d", channel_info.bcal.sector, BCAL_NUM_SECTORS);
      cerr << str << endl;
      throw JException(str);
   }
   if ( (channel_info.bcal.end != DBCALGeometry::kUpstream) 
       && (channel_info.bcal.end != DBCALGeometry::kDownstream) )
   {
      sprintf(str, "Bad end # requested in DBCALHit_factory::GetConstant()!"
                   " requested=%d , should be 0-1", channel_info.bcal.end);
      cerr << str << endl;
      throw JException(str);
   }

   int the_cell = DBCALGeometry::cellId(channel_info.bcal.module,
                                        channel_info.bcal.layer,
                                        channel_info.bcal.sector);
   
   if (channel_info.bcal.end == DBCALGeometry::kUpstream) {
       // handle the upstream end
       return the_table.at(the_cell).first;
   } else {
       // handle the downstream end
       return the_table.at(the_cell).second;
   }
}
*/
