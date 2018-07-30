// $Id$
//
//    File: DCDCHit_factory_Calib.h
// Created: Fri Mar  9 16:57:04 EST 2018
// Creator: B. Zihlmann derived/copyed from DCDCHit_factory.h
//

#ifndef _DCDCHit_factory_Calib_
#define _DCDCHit_factory_Calib_

#include <vector>
using namespace std;

#include <JANA/JFactory.h>
#include <HDGEOMETRY/DGeometry.h>
#include <TTAB/DTranslationTable.h>
#include <DAQ/Df125PulseIntegral.h>
#include <DAQ/Df125Config.h>
#include <DAQ/Df125CDCPulse.h>

#include "DCDCHit.h"
#include "DCDCDigiHit.h"
#include "DCDCWire.h"

// store constants indexed by ring/straw number
typedef  vector< vector<double> >  cdc_digi_constants_t;

class DCDCHit_factory_Calib:public jana::JFactory<DCDCHit>{
 public:
  DCDCHit_factory_Calib(){};
  ~DCDCHit_factory_Calib(){};
  const char* Tag(void){return "Calib";}

  // overall scale factors.
  double a_scale, amp_a_scale;
  double t_scale;
  double t_base;
  
  // calibration constant tables
  cdc_digi_constants_t gains;
  cdc_digi_constants_t pedestals;
  cdc_digi_constants_t time_offsets;
  
  const double GetConstant(const cdc_digi_constants_t &the_table,
			   const int in_ring, const int in_straw) const;
  const double GetConstant(const cdc_digi_constants_t &the_table,
			   const DCDCDigiHit *the_digihit) const;
  const double GetConstant(const cdc_digi_constants_t &the_table,
			   const DCDCHit *the_hit) const;
  //const double GetConstant(const cdc_digi_constants_t &the_table,
  //			 const DTranslationTable *ttab,
  //			 const int in_rocid, const int in_slot, const int in_channel) const;
  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  void CalcNstraws(jana::JEventLoop *eventLoop, int32_t runnumber, vector<unsigned int> &Nstraws);
  void FillCalibTable(vector< vector<double> > &table, vector<double> &raw_table, 
		      vector<unsigned int> &Nstraws);
  
  // Geometry information
  unsigned int maxChannels;
  unsigned int Nrings; // number of rings (layers)
  vector<unsigned int> Nstraws; // number of straws for each ring
  
};

#endif // _DCDCHit_factory_Calib_

