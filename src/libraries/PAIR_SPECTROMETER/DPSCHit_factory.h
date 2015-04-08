// $Id$
//
//    File: DPSCHit_factory.h
// Created: Wed Oct 15 16:45:33 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSCHit_factory_
#define _DPSCHit_factory_

#include <vector>
#include <utility>

#include <JANA/JFactory.h>
#include "TTAB/DTTabUtilities.h"

#include "DPSCHit.h"
#include "DPSCDigiHit.h"
#include "DPSCTDCDigiHit.h"
#include "DPSGeometry.h"

typedef vector< pair<double,double> > psc_digi_constants_t;

class DPSCHit_factory:public jana::JFactory<DPSCHit>{
 public:
  DPSCHit_factory(){};
  ~DPSCHit_factory(){};

  // config. parameters
  double DELTA_T_ADC_TDC_MAX;  
  double ADC_THRESHOLD;

  // overall scale factors
  double a_scale;
  double t_scale;
  double t_base;
  double t_tdc_base;

  // calibration constants stored by channel
  psc_digi_constants_t  adc_gains;
  psc_digi_constants_t  adc_pedestals;
  psc_digi_constants_t  adc_time_offsets;
  psc_digi_constants_t  tdc_time_offsets;

  const DPSGeometry::Arm GetArm(const int counter_id,const int num_counters_per_arm) const;
  const int GetModule(const int counter_id,const int num_counters_per_arm) const;
  DPSCHit* FindMatch(DPSGeometry::Arm arm, int module, double T);

  const double GetConstant( const psc_digi_constants_t &the_table,
			    const DPSGeometry::Arm in_arm, const int in_module,
			    const DPSGeometry &psGeom ) const; 
  const double GetConstant( const psc_digi_constants_t &the_table,
			    const DPSCDigiHit *the_digihit, const DPSGeometry &psGeom ) const;
  const double GetConstant( const psc_digi_constants_t &the_table,
			    const DPSCTDCDigiHit *the_hit, const DPSGeometry &psGeom ) const;
  const double GetConstant( const psc_digi_constants_t &the_table,
			    const DPSCHit *the_hit, const DPSGeometry &psGeom ) const;

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

  void FillCalibTable(psc_digi_constants_t &table, vector<double> &raw_table,
		      const DPSGeometry &tofGeom);

};

#endif // _DPSCHit_factory_

