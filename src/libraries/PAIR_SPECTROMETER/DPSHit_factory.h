// $Id$
//
//    File: DPSHit_factory.h
// Created: Wed Oct 15 16:45:01 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//


#ifndef _DPSHit_factory_
#define _DPSHit_factory_

#include <JANA/JFactory.h>
#include "DPSHit.h"
#include "DPSDigiHit.h"
#include "DPSGeometry.h"

#include <vector>
#include <utility>

typedef vector< vector<double> > ps_digi_constants_t;

class DPSHit_factory:public jana::JFactory<DPSHit>{
	public:
		DPSHit_factory(){};
		~DPSHit_factory(){};


                // overall scale factors
                double a_scale;
                double t_scale;
                double t_base;

                // calibration constants stored by channel
                ps_digi_constants_t  adc_gains;
                ps_digi_constants_t  adc_pedestals;
                ps_digi_constants_t  adc_time_offsets;

                const double GetConstant( const ps_digi_constants_t &the_table,
                                          const DPSGeometry::Arm in_arm, const int in_column,
					  const DPSGeometry &psGeom ) const; 
                const double GetConstant( const ps_digi_constants_t &the_table,
                                          const DPSDigiHit *the_digihit, const DPSGeometry &psGeom ) const;
                const double GetConstant( const ps_digi_constants_t &the_table,
                                          const DPSHit *the_hit, const DPSGeometry &psGeom ) const;


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

                void FillCalibTable(ps_digi_constants_t &table, string table_name,
                                    const DPSGeometry &tofGeom);

};

#endif // _DPSHit_factory_

