// $Id$
//
//    File: DSCHit_factory.h
// Created: Tue Aug  6 12:53:32 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DSCHit_factory_
#define _DSCHit_factory_

#include <JANA/JFactory.h>
#include "DSCHit.h"

class DSCHit_factory:public jana::JFactory<DSCHit>{
	public:
		DSCHit_factory(){};
		~DSCHit_factory(){};

		// overall scale factors
		double a_scale;
		double t_scale;
		double tdc_scale;

		// calibration constants stored by channel
		vector<double>  a_gains;
		vector<double>  a_pedestals;
		vector<double>  adc_time_offsets;
		vector<double>  tdc_time_offsets;

		//map<string,double>  propogation_corr_factors;
		//double<string,double>  attenuation_corr_factors;
		
		double DELTA_T_ADC_TDC_MAX;

		DSCHit* FindMatch(int sector, double T);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.


};

#endif // _DSCHit_factory_

