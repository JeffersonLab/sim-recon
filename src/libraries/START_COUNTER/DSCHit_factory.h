// $Id$
//
//    File: DSCHit_factory.h
// Created: Tue Aug  6 12:53:32 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DSCHit_factory_
#define _DSCHit_factory_

#include <JANA/JFactory.h>
#include "TTab/DTranslationTable.h"
#include "DSCHit.h"


class DSCHit_factory:public jana::JFactory<DSCHit>{
	public:
		DSCHit_factory(){};
		~DSCHit_factory(){};

		// overall scale factors
		double a_scale;
		double t_scale;
                double t_min;
		double tdc_scale;

		// calibration constants stored by channel
		vector<double>  a_gains;
		vector<double>  a_pedestals;
		vector<double>  adc_time_offsets;
		vector<double>  tdc_time_offsets;

		//map<string,double>  propogation_corr_factors;
		//double<string,double>  attenuation_corr_factors;
		
		double DELTA_T_ADC_TDC_MAX;

		// geometry information
		static const int MAX_SECTORS = 30.;

		DSCHit* FindMatch(int sector, double T);

		const double GetConstant(const vector<double>  &the_table,
					 const int in_sector) const;
		const double GetConstant(const vector<double>  &the_table,
					 const DSCDigiHit *the_digihit) const;
		const double GetConstant(const vector<double>  &the_table,
					 const DSCHit *the_hit) const;
		//const double GetConstant(const vector<double>  &the_table,
		//			 const DTranslationTable *ttab,
		//			 const int in_rocid, const int in_slot, const int in_channel) const;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.


};

#endif // _DSCHit_factory_

