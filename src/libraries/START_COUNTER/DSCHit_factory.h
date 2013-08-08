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

		// Theses are placeholders for calibration constants.
		// In the "real" implmentation, we will probably need
		// a map of these indexed by the spcific channel
		double a_scale;
		double a_pedestal;

		double t_scale;
		double t_offset;

		double tdc_scale;
		double tdc_offset;
		
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

