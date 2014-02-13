// $Id$
//
//    File: DTOFHit_factory.h
// Created: Wed Aug  7 09:30:17 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DTOFHit_factory_
#define _DTOFHit_factory_

#include <JANA/JFactory.h>
#include "DTOFHit.h"

class DTOFHit_factory:public jana::JFactory<DTOFHit>{
	public:
		DTOFHit_factory(){};
		~DTOFHit_factory(){};

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

		DTOFHit* FindMatch(int plane, int bar, int end, double T);

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DTOFHit_factory_

