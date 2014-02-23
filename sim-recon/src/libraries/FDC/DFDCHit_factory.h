// $Id$
//
//    File: DFDCHit_factory.h
// Created: Wed Aug  7 11:55:02 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DFDCHit_factory_
#define _DFDCHit_factory_

#include <JANA/JFactory.h>
#include "DFDCHit.h"

class DFDCHit_factory:public jana::JFactory<DFDCHit>{
	public:
		DFDCHit_factory(){};
		~DFDCHit_factory(){};

		// Theses are placeholders for calibration constants.
		// In the "real" implmentation, we will probably need
		// a map of these indexed by the spcific channel
		double a_scale;
		double a_pedestal;

		double t_scale;
		double t_offset;

		double tdc_scale;
		double tdc_offset;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DFDCHit_factory_

