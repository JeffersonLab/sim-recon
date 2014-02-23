// $Id$
//
//    File: DBCALHit_factory.h
// Created: Tue Aug  6 09:26:13 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DBCALHit_factory_
#define _DBCALHit_factory_

#include <JANA/JFactory.h>
#include "DBCALHit.h"

class DBCALHit_factory:public jana::JFactory<DBCALHit>{
	public:
		DBCALHit_factory(){};
		~DBCALHit_factory(){};

		// Theses are placeholders for calibration constants.
		// In the "real" implmentation, we will probably need
		// a map of these indexed by the spcific channel
		double a_scale;
		double a_pedestal;

		double t_scale;
		double t_offset;
		

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DBCALHit_factory_

