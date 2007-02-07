// $Id$
//
//    File: DSCHit_factory.h
// Created: Wed Feb  7 10:46:20 EST 2007
// Creator: davidl (on Linux megrez.jlab.org 2.6.9-42.0.2.ELsmp x86_64)
//

#ifndef _DSCHit_factory_
#define _DSCHit_factory_

#include <JANA/JFactory.h>
#include "DSCHit.h"

class DSCHit_factory:public JFactory<DSCHit>{
	public:
		DSCHit_factory(){};
		~DSCHit_factory(){};
		const string toString(void);


	private:
		//jerror_t init(void);						///< Called once at program start.
		//jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		//jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DSCHit_factory_

