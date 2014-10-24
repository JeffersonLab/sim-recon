// $Id$
//
//    File: DPSCHit_factory.h
// Created: Wed Oct 15 16:45:33 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSCHit_factory_
#define _DPSCHit_factory_

#include <JANA/JFactory.h>
#include "DPSCHit.h"

class DPSCHit_factory:public jana::JFactory<DPSCHit>{
	public:
		DPSCHit_factory(){};
		~DPSCHit_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPSCHit_factory_

