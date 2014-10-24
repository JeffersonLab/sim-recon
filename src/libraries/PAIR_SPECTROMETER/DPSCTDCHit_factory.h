// $Id$
//
//    File: DPSCTDCHit_factory.h
// Created: Wed Oct 15 16:48:32 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _DPSCTDCHit_factory_
#define _DPSCTDCHit_factory_

#include <JANA/JFactory.h>
#include "DPSCTDCHit.h"

class DPSCTDCHit_factory:public jana::JFactory<DPSCTDCHit>{
	public:
		DPSCTDCHit_factory(){};
		~DPSCTDCHit_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPSCTDCHit_factory_

