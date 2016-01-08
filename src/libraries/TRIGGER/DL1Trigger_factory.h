// $Id$
//
//    File: DL1Trigger_factory.h
// Created: Fri Jan  8 10:57:58 EST 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DL1Trigger_factory_
#define _DL1Trigger_factory_

#include <JANA/JFactory.h>
#include "DL1Trigger.h"

class DL1Trigger_factory:public jana::JFactory<DL1Trigger>{
	public:
		DL1Trigger_factory(){};
		~DL1Trigger_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DL1Trigger_factory_

