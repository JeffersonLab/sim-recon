// $Id$
//
//    File: DL3Trigger_factory.h
// Created: Wed Jul 31 14:34:24 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#ifndef _DL3Trigger_factory_
#define _DL3Trigger_factory_

#include <JANA/JFactory.h>
#include "DL3Trigger.h"

class DL3Trigger_factory:public jana::JFactory<DL3Trigger>{
	public:
		DL3Trigger_factory(){};
		~DL3Trigger_factory(){};

		double FRACTION_TO_KEEP;
		bool DO_WIRE_BASED_TRACKING;
		bool DO_BCAL_CLUSTER;
		uint32_t L1_TRIG_MASK;
		uint32_t L1_FP_TRIG_MASK;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DL3Trigger_factory_

