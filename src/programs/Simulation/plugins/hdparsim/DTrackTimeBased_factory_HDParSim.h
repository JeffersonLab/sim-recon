// $Id$
//
//    File: DTrackTimeBased_factory_HDParSim.h
// Created: Fri Feb 19 16:08:15 EST 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _DTrackTimeBased_factory_HDParSim_
#define _DTrackTimeBased_factory_HDParSim_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackTimeBased.h>

#include "DTrackingResolution.h"

class DTrackTimeBased_factory_HDParSim:public jana::JFactory<DTrackTimeBased>{
	public:
		DTrackTimeBased_factory_HDParSim();
		~DTrackTimeBased_factory_HDParSim(){};
		const char* Tag(void){return "HDParSim";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		bool APPLY_EFFICIENCY_CHARGED;

		DTrackingResolution *res;
};

#endif // _DTrackTimeBased_factory_HDParSim_

