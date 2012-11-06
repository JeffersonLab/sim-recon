// $Id$
//
//    File: DPhoton_factory_HDParSim.h
// Created: Tue Feb  3 11:29:30 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DPhoton_factory_HDParSim_
#define _DPhoton_factory_HDParSim_

#include <JANA/JFactory.h>
#include "PID/DPhoton.h"

#include "DTrackingResolution.h"

class DPhoton_factory_HDParSim:public jana::JFactory<DPhoton>{
	public:
		DPhoton_factory_HDParSim();
		~DPhoton_factory_HDParSim(){};
		const char* Tag(void){return "HDParSim";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		bool APPLY_EFFICIENCY_PHOTON;

		DTrackingResolution *res;
};

#endif // _DPhoton_factory_HDParSim_

