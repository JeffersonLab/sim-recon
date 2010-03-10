// $Id$
//
//    File: DParticle_factory_HDParSim.h
// Created: Tue Feb  3 11:27:46 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DParticle_factory_HDParSim_
#define _DParticle_factory_HDParSim_

#include <JANA/JFactory.h>
#include "PID/DParticle.h"

#include "DTrackingResolution.h"

class DParticle_factory_HDParSim:public jana::JFactory<DParticle>{
	public:
		DParticle_factory_HDParSim();
		~DParticle_factory_HDParSim(){};
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

#endif // _DParticle_factory_HDParSim_

