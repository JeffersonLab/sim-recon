// $Id$
//
//    File: DParticleSet_factory.h
// Created: Tue Mar 15 11:17:35 EDT 2011
// Creator: staylor (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DParticleSet_factory_
#define _DParticleSet_factory_

#include <JANA/JFactory.h>
#include "DParticleSet.h"

class DParticleSet_factory:public jana::JFactory<DParticleSet>{
	public:
		DParticleSet_factory(){};
		~DParticleSet_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DParticleSet_factory_

