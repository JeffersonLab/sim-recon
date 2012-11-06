// $Id$
//
//    File: DVertex_factory_THROWN.h
// Created: Tue Apr  6 17:01:54 EDT 2010
// Creator: pmatt (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DVertex_factory_THROWN_
#define _DVertex_factory_THROWN_

#include <JANA/JFactory.h>
#include <DVertex.h>
#include <TRACKING/DMCThrown.h>
#include <PID/DChargedTrack.h>

class DVertex_factory_THROWN : public jana::JFactory<DVertex>{
	public:
		DVertex_factory_THROWN(){};
		~DVertex_factory_THROWN(){};
		const char* Tag(void){return "THROWN";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double dTargetCenter;
};

#endif // _DVertex_factory_THROWN_

