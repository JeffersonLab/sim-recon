// $Id$
//
//    File: DParticle_factory_THROWN.h
// Created: Sat Oct  4 22:04:56 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DParticle_factory_THROWN_
#define _DParticle_factory_THROWN_

#include <JANA/JFactory.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <HDGEOMETRY/DRootGeom.h>
#include <HDGEOMETRY/DGeometry.h>
#include "DParticle.h"

class DTrackFitter;
class DTrackHitSelector;

class DParticle_factory_THROWN:public jana::JFactory<DParticle>{
	public:
		DParticle_factory_THROWN();
		~DParticle_factory_THROWN(){};
		const char* Tag(void){return "THROWN";}

	private:

		
		//jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DTrackFitter *fitter;
		const DTrackHitSelector *hitselector;
		vector<DReferenceTrajectory*> rt_pool;
		
		DRootGeom *RootGeom;
		DGeometry *geom;
		string MATERIAL_MAP_MODEL;
};

#endif // _DParticle_factory_THROWN_

