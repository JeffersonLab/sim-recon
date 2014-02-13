// $Id$
//
//    File: DTrackWireBased_factory_THROWN.h
// Created: Mon Sep  3 19:57:11 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#ifndef _DTrackWireBased_factory_THROWN_
#define _DTrackWireBased_factory_THROWN_

#include <JANA/JFactory.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <HDGEOMETRY/DRootGeom.h>
#include <HDGEOMETRY/DGeometry.h>
#include <HDGEOMETRY/DMagneticFieldMap.h>
#include "DTrackWireBased.h"

class DTrackFitter;
class DTrackHitSelector;

class DTrackWireBased_factory_THROWN:public JFactory<DTrackWireBased>{
	public:
		DTrackWireBased_factory_THROWN();
		~DTrackWireBased_factory_THROWN(){};
		const char* Tag(void){return "THROWN";}

	private:
		//jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		DTrackFitter *fitter;
		const DTrackHitSelector *hitselector;
		vector<DReferenceTrajectory*> rt_pool;
	
		DRootGeom *RootGeom;
		DGeometry *geom;
		DMagneticFieldMap *bfield;
};

#endif // _DTrackWireBased_factory_THROWN_

