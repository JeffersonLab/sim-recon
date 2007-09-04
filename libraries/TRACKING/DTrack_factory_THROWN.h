// $Id$
//
//    File: DTrack_factory_THROWN.h
// Created: Mon Sep  3 19:57:11 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#ifndef _DTrack_factory_THROWN_
#define _DTrack_factory_THROWN_

#include <JANA/JFactory.h>
#include "DTrack.h"

class DTrack_factory_THROWN:public JFactory<DTrack>{
	public:
		DTrack_factory_THROWN(){};
		~DTrack_factory_THROWN(){};
		const string toString(void);
		const char* Tag(void){return "THROWN";}

	private:
		//jerror_t init(void);						///< Called once at program start.
		//jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.

		double SampleGaussian(double sigma);
		
		vector<DReferenceTrajectory*> rt;
		vector<DMatrixDSym*> cov;
};

#endif // _DTrack_factory_THROWN_

