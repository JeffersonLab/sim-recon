// $Id$
//
//    File: DTrackTimeBased_factory_Kalman.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DTrackTimeBased_factory_Kalman_
#define _DTrackTimeBased_factory_Kalman_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitter.h>

class DTrackWireBased;
class DTrackHitSelector;

#include "DTrackTimeBased.h"

/// Time based tracks

class DTrackTimeBased_factory_Kalman:public jana::JFactory<DTrackTimeBased>{
	public:
		DTrackTimeBased_factory_Kalman(){};
		~DTrackTimeBased_factory_Kalman(){};
		const char* Tag(void){return "Kalman";}


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *loop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		int DEBUG_LEVEL;
		double MOMENTUM_CUT_FOR_DEDX;
		DTrackFitter *fitter;
		vector<DReferenceTrajectory*> rtv;
		vector<double> mass_hypotheses;

		DTrackTimeBased* MakeDTrackTimeBased(const DTrackWireBased *track);
		double GetFOM(DTrackTimeBased *dtrack);
		double GetRangeOutFOM(DTrackTimeBased *dtrack);
};

#endif // _DTrackTimeBased_factory_Kalman_

