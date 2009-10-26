// $Id$
//
//    File: DParticle_factory_Kalman.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DParticle_factory_Kalman_
#define _DParticle_factory_Kalman_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitter.h>
#include <BCAL/DBCALPhoton.h>
#include <FCAL/DFCALPhoton.h>
#include <TOF/DTOFPoint.h>

class DTrackTimeBased;
class DTrackHitSelector;

#include "DParticle.h"

/// Time based tracks

class DParticle_factory_Kalman:public jana::JFactory<DParticle>{
	public:
		DParticle_factory_Kalman(){};
		~DParticle_factory_Kalman(){};
		const char* Tag(void){return "Kalman";}


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *loop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DParticle* MakeDParticle(const DTrackTimeBased *track,
					 double mass);

		jerror_t MatchToBCAL(const DTrackTimeBased *track,
				     DReferenceTrajectory *rt, 
				     vector<const DBCALPhoton*>bcal_clusters,
				     vector<bool>&bcal_matches,
				     double &mass);
		jerror_t MatchToFCAL(const DTrackTimeBased *track,
				     DReferenceTrajectory *rt,
				     vector<const DFCALPhoton*>fcal_clusters,
				     vector<bool>&fcal_matches,
				     double &mass); 
		jerror_t MatchToTOF(const DTrackTimeBased *track,
				    DReferenceTrajectory *rt,	
				    vector<const DTOFPoint*>tof_points,
				    double &mass); 
		

		int DEBUG_LEVEL;
		DTrackFitter *fitter;
		vector<DReferenceTrajectory*> rtv;
		double mPathLength,mEndTime,mStartTime;
		DetectorSystem_t mDetector, mStartDetector;
};

#endif // _DParticle_factory_Kalman_

