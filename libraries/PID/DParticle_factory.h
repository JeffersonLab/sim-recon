// $Id$
//
//    File: DParticle_factory.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DParticle_factory_
#define _DParticle_factory_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitter.h>
#include <BCAL/DBCALPhoton.h>
#include <FCAL/DFCALPhoton.h>
#include <TOF/DTOFPoint.h>

class DTrackTimeBased;
class DTrackHitSelector;
#include "DParticle.h"

class DParticle_factory:public jana::JFactory<DParticle>{
	public:
		DParticle_factory(){};
		~DParticle_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *loop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		DParticle* MakeDParticle(const DTrackTimeBased *track,
					 double mass);

		jerror_t MatchToBCAL(const DTrackTimeBased *track,
				     vector<const DBCALPhoton*>bcal_clusters,
				     vector<bool>&bcal_matches,
				     double &mass);
		jerror_t MatchToFCAL(const DTrackTimeBased *track,
				     vector<const DFCALPhoton*>fcal_clusters,
				     vector<bool>&fcal_matches,
				     double &mass); 
		jerror_t MatchToTOF(const DTrackTimeBased *track,
				    vector<const DTOFPoint*>tof_points,
				    double &mass); 
		  
		int DEBUG_LEVEL;
		DTrackFitter *fitter;
		vector<DReferenceTrajectory*> rtv;
		double mPathLength,mEndTime,mStartTime;
		DetectorSystem_t mDetector, mStartDetector;
};

#endif // _DParticle_factory_

