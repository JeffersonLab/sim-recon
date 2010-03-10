// $Id$
//
//    File: DTrackTimeBased_factory.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DTrackTimeBased_factory_
#define _DTrackTimeBased_factory_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitter.h>
#include <BCAL/DBCALPhoton.h>
#include <FCAL/DFCALPhoton.h>
#include <TOF/DTOFPoint.h>

class DTrackWireBased;
class DTrackHitSelector;

#include "DTrackTimeBased.h"

/// Time based tracks

class DTrackTimeBased_factory:public jana::JFactory<DTrackTimeBased>{
 public:
  DTrackTimeBased_factory(){};
  ~DTrackTimeBased_factory(){};
  
  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *loop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *loop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  int DEBUG_LEVEL;
  double MOMENTUM_CUT_FOR_DEDX;
  double MOMENTUM_CUT_FOR_PROTON_ID;
  DTrackFitter *fitter;
  vector<DReferenceTrajectory*> rtv;
 
  double GetFOM(DTrackTimeBased *dtrack,
		vector<const DBCALPhoton*>bcal_clusters,
		vector<const DFCALPhoton*>fcal_clusters,
		vector<const DTOFPoint*>tof_points);
  double MatchToTOF(DTrackTimeBased *track,
		    vector<const DTOFPoint*>tof_points);
  double MatchToBCAL(DTrackTimeBased *track,
		     vector<const DBCALPhoton*>bcal_clusters);
  // The following routine is not currently used
  double GetRangeOutFOM(DTrackTimeBased *dtrack);
  

  double mPathLength,mEndTime,mStartTime;
  DetectorSystem_t mDetector, mStartDetector;
 
};

#endif // _DTrackTimeBased_factory_

