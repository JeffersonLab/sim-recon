// $Id$
//
//    File: DTrackTimeBased_factory_Kalman.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DTrackTimeBased_factory_Kalman_
#define _DTrackTimeBased_factory_Kalman_
#include <TH2F.h>
#include <TROOT.h>

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitter.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALPhoton.h>
#include <TOF/DTOFPoint.h>
#include <START_COUNTER/DSCHit.h>

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
  bool DEBUG_HISTS;
  double MOMENTUM_CUT_FOR_DEDX;
  double MOMENTUM_CUT_FOR_PROTON_ID;
  DTrackFitter *fitter;
  vector<DReferenceTrajectory*> rtv;

  void FilterDuplicates(void);  
  double GetFOM(DTrackTimeBased *dtrack,
		vector<const DBCALShower*>bcal_clusters,
		vector<const DFCALPhoton*>fcal_clusters,
		vector<const DTOFPoint*>tof_points,
		vector<const DSCHit *>sc_hits);
  double MatchToSC(DTrackTimeBased *track,
		   vector<const DSCHit *>sc_hits);
  double MatchToTOF(DTrackTimeBased *track,
		    vector<const DTOFPoint*>tof_points);
  double MatchToBCAL(DTrackTimeBased *track,
		     vector<const DBCALShower*>bcal_clusters);
  
  double mPathLength,mEndTime,mStartTime,mFlightTime;
  DetectorSystem_t mDetector, mStartDetector;

  // Geometry
  const DGeometry *geom;
  // Start counter geometry 
  double sc_light_guide_length;
  double sc_costheta; // cos(theta) of bent part
  vector<DVector3>sc_pos;
  vector<DVector3>sc_norm;

  // Debug histograms
  TH2F *HBCALdTime_vs_E,*HBCALdTime_vs_E_scaled,*HBCALdTime,*HBCALPull;
  TH2F *HTOFdTime,*HTOFPull,*HTOFdTime_vs_beta,*HBCALdTime_vs_beta;
  TH2F *HdEdxDiff,*HdEdxPull,*HdEdxDiff_vs_beta;
};

#endif // _DTrackTimeBased_factory_Kalman_

