// $Id$
//
//    File: DTrackTimeBased_factory.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DTrackTimeBased_factory_
#define _DTrackTimeBased_factory_

#include <TH2.h>

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitter.h>
#include <BCAL/DBCALShower.h>
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
  
  bool DEBUG_HISTS;
  int DEBUG_LEVEL;
  double MOMENTUM_CUT_FOR_DEDX;
  double MOMENTUM_CUT_FOR_PROTON_ID;
  DTrackFitter *fitter;
  vector<DReferenceTrajectory*> rtv;
 
	// Optional debugging histograms
	TH1F *fom_tdiff_bcal;
	TH1F *fom_tdiff_tof;
	TH1F *fom_chi2_trk;
	TH1F *fom_chi2_dedx;
	TH1F *fom_chi2_tof;
	TH1F *fom_chi2_bcal;
	TH1F *time_based_start;
  
  void FilterDuplicates(void);  
  double GetFOM(DTrackTimeBased *dtrack,
		vector<const DBCALShower*>bcal_clusters,
		vector<const DFCALPhoton*>fcal_clusters,
		vector<const DTOFPoint*>tof_points);
  double MatchToTOF(DTrackTimeBased *track,
		    vector<const DTOFPoint*>tof_points);
  double MatchToBCAL(DTrackTimeBased *track,
		     vector<const DBCALShower*>bcal_clusters);

  // Geometry
  const DGeometry *geom;

  double mPathLength,mEndTime,mStartTime,mFlightTime;
  DetectorSystem_t mDetector, mStartDetector;
 
};

#endif // _DTrackTimeBased_factory_

