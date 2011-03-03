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
#include <PID/DParticleID.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALPhoton.h>
#include <TOF/DTOFPoint.h>
#include <TRACKING/DMCThrown.h>
#include <START_COUNTER/DSCHit.h>

class DTrackWireBased;
class DTrackHitSelector;
class DParticleID;

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
  bool PID_FORCE_TRUTH;
  DTrackFitter *fitter;
  DParticleID *pid_algorithm;
  vector<DReferenceTrajectory*> rtv;
 
  // Optional debugging histograms
  TH1F *fom_chi2_trk;
  TH1F *fom_chi2_dedx;
  TH1F *fom;
  TH1F *hitMatchFOM;
  TH2F *chi2_trk_mom;
 
  void FilterDuplicates(void);  
  double GetTruthMatchingFOM(int trackIndex,DTrackTimeBased *dtrack,vector<const DMCThrown*>mcthrowns);
  void GetThrownIndex(const DKinematicData *kd, int &MAX_TRACKS, double &f, int &track);

  // Geometry
  const DGeometry *geom;

  double mPathLength,mEndTime,mStartTime,mFlightTime;
  DetectorSystem_t mDetector, mStartDetector;
  
  // start counter geometry
  double sc_light_guide_length_cor;
  double sc_angle_cor;
  vector<DVector3>sc_pos;
  vector<DVector3>sc_norm;

 
};

#endif // _DTrackTimeBased_factory_

