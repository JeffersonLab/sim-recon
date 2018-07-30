// $Id$
//
//    File: DTrackTimeBased_factory.h
// Created: Thu Sep  4 14:02:44 EDT 2008
// Creator: davidl (on Darwin harriet.jlab.org 8.11.1 i386)
//

#ifndef _DTrackTimeBased_factory_
#define _DTrackTimeBased_factory_

#include <memory>

#include <TH2.h>

#include <JANA/JFactory.h>
#include <PID/DParticleID.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <TOF/DTOFPoint.h>
#include <CDC/DCDCHit.h>
#include <START_COUNTER/DSCHit.h>
#include "PID/DParticleID.h"

#include "DMCThrown.h"
#include "DTrackFitter.h"
#include "DTrackTimeBased.h"
#include "DReferenceTrajectory.h"

class DTrackWireBased;
class DTrackHitSelector;
class DParticleID;

using namespace jana;
using namespace std;

class DTrackTimeBased_factory:public jana::JFactory<DTrackTimeBased>{
  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *loop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *loop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.

 
  
  bool DEBUG_HISTS;
  int DEBUG_LEVEL;
  double MOMENTUM_CUT_FOR_DEDX;
  double MOMENTUM_CUT_FOR_PROTON_ID;
  bool PID_FORCE_TRUTH;
  unsigned int MIN_CDC_HITS_FOR_TB_FORWARD_TRACKING;
  bool BYPASS_TB_FOR_FORWARD_TRACKS;
//  bool SKIP_MASS_HYPOTHESES_TIMEBASED;
  bool USE_HITS_FROM_WIREBASED_FIT;

  DTrackFitter *fitter;
  const DParticleID* pid_algorithm;
  vector<int> mass_hypotheses_positive;
  vector<int> mass_hypotheses_negative;
 
  // Optional debugging histograms
  TH1F *fom_chi2_trk;
  TH1F *fom;
  TH1F *hitMatchFOM;
  TH2F *chi2_trk_mom;
  TH2F *Hstart_time;
 
  void FilterDuplicates(void);  
  double GetTruthMatchingFOM(int trackIndex,DTrackTimeBased *dtrack,vector<const DMCThrown*>mcthrowns);
  int GetThrownIndex(vector<const DMCThrown*>& locMCThrowns, const DKinematicData *kd, double &f);

  void CreateStartTimeList(const DTrackWireBased *track,
			   vector<const DSCHit*>&sc_hits,
			   vector<const DTOFPoint*>&tof_points,
			   vector<const DBCALShower*>&bcal_showers,	  
			   vector<const DFCALShower*>&fcal_showers,
			   vector<DTrackTimeBased::DStartTime_t>&start_times);
  bool DoFit(const DTrackWireBased *track,
	     vector<DTrackTimeBased::DStartTime_t>&start_times,
	     JEventLoop *loop,double mass);  

  void AddMissingTrackHypothesis(vector<DTrackTimeBased*>&tracks_to_add,
				 const DTrackTimeBased *src_track,
				 double my_mass,double q,
				 JEventLoop *loop);
  bool InsertMissingHypotheses(JEventLoop *loop);
  void CorrectForELoss(DVector3 &position,DVector3 &momentum,double q,double my_mass);
  void AddMissingTrackHypotheses(unsigned int mass_bits,
				 vector<DTrackTimeBased*>&tracks_to_add,
				 vector<DTrackTimeBased *>&hypotheses,
				 double q,bool flipped_charge,JEventLoop *loop);

  // Geometry
  const DGeometry *geom;

//  double mPathLength,mEndTime,mStartTime,mFlightTime;
  double mStartTime;
//  DetectorSystem_t mDetector, mStartDetector;
  DetectorSystem_t mStartDetector;
  int mNumHypPlus,mNumHypMinus;
  bool dIsNoFieldFlag;
  bool USE_SC_TIME; // use start counter hits for t0
  bool USE_FCAL_TIME; // use fcal hits for t0
  bool USE_BCAL_TIME; // use bcal hits for t0
  bool USE_TOF_TIME; // use tof hits for t0
  bool SKIP_MASS_HYPOTHESES_WIRE_BASED;
//  double SC_DPHI_CUT_WB;

  // start counter geometry
//  double sc_light_guide_length_cor;
//  double sc_angle_cor;
  vector<DVector3>sc_pos;
  vector<DVector3>sc_norm;

  int myevt;
};

#endif // _DTrackTimeBased_factory_

