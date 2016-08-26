// $Id$
//
//    File: DTrackCandidate_factory_StraightLine.h
// Created: Fri Aug 15 09:14:04 EDT 2014
// Creator: staylor (on Linux ifarm1102 2.6.32-220.7.1.el6.x86_64 x86_64)
//

#ifndef _DTrackCandidate_factory_StraightLine_
#define _DTrackCandidate_factory_StraightLine_

#include <JANA/JFactory.h>
#include "DTrackCandidate.h"
#include "TRACKING/DTrackFinder.h"
#include <TRACKING/DTrackFitter.h>
#include "CDC/DCDCTrackHit.h"
#include "FDC/DFDCPseudo.h"
#include "DMatrixSIMD.h"
#include <deque>
#include <TROOT.h>
#include <TH2F.h>
#include "PID/DParticleID.h"

class DParticleID;

class DTrackCandidate_factory_StraightLine:public jana::JFactory<DTrackCandidate>{
 public:
  DTrackCandidate_factory_StraightLine(){};
  ~DTrackCandidate_factory_StraightLine(){};
  const char* Tag(void){return "StraightLine";}

  enum state_vector{
    state_x,
    state_y,
    state_tx,
    state_ty,
  };

  class trajectory_t{
  public:
  trajectory_t(double z,double t,DMatrix4x1 S,DMatrix4x4 J,DMatrix4x1 Skk,
	       DMatrix4x4 Ckk,unsigned int id=0,unsigned int numhits=0)
      :z(z),t(t),S(S),J(J),Skk(Skk),Ckk(Ckk),id(id),numhits(numhits){}
    double z,t; 
    DMatrix4x1 S;
    DMatrix4x4 J;
    DMatrix4x1 Skk;
    DMatrix4x4 Ckk;
    unsigned int id,numhits;

  };

  typedef struct{
    double resi,err,d,tdrift,s;
  }update_t;

  typedef struct{
    double d,tdrift,s;
    DMatrix4x1 S;
    DMatrix4x4 C;
    DMatrix2x2 V;
  }fdc_update_t;
  

 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  double fdc_drift_distance(double time);
  jerror_t DoFilter(double t0,double start_z,DMatrix4x1 &S,
		    vector<const DFDCPseudo *>&hits);
  jerror_t DoFilter(double t0,double OuterZ,DMatrix4x1 &S,
		    vector<const DCDCTrackHit *>&hits,double dzsign);

  jerror_t SetReferenceTrajectory(double t0,double z,DMatrix4x1 &S,
				  deque<trajectory_t>&trajectory,
				  vector<const DFDCPseudo *>&pseudos);
  jerror_t SetReferenceTrajectory(double t0,double z,DMatrix4x1 &S,
				  deque<trajectory_t>&trajectory,
				  const DCDCTrackHit *last_cdc,double dzsign); 

  jerror_t KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
			vector<const DFDCPseudo *>&hits,
			vector<int>&used_fdc_hits,
			deque<trajectory_t>&trajectory,
			vector<fdc_update_t>&pulls,
			double &chi2,unsigned int &ndof);
  jerror_t KalmanFilter(DMatrix4x1 &S,DMatrix4x4 &C,
			vector<const DCDCTrackHit *>&hits,
			vector<int>&used_hits,
			deque<trajectory_t>&trajectory,
			vector<update_t>&pulls,
			double &chi2,unsigned int &ndof,bool timebased=false);
  jerror_t Smooth(deque<trajectory_t>&trajectory,
		  vector<fdc_update_t>&updates,
		  vector<const DFDCPseudo *>&hits,
		  DTrackCandidate *cand); 

  double CDCDriftDistance(double t);
  double CDCDriftDistance(double dphi, double delta,double t);
  double CDCDriftVariance(double t);
  unsigned int Locate(vector<double>&xx,double x);

  bool COSMICS,DEBUG_HISTS,USE_FDC_DRIFT_TIMES;
  float CHI2CUT;
  int DO_PRUNING;
  int PLANE_TO_SKIP;

  DTrackFinder *finder;
  const DParticleID* dParticleID;

 // drift time tables
  vector<double>cdc_drift_table;
  vector<double>fdc_drift_table;
  
  // Resolution parameters
  double CDC_RES_PAR1,CDC_RES_PAR2;

  // variables to deal with CDC straw sag
  vector<vector<double> >max_sag;
  vector<vector<double> >sag_phi_offset;
  double long_drift_func[3][3];
  double short_drift_func[3][3];

  // Diagnostic histograms
  TH2F *Hvres;
  
};

#endif // _DTrackCandidate_factory_StraightLine_

