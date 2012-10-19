// $Id$
//
//    File: DEventProcessor_dc_alignment.h
// Created: Thu Oct 18 17:15:41 EDT 2012
// Creator: staylor (on Linux ifarm1102 2.6.18-274.3.1.el5 x86_64)
//

#ifndef _DEventProcessor_dc_alignment_
#define _DEventProcessor_dc_alignment_

#include <pthread.h>
#include <map>
#include <vector>
#include <deque>
using std::map;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JCalibration.h>

#include <PID/DKinematicData.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCGeometry.h>
#include <BCAL/DBCALShower.h>
#include <FCAL/DFCALShower.h>
#include <FDC/DFDCPseudo.h>
#include <DMatrixSIMD.h>
#include <DVector2.h>

#include "FDC_branch.h"

typedef struct{
  DMatrix4x1 S;
  DMatrix4x4 J;
  DMatrix4x1 Skk;
  DMatrix4x4 Ckk;
  double z,t;
  int num_hits;
  unsigned int h_id;
}trajectory_t;

typedef struct{
  unsigned int id;
  double ures,vres;
  DMatrix4x1 S;
  DMatrix4x4 C;
  DMatrix2x2 R;
  DMatrix4x2 H_T;
  DMatrix2x4 H;
  double doca;
  double drift,drift_time;
}update_t;

typedef struct{
  unsigned int id;
  DMatrix4x1 S;
  DMatrix4x4 C; 
  DMatrix4x1 H_T;
  DMatrix1x4 H;
  double ures,vres;
  double R;
  double drift,drift_time;
}strip_update_t;



typedef struct{
  DMatrix3x1 A;
  DMatrix3x3 E;  
}align_t;


typedef struct{
  bool matched;
  DMatrix4x1 S;
  vector<const DFDCPseudo *>hits;
}segment_t;

#define EPS 1e-3
#define ITER_MAX 20
#define ADJACENT_MATCH_RADIUS 1.0
#define MATCH_RADIUS 2.0

class DEventProcessor_dc_alignment:public jana::JEventProcessor{
 public:
  DEventProcessor_dc_alignment();
  ~DEventProcessor_dc_alignment();
  const char* className(void){return "DEventProcessor_dc_alignment";}

  TDirectory *dir;
  TTree *fdctree;
  FDC_branch fdc;
  FDC_branch *fdc_ptr;
  TBranch *fdcbranch;

  enum track_type{
    kWireBased,
    kTimeBased,
  };
  enum state_vector{
    state_x,
    state_y,
    state_tx,
    state_ty,
  };
  enum align_parms{
    kDx,
    kDy,
    kDPhi,
  };

  
 private:
  jerror_t init(void);						///< Called once at program start.
  jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
  jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
  jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
  jerror_t fini(void);						///< Called after last event of last event source has been processed.
  
  DMatrix4x1 FitLine(vector<const DFDCPseudo*> &fdchits);
  DMatrix4x1 FitLine(vector<const DFDCPseudo*> &fdchits,
		     double &var_x,double &cov_x_tx,
		     double &var_tx,double &chi2x,
		     double &var_y,double &cov_y_ty,
		     double &var_ty,double &chi2y);
  jerror_t DoFilter(DMatrix4x1 &S,vector<const DFDCPseudo*> &fdchits);
  jerror_t KalmanFilter(double anneal_factor,DMatrix4x1 &S,DMatrix4x4 &C,
			vector<const DFDCPseudo *>&hits,
			deque<trajectory_t>&trajectory,
			vector<strip_update_t>&updates,
			double &chi2,unsigned int &ndof);
  jerror_t KalmanFilter(double anneal_factor,
			DMatrix4x1 &S,DMatrix4x4 &C,
			vector<const DFDCPseudo *>&hits,
			deque<trajectory_t>&trajectory,
			vector<update_t>&updates,
			double &chi2,unsigned int &ndof);	
  jerror_t Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
		  deque<trajectory_t>&trajectory,
		  vector<const DFDCPseudo *>&hits,
		  vector<strip_update_t>updates,
		  vector<strip_update_t>&smoothed_updates);
  jerror_t Smooth(DMatrix4x1 &Ss,DMatrix4x4 &Cs,
		  deque<trajectory_t>&trajectory,
		  vector<const DFDCPseudo *>&hits,
		  vector<update_t>updates,
		  vector<update_t>&smoothed_updates);
  jerror_t SetReferenceTrajectory(double z,DMatrix4x1 &S,
				  deque<trajectory_t>&trajectory,
				  vector<const DFDCPseudo *>&wires);
  jerror_t FindSegments(vector<const DFDCPseudo*>&pseudos,
			vector<segment_t>&segments);
  jerror_t LinkSegments(vector<segment_t>segments[4], 
			vector<vector<const DFDCPseudo *> >&LinkedSegments);
  jerror_t FindOffsets(vector<const DFDCPseudo *>&hits,
		       vector<update_t>smoothed_updates);
  jerror_t FindOffsets(vector<const DFDCPseudo *>&hits,
		       vector<strip_update_t>smoothed_updates);

  double GetDriftDistance(double t);
  double GetDriftVariance(double t);
  
  pthread_mutex_t mutex;

  TH1F *Hprob;
  TH2F *Hures_vs_layer;	
  TH2F *Hdrift_time;
  TH2F *Hres_vs_drift_time,*Hvres_vs_layer;
  TH2F *Hdv_vs_dE;

  double mT0;
  double target_to_fcal_distance;
  double fdc_drift_table[140];
  DMatrix4x1 Zero4x1;
  DMatrix4x4 Zero4x4;

  double endplate_z;
  int myevt;

  vector<align_t>alignments;
  vector<vector<DFDCWire*> >fdcwires;
};

#endif // _DEventProcessor_dc_alignment_

