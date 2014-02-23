// $Id$
//
//    File: DEventProcessor_fdc_hists.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DEventProcessor_fdc_hists_
#define _DEventProcessor_fdc_hists_

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

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCGeometry.h>
#include <FCAL/DFCALShower.h>
#include <FDC/DFDCPseudo.h>
#include <FDC/DFDCIntersection.h>
#include <DMatrixSIMD.h>

#include "FDC_branch.h"
#include "FDChit_branch.h"

class DCDCTrackHit;

typedef struct{
  double w,s,cosa,sina;
  double t,z;
  int wire,layer;
}wire_t;

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


class DEventProcessor_fdc_hists:public JEventProcessor{

	public:
		DEventProcessor_fdc_hists();
		~DEventProcessor_fdc_hists();

		TDirectory *dir;
		TTree *fdctree;
		FDC_branch fdc;
		FDC_branch *fdc_ptr;
		TTree *fdchittree;
		FDChit_branch fdchit;
		FDChit_branch *fdchit_ptr;
		TBranch *fdcbranch, *fdchitbranch;

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
		vector<vector<DFDCWire*> >fdcwires;

		
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method
	
		DMatrix4x1 FitLine(vector<const DFDCPseudo*> &fdchits);
		DMatrix4x1 FitLine(vector<const DFDCPseudo*> &fdchits,
				   double &var_x,double &cov_x_tx,
				   double &var_tx,double &chi2x,
				   double &var_y,double &cov_y_ty,
				   double &var_ty,double &chi2y);
		jerror_t DoFilter(DMatrix4x1 &S,
				  vector<const DFDCPseudo*> &fdchits);
		jerror_t KalmanFilter(double anneal_factor,
				      DMatrix4x1 &S,DMatrix4x4 &C,
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
		TH1F *Hprob,*Htime_prob;
		TH2F *Hures_vs_layer;	
		TH2F *Hcand_ty_vs_tx,*Htime_ty_vs_tx,*Hty_vs_tx;
		TH2F *Hdrift_time;
		TH1F *Hdrift_integral;
		TH2F *Hres_vs_drift_time,*Hvres_vs_layer;
		TH3F *Htime_y_vs_x;
		TH2F *Hqratio_vs_wire,*Hdelta_z_vs_wire;
		TH1F *Hxshift,*Hyshift,*Hphishift;
		TH1F *Hxcand_prob,*Hycand_prob;
		TH1F *Hreduced_chi2,*Hbcal_match;
		TH2F *Hdv_vs_dE;
		TH1F *Hz_target,*Hfcal_match,*Htheta,*HdEdx;

		double mT0;
		double target_to_fcal_distance;
		double fdc_drift_table[140];
		DMatrix4x1 Zero4x1;
		DMatrix4x4 Zero4x4;

		double endplate_z;
		int myevt;

		bool DoAlign;

		vector<align_t>alignments;
};

#endif // _DEventProcessor_fdc_hists_

