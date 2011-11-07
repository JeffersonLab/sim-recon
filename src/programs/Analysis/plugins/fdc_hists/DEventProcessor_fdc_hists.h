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

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>

#include <PID/DKinematicData.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DMCTrackHit.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>
#include <FDC/DFDCGeometry.h>

#include <FDC/DFDCPseudo.h>
#include <FDC/DFDCIntersection.h>
#include <DMatrixSIMD.h>

#include "FDC_branch.h"
#include "FDChit_branch.h"

class DCDCTrackHit;

typedef struct{
  DMatrix4x1 S;
  DMatrix4x4 J;
  DMatrix4x1 Skk;
  DMatrix4x4 Ckk;
  DMatrix4x1 Skkp;
  DMatrix4x4 Ckkp;
  double z,t;
  int num_hits;
  unsigned int h_id;
}trajectory_t;


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
		
		
	private:
		vector<vector<DFDCWire*> >fdcwires;

		
		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		DMatrix4x1 SetSeed(vector<const DFDCIntersection*> fdchits);
		jerror_t Fit(int fit_type,DMatrix4x1 &S,DMatrix4x4 &C,
			     double &chi2,unsigned int &ndof);
		jerror_t Smooth(int fit_type,DMatrix4x1 &Ss,DMatrix4x4 &Cs);
		jerror_t SetReferenceTrajectory(DMatrix4x1 &S,
						int layer_to_skip);
		double GetDriftDistance(double t);
		double GetDriftVariance(double t);

		vector<DFDCPseudo *>points;
		deque<trajectory_t>trajectory;

		pthread_mutex_t mutex;
		TH1F *Hwire_prob,*Htime_prob;
		TH2F *Hwire_res_vs_layer;	
		TH2F *Htime_res_vs_layer;
		TH2F *Hcand_ty_vs_tx,*Htime_ty_vs_tx,*Hwire_ty_vs_tx;
		TH2F *Hres_vs_drift_time;

		double mT0;
		double fdc_drift_table[140];
};

#endif // _DEventProcessor_fdc_hists_

