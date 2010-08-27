// $Id$
//
//    File: DEventProcessor_trackres_tree.h
// Created: Tue Apr  7 14:54:33 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DEventProcessor_trkres_tree_
#define _DEventProcessor_trkres_tree_

#include <vector>
using namespace std;

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
#include <TRACKING/DMCThrown.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCHit.h>

#include "trackres.h"

class DMCTrajectoryPoint;
class DCoordinateSystem;

class DEventProcessor_trkres_tree:public JEventProcessor{
	public:
		DEventProcessor_trkres_tree();
		~DEventProcessor_trkres_tree();
		const char* className(void){return "DEventProcessor_trackres_tree";}

		trackres *trkres_ptr, trkres;
		TTree *ttrkres;
		const DMagneticFieldMap *bfield;

		pthread_mutex_t mutex;
		
		double SIGMA_CDC;
		double SIGMA_FDC_ANODE;
		double SIGMA_FDC_CATHODE;
		
		class meas_t{
			public:
				double s;
				double err;
				double errc;
				double radlen;
				double B;
				const DMCTrajectoryPoint* traj;
		};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		
		const DMCTrajectoryPoint* FindTrajectoryPoint(const DCoordinateSystem *wire, double &radlen, double &s, vector<const DMCTrajectoryPoint*> trajpoints);
		void GetPtRes(vector<meas_t> &meas, double &deltak, double &pt_res);
		void GetThetaRes(vector<meas_t> &meas, double &theta_res);

};

#endif // _DEventProcessor_trackres_tree_

