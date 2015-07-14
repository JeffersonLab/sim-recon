// $Id$
//
//    File: DEventProcessor_trackeff_hists2.h
// Created: Wed Oct 10 13:30:37 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DEventProcessor_trackeff_hists2_
#define _DEventProcessor_trackeff_hists2_

#include <pthread.h>
#include <map>
#include <iostream>
#include <cmath>

using std::map;

#include <TThread.h>
#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>
#include <DANA/DApplication.h>

#include <PID/DKinematicData.h>
#include <PID/DParticleID.h>
#include <TRACKING/DTrackWireBased.h>
#include <TRACKING/DTrackTimeBased.h>
#include <TRACKING/DTrackCandidate.h>
#include <TRACKING/DMCThrown.h>
#include <TRACKING/DMCTrajectoryPoint.h>

#include <FDC/DFDCGeometry.h>
#include <CDC/DCDCTrackHit.h>
#include <FDC/DFDCPseudo.h>
#include <FDC/DFDCHit.h>

#include <DVector2.h>
#include <particleType.h>

#include "track2.h"
#include "DTrackingResolution.h"
#include "DTrackingResolutionGEANT.h"

class DReferenceTrajectory;
class DCoordinateSystem;

class DEventProcessor_trackeff_hists2:public JEventProcessor{

	public:
		DEventProcessor_trackeff_hists2();
		~DEventProcessor_trackeff_hists2();

		TTree *trkeff;
		track2 trk;
		track2 *trk_ptr;

		bool Search_ChargedTrackHypotheses(JEventLoop *loop, int eventnumber, const DMCThrown *mcthrown);
		bool Search_WireBasedTracks(JEventLoop *loop, int eventnumber, const DMCThrown *mcthrown);
		bool Search_TrackCandidates(JEventLoop *loop, int eventnumber, const DMCThrown *mcthrown);


	private:
		DParticleID *dPIDAlgorithm;
		bool use_rt_thrown;

		jerror_t init(void);	///< Invoked via DEventProcessor virtual method
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
		jerror_t erun(void);					///< Invoked via DEventProcessor virtual method
		jerror_t fini(void);					///< Invoked via DEventProcessor virtual method

		bool isReconstructable(const DMCThrown *mcthrown, vector<const DMCTrajectoryPoint*> &mctrajpoints);

		DTrackingResolution *trkres;
		pthread_mutex_t mutex;
		DReferenceTrajectory *rt_thrown;
		
		double CDCZmin, CDCZmax;
		
		int DEBUG;
		
		void FindLR(vector<const DCoordinateSystem*> &wires, const DReferenceTrajectory *crt, vector<int> &LRhits);
		void FindLR(vector<const DCoordinateSystem*> &wires, vector<const DMCTrajectoryPoint*> &trajpoints, vector<int> &LRhits);
};

#endif // _DEventProcessor_trackeff_hists2_

