// $Id$
//
//    File: JEventProcessor_CDC_Efficiency.h
// Created: Tue Sep  9 15:41:38 EDT 2014
// Creator: hdcdcops (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_CDC_Efficiency_
#define _JEventProcessor_CDC_Efficiency_

//#include <pthread.h>
#include <map>
#include <vector>
#include <deque>
using namespace std;

#include <TTree.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>

#include <JANA/JFactory.h>
#include <JANA/JEventProcessor.h>
#include <JANA/JEventLoop.h>
#include <JANA/JCalibration.h>

#include <HDGEOMETRY/DGeometry.h>
#include <TRACKING/DTrackCandidate_factory_StraightLine.h>
#include <TRACKING/DReferenceTrajectory.h>
#include <TRACKING/DTrackWireBased.h>
#include <PID/DChargedTrack.h>
#include <PID/DDetectorMatches.h>
#include <CDC/DCDCHit.h>

class JEventProcessor_CDC_Efficiency:public jana::JEventProcessor{
	public:
		JEventProcessor_CDC_Efficiency();
		~JEventProcessor_CDC_Efficiency();
		const char* className(void){return "JEventProcessor_CDC_Efficiency";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
        DGeometry * dgeom;
        bool dIsNoFieldFlag;
		vector< vector< DCDCWire * > > cdcwires; // CDC Wires Referenced by [ring][straw]
        vector<vector<double> >max_sag;
        vector<vector<double> >sag_phi_offset;
        int ChannelFromRingStraw[28][209];
        double DOCACUT;
};

#endif // _JEventProcessor_CDC_Efficiency_

