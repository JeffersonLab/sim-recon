// $Id$
//
//    File: JEventProcessor_cdc_efficiency.h
// Created: Tue Sep  9 15:41:38 EDT 2014
// Creator: hdcdcops (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_cdc_efficiency_
#define _JEventProcessor_cdc_efficiency_

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
#include <PID/DDetectorMatches.h>
#include <CDC/DCDCHit.h>
#include <DAQ/DEPICSvalue.h>

class JEventProcessor_cdc_efficiency:public jana::JEventProcessor{
	public:
		JEventProcessor_cdc_efficiency();
		~JEventProcessor_cdc_efficiency();
		const char* className(void){return "JEventProcessor_cdc_efficiency";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
        void TimeBasedAnalysis(JEventLoop *);
        DGeometry * dgeom;
		vector< vector< DCDCWire * > > cdcwires; // CDC Wires Referenced by [ring][straw]
        //pthread_mutex_t mutex;
		int BFIELD, COSMICS, USE_TIMEBASEDTRACKS, REQUIRE_BEAM, BEAM_EVENTS_TO_KEEP;
        double DOCACUT, BEAM_CURRENT;
        //const DMagneticFieldMap *bfield;
       // DReferenceTrajectory * rt;
};

#endif // _JEventProcessor_cdc_efficiency_

