// $Id$
//
//    File: JEventProcessor_ST_online_tracking.h
// Created: Fri Jun 19 13:22:21 EDT 2015
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_ST_online_tracking_
#define _JEventProcessor_ST_online_tracking_

#include <JANA/JEventProcessor.h>
#include <START_COUNTER/DSCHit.h>
#include <RF/DRFTDCDigiTime.h>
#include <RF/DRFTime_factory.h>
#include <PID/DEventRFBunch.h>
#include <PID/DParticleID.h>
#include <TRACKING/DTrackFitter.h>

#include "TF1.h"
#include "TH1D.h"
#include "TObjArray.h"
#include "TMath.h"
class JEventProcessor_ST_online_tracking:public jana::JEventProcessor{
	public:
		JEventProcessor_ST_online_tracking();
		~JEventProcessor_ST_online_tracking();
		const char* className(void){return "JEventProcessor_ST_online_tracking";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
const DParticleID* dParticleID;
};

#endif // _JEventProcessor_ST_online_tracking_

