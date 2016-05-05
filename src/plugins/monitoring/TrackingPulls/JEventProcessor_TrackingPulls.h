// $Id$
//
//    File: JEventProcessor_TrackingPulls.h
// Created: Tue Apr 26 09:29:02 EDT 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TrackingPulls_
#define _JEventProcessor_TrackingPulls_

#include <JANA/JEventProcessor.h>
#include "PID/DChargedTrack.h"
#include "TRACKING/DTrackTimeBased.h"

class JEventProcessor_TrackingPulls:public jana::JEventProcessor{
	public:
		JEventProcessor_TrackingPulls();
		~JEventProcessor_TrackingPulls();
		const char* className(void){return "JEventProcessor_TrackingPulls";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
        void FillCDCPulls( DTrackFitter::pull_t , const DTrackTimeBased*);
        void FillFDCPulls( DTrackFitter::pull_t );
        int EXCLUDERING;
        int EXCLUDEPLANE;
        vector<vector<double> >max_sag;
        vector<vector<double> >sag_phi_offset;
};

#endif // _JEventProcessor_TrackingPulls_

