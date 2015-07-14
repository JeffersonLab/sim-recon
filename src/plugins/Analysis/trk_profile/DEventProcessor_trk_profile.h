// $Id$
//
//    File: DEventProcessor_trk_profile.h
// Created: Wed Jan 12 08:02:32 EST 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.6.0 i386)
//

#ifndef _DEventProcessor_trk_profile_
#define _DEventProcessor_trk_profile_

#include <JANA/JEventProcessor.h>

class DEventProcessor_trk_profile:public jana::JEventProcessor{
	public:
		DEventProcessor_trk_profile();
		~DEventProcessor_trk_profile();
		const char* className(void){return "DEventProcessor_trk_profile";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DEventProcessor_trk_profile_

