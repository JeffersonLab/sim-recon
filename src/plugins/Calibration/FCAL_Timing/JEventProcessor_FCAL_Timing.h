// $Id$
//
//    File: JEventProcessor_FCAL_Timing.h
// Created: Mon Jan  4 13:03:20 EST 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FCAL_Timing_
#define _JEventProcessor_FCAL_Timing_

#include <JANA/JEventProcessor.h>

class JEventProcessor_FCAL_Timing:public jana::JEventProcessor{
	public:
		JEventProcessor_FCAL_Timing();
		~JEventProcessor_FCAL_Timing();
		const char* className(void){return "JEventProcessor_FCAL_Timing";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_FCAL_Timing_

