// $Id$
//
//    File: JEventProcessor_BCAL_TDC_Timing.h
// Created: Tue Jul 28 10:55:56 EDT 2015
// Creator: mstaib (on Linux egbert 2.6.32-504.30.3.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_TDC_Timing_
#define _JEventProcessor_BCAL_TDC_Timing_

#include <JANA/JEventProcessor.h>

class JEventProcessor_BCAL_TDC_Timing:public jana::JEventProcessor{
	public:
		JEventProcessor_BCAL_TDC_Timing();
		~JEventProcessor_BCAL_TDC_Timing();
		const char* className(void){return "JEventProcessor_BCAL_TDC_Timing";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_BCAL_TDC_Timing_

