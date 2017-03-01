// $Id$
//
//    File: JEventProcessor_TAGM_thresh.h
// Created: Tue Feb 28 17:43:55 EST 2017
// Creator: barnes (on Linux gluey.phys.uconn.edu 2.6.32-642.13.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TAGM_thresh_
#define _JEventProcessor_TAGM_thresh_

#include <JANA/JEventProcessor.h>

class JEventProcessor_TAGM_thresh:public jana::JEventProcessor{
	public:
		JEventProcessor_TAGM_thresh();
		~JEventProcessor_TAGM_thresh();
		const char* className(void){return "JEventProcessor_TAGM_thresh";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TAGM_thresh_

