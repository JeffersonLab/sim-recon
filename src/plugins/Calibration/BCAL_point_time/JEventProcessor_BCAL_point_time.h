// $Id$
//
//    File: JEventProcessor_BCAL_point_time.h
// Created: Fri Apr  8 12:59:18 EDT 2016
// Creator: dalton (on Linux gluon109.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_point_time_
#define _JEventProcessor_BCAL_point_time_

#include <TDirectory.h>

#include <JANA/JEventProcessor.h>

class JEventProcessor_BCAL_point_time:public jana::JEventProcessor{
	public:
		JEventProcessor_BCAL_point_time();
		~JEventProcessor_BCAL_point_time();
		const char* className(void){return "JEventProcessor_BCAL_point_time";}

		TDirectory *maindir;
		TDirectory *peddir;

	private:
		uint32_t VERBOSE;

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_BCAL_point_time_

