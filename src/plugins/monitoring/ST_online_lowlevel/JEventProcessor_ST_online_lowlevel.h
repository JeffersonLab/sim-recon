// $Id$
//
//    File: JEventProcessor_ST_online_lowlevel.h
// Created: Fri Jun 19 13:21:45 EDT 2015
// Creator: mkamel (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_ST_online_lowlevel_
#define _JEventProcessor_ST_online_lowlevel_

#include <JANA/JEventProcessor.h>

class JEventProcessor_ST_online_lowlevel:public jana::JEventProcessor{
	public:
		JEventProcessor_ST_online_lowlevel();
		~JEventProcessor_ST_online_lowlevel();
		const char* className(void){return "JEventProcessor_ST_online_lowlevel";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_ST_online_lowlevel_

