// $Id$
//
//    File: JEventProcessor_pedestal_online.h
// Created: Thu Aug  7 09:37:01 EDT 2014
// Creator: dalton (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_pedestal_online_
#define _JEventProcessor_pedestal_online_

#include <TDirectory.h>

#include <JANA/JEventProcessor.h>
#include <JANA/JEvent.h>

class JEventProcessor_pedestal_online:public jana::JEventProcessor{
	public:

		JEventProcessor_pedestal_online();
		~JEventProcessor_pedestal_online();
		const char* className(void){return "JEventProcessor_pedestal_online";}

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

#endif // _JEventProcessor_pedestal_online_

