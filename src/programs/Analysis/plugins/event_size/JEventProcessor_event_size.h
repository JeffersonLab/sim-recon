// $Id$
//
//    File: JEventProcessor_event_size.h
// Created: Tue Jun  7 10:07:36 EDT 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.7.0 i386)
//

#ifndef _JEventProcessor_event_size_
#define _JEventProcessor_event_size_

#include <TTree.h>

#include <JANA/JEventProcessor.h>

#include "Event.h"

class JEventProcessor_event_size:public jana::JEventProcessor{
	public:
		JEventProcessor_event_size();
		~JEventProcessor_event_size();
		const char* className(void){return "JEventProcessor_event_size";}
		
		TTree *evt_tree;
		Event *evt;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		pthread_mutex_t mutex;
};

#endif // _JEventProcessor_event_size_

