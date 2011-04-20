// $Id$
//
//    File: DPhysicsEvent_factory.h
// Created: Wed Aug  4 10:37:55 EDT 2010
// Creator: davidl (on Darwin eleanor.jlab.org 10.4.0 i386)
//

#ifndef _DPhysicsEvent_factory_
#define _DPhysicsEvent_factory_

#include <JANA/JFactory.h>

#include <TDirectoryFile.h>

#include <PID/DPhysicsEvent.h>
#include <TRACKING/DHoughFind.h>

class DPhysicsEvent_factory:public jana::JFactory<DPhysicsEvent>{
	public:
		DPhysicsEvent_factory(){};
		~DPhysicsEvent_factory(){};
		
	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPhysicsEvent_factory_

