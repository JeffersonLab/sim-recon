// $Id$
//
//    File: DNeutralTrack_factory.h
// Created: Tue Aug  9 14:29:24 EST 2011
// Creator: pmatt (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DNeutralTrack_factory_
#define _DNeutralTrack_factory_

#include <JANA/JFactory.h>
#include <PID/DNeutralTrack.h>
#include <PID/DNeutralShowerCandidate.h>

class DNeutralTrack_factory:public jana::JFactory<DNeutralTrack>{
	public:
		DNeutralTrack_factory(){};
		~DNeutralTrack_factory(){};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *locEventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *locEventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DNeutralTrack_factory_

