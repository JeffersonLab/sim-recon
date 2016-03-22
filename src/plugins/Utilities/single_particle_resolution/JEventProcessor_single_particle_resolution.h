// $Id$
//
//    File: JEventProcessor_single_particle_resolution.h
// Created: Mon Feb  8 15:12:19 EST 2016
// Creator: dalton (on Linux gluon02.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_single_particle_resolution_
#define _JEventProcessor_single_particle_resolution_

#include <JANA/JEventProcessor.h>

class JEventProcessor_single_particle_resolution:public jana::JEventProcessor{
	public:
		JEventProcessor_single_particle_resolution();
		~JEventProcessor_single_particle_resolution();
		const char* className(void){return "JEventProcessor_single_particle_resolution";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		int VERBOSE;
};

#endif // _JEventProcessor_single_particle_resolution_

