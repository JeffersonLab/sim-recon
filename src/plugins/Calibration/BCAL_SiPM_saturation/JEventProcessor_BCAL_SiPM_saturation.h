// $Id$
//
//    File: JEventProcessor_BCAL_SiPM_saturation.h
//          Modified file from BCAL_neutron_discriminator.h   ES 5/10/2018
// Created: Thu Apr  5 16:36:00 EDT 2018
// Creator: dalton (on Linux gluon119.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_SiPM_saturation_
#define _JEventProcessor_BCAL_SiPM_saturation_

#include <JANA/JEventProcessor.h>

class JEventProcessor_BCAL_SiPM_saturation:public jana::JEventProcessor{
	public:
		JEventProcessor_BCAL_SiPM_saturation();
		~JEventProcessor_BCAL_SiPM_saturation();
		const char* className(void){return "JEventProcessor_BCAL_SiPM_saturation";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		uint32_t VERBOSE;


		vector< vector<double> > attenuation_parameters;    // store attenuation length
};

#endif // _JEventProcessor_BCAL_SiPM_saturation_

