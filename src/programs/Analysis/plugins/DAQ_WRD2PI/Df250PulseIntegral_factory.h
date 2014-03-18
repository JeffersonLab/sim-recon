// $Id$
//
//    File: Df250PulseIntegral_factory.h
// Created: Thu Feb 13 12:49:12 EST 2014
// Creator: dalton (on Linux gluon104.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _Df250PulseIntegral_factory_
#define _Df250PulseIntegral_factory_

#include <JANA/JFactory.h>

#include <DAQ/Df250PulseIntegral.h>

class Df250PulseIntegral_factory:public jana::JFactory<Df250PulseIntegral>{
	public:
                Df250PulseIntegral_factory(){
		  use_factory=1;
		  ped_samples=5;
		};
		~Df250PulseIntegral_factory(){};
		/* uint32_t mypulse_number; */
		/* uint32_t myquality_factor; */
		/* uint32_t myintegral; */
		/* uint32_t pedestalsum; */
		uint32_t ped_samples;
		/* uint32_t nsamples; */

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _Df250PulseIntegral_factory_

