// $Id$
//
//    File: JEventProcessor_EventTagPi0.h
// Created: Fri Feb  5 23:23:22 EST 2016
// Creator: davidl (on Darwin harriet 13.4.0 i386)
//

#ifndef _JEventProcessor_EventTagPi0_
#define _JEventProcessor_EventTagPi0_

#include <JANA/JEventProcessor.h>

class JEventProcessor_EventTagPi0:public jana::JEventProcessor{
	public:
		JEventProcessor_EventTagPi0();
		~JEventProcessor_EventTagPi0();
		const char* className(void){return "JEventProcessor_EventTagPi0";}

		double Emin_MeV;
		double Rmin_cm;
		double Rmin_cm_2;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_EventTagPi0_

