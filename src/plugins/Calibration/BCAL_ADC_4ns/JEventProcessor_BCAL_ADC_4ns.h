// $Id$
//
//    File: JEventProcessor_BCAL_ADC_4ns.h
// Created: Fri Jul 21 10:41:38 EDT 2017
// Creator: dalton (on Linux gluon106.jlab.org 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_BCAL_ADC_4ns_
#define _JEventProcessor_BCAL_ADC_4ns_

#include <JANA/JEventProcessor.h>

class JEventProcessor_BCAL_ADC_4ns:public jana::JEventProcessor{
	public:
		JEventProcessor_BCAL_ADC_4ns();
		~JEventProcessor_BCAL_ADC_4ns();
		const char* className(void){return "JEventProcessor_BCAL_ADC_4ns";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_BCAL_ADC_4ns_

