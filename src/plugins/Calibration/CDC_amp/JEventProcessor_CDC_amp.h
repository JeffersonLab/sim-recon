// $Id$
//
//    File: JEventProcessor_CDC_amp.h
// Created: Tue Sep  6 10:13:02 EDT 2016
// Creator: njarvis (on Linux egbert 2.6.32-642.3.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_CDC_amp_
#define _JEventProcessor_CDC_amp_

#include <JANA/JEventProcessor.h>

class JEventProcessor_CDC_amp:public jana::JEventProcessor{
	public:
		JEventProcessor_CDC_amp();
		~JEventProcessor_CDC_amp();
		const char* className(void){return "JEventProcessor_CDC_amp";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_CDC_amp_
