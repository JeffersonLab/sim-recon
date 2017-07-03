// $Id$
//
//    File: JEventProcessor_eventWriterHDDM.h
// Created: Mon Mar  6 14:11:50 EST 2017
// Creator: tbritton (on Linux halld03.jlab.org 3.10.0-514.6.1.el7.x86_64 x86_64)
//

#ifndef _JEventProcessor_eventWriterHDDM_
#define _JEventProcessor_eventWriterHDDM_

#include <JANA/JEventProcessor.h>
#include <HDDM/DEventWriterHDDM.h>

class JEventProcessor_eventWriterHDDM:public jana::JEventProcessor{
	public:
		JEventProcessor_eventWriterHDDM();
		~JEventProcessor_eventWriterHDDM();
		const char* className(void){return "JEventProcessor_eventWriterHDDM";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_eventWriterHDDM_

