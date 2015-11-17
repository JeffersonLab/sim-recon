// $Id$
//
//    File: JEventProcessor_CDC_TimeToDistance.h
// Created: Mon Nov  9 12:37:01 EST 2015
// Creator: mstaib (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_CDC_TimeToDistance_
#define _JEventProcessor_CDC_TimeToDistance_

#include <JANA/JEventProcessor.h>

class JEventProcessor_CDC_TimeToDistance:public jana::JEventProcessor{
	public:
		JEventProcessor_CDC_TimeToDistance();
		~JEventProcessor_CDC_TimeToDistance();
		const char* className(void){return "JEventProcessor_CDC_TimeToDistance";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
        vector<vector<double> >max_sag;
        vector<vector<double> >sag_phi_offset;
};

#endif // _JEventProcessor_CDC_TimeToDistance_

