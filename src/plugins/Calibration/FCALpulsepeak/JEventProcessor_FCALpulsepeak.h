// $Id$
//
//    File: JEventProcessor_FCALpulsepeak.h
// Created: Tue Sep 27 11:18:28 EDT 2016
// Creator: asubedi (on Linux stanley.physics.indiana.edu 2.6.32-573.18.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_FCALpulsepeak_
#define _JEventProcessor_FCALpulsepeak_

#include <JANA/JEventProcessor.h>

class JEventProcessor_FCALpulsepeak:public jana::JEventProcessor{
	public:
		JEventProcessor_FCALpulsepeak();
		~JEventProcessor_FCALpulsepeak();
		const char* className(void){return "JEventProcessor_FCALpulsepeak";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
		int m_x;
  		int m_y;
  		int m_chan;
		double m_peak;


};

#endif // _JEventProcessor_FCALpulsepeak_

