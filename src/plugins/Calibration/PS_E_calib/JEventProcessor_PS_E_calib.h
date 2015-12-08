// $Id$
//
//    File: JEventProcessor_PS_E_calib.h
// Created: Thu Jul  9 17:44:32 EDT 2015
// Creator: aebarnes (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_PS_E_calib_
#define _JEventProcessor_PS_E_calib_

#include <JANA/JEventProcessor.h>

class JEventProcessor_PS_E_calib:public jana::JEventProcessor{
	public:
		JEventProcessor_PS_E_calib();
		~JEventProcessor_PS_E_calib();
		const char* className(void){return "JEventProcessor_PS_E_calib";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_PS_E_calib_

