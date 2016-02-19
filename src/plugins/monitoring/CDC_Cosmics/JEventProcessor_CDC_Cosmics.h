// $Id$
//
//    File: JEventProcessor_CDC_Cosmics.h
// Created: Mon Jul  6 13:00:51 EDT 2015
// Creator: mstaib (on Linux egbert 2.6.32-504.16.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_CDC_Cosmics_
#define _JEventProcessor_CDC_Cosmics_

#include <JANA/JEventProcessor.h>

class JEventProcessor_CDC_Cosmics:public jana::JEventProcessor{
	public:
		JEventProcessor_CDC_Cosmics();
		~JEventProcessor_CDC_Cosmics();
		const char* className(void){return "JEventProcessor_CDC_Cosmics";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
        int EXCLUDERING;
        vector<vector<double> >max_sag;
        vector<vector<double> >sag_phi_offset;
};

#endif // _JEventProcessor_CDC_Cosmics_

