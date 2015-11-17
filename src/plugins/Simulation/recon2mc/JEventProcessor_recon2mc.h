// $Id$
//
//    File: JEventProcessor_recon2mc.h
// Created: Tue Nov 10 13:07:57 EST 2015
// Creator: davidl (on Linux gluon47.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_recon2mc_
#define _JEventProcessor_recon2mc_

#include <JANA/JEventProcessor.h>
#include <HDDM/hddm_s.hpp>
#include <fstream>

using namespace std;

class JEventProcessor_recon2mc:public jana::JEventProcessor{
	public:
		JEventProcessor_recon2mc();
		~JEventProcessor_recon2mc();
		const char* className(void){return "JEventProcessor_recon2mc";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.


		int runNumber;
		hddm_s::ostream *ostr_s;
		ofstream ofs;
		pthread_mutex_t mutex;
};

#endif // _JEventProcessor_recon2mc_

