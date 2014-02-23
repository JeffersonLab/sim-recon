// $Id$
//
//    File: JEventProcessor_extract_ptype_hddm.h
// Created: Mon Sep  5 12:29:45 EDT 2011
// Creator: davidl (on Linux ifarm1101 2.6.18-128.7.1.el5 x86_64)
//

#ifndef _JEventProcessor_extract_ptype_hddm_
#define _JEventProcessor_extract_ptype_hddm_

#include <stdlib.h>

#include <JANA/JEventProcessor.h>
#include <HDDM/hddm_s.h>

class JEventProcessor_extract_ptype_hddm:public jana::JEventProcessor{
	public:
		JEventProcessor_extract_ptype_hddm();
		~JEventProcessor_extract_ptype_hddm();
		const char* className(void){return "JEventProcessor_extract_ptype_hddm";}

		double randm(double low, double high){return ((high - low) * drand48() + low);}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.

		pthread_mutex_t mutex;
		s_iostream_t *hddmout;
		string OUTFILENAME;
		unsigned long Nevents;
		unsigned int PTYPE;
};

#endif // _JEventProcessor_extract_ptype_hddm_

