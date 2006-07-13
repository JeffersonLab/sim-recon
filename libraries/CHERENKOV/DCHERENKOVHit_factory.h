// $Id$
//
//    File: DCHERENKOVHit_factory.h
// Created: Thu Jun  9 10:32:49 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DCHERENKOVHit_factory_
#define _DCHERENKOVHit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"

#include "HDDM/hddm_s.h"
#include "DCHERENKOVHit.h"

class DCHERENKOVHit_factory:public JFactory<DCHERENKOVHit>{
	public:
		DCHERENKOVHit_factory(){};
		~DCHERENKOVHit_factory(){};
		jerror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);
	
	protected:
		//jerror_t init(void);						///< Called once at program start.
		//jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DCHERENKOVHit_factory_

