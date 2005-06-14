// $Id$
//
//    File: DFactory_DBCALHit.h
// Created: Thu Jun  9 10:14:35 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFactory_DBCALHit_
#define _DFactory_DBCALHit_

#include "DFactory.h"
#include "DEventLoop.h"
#include "DBCALHit.h"

class DFactory_DBCALHit:public DFactory<DBCALHit>{
	public:
		DFactory_DBCALHit(){};
		~DFactory_DBCALHit(){};
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);
	
	protected:
		//derror_t init(void);						///< Called once at program start.
		//derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//derror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DFactory_DBCALHit_

