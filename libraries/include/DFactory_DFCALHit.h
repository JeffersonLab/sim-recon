// $Id$
//
//    File: DFactory_DFCALHit.h
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFactory_DFCALHit_
#define _DFactory_DFCALHit_

#include "DFactory.h"
#include "DEventLoop.h"
#include "DFCALHit.h"

class DFactory_DFCALHit:public DFactory<DFCALHit>{
	public:
		DFactory_DFCALHit(){};
		~DFactory_DFCALHit(){};
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);
	
	protected:
		//derror_t init(void);						///< Called once at program start.
		//derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//derror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DFactory_DFCALHit_

