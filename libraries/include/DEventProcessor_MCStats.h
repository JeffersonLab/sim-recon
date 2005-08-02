// $Id$
//
//    File: DEventProcessor_MCStats.h
// Created: Tue Jul 19 09:29:39 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <TH1.h>
#include <TH2.h>

#ifndef _DEventProcessor_MCStats_
#define _DEventProcessor_MCStats_

#include "DEventProcessor.h"

class DEventProcessor_MCStats:public DEventProcessor{
	public:
		DEventProcessor_MCStats(){};
		~DEventProcessor_MCStats(){};
		const char* className(void){return "DEventProcessor_MCStats";}

	private:
		derror_t init(void);						///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		derror_t evnt(DEventLoop *eventLoop, int eventnumber);	///< Called every event.
		derror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		derror_t fini(void);						///< Called after last event of last event source has been processed.


		TH1F *dist_same, *dist_diff, *r0;
		TH1F *stats;
		TH2F *r0_vs_dist_same, *r0_vs_dist_diff;
};

#endif // _DEventProcessor_MCStats_

