// $Id$
//
//    File: JEventProcessor_syncskim.h
// Created: Wed Feb 22 20:04:25 EST 2017
// Creator: davidl (on Linux gluon48.jlab.org 2.6.32-431.20.3.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_syncskim_
#define _JEventProcessor_syncskim_

#include <JANA/JEventProcessor.h>

#include <TTree.h>
#include <SyncEvent.h>

class JEventProcessor_syncskim:public jana::JEventProcessor{
	public:
		JEventProcessor_syncskim();
		~JEventProcessor_syncskim();
		const char* className(void){return "JEventProcessor_syncskim";}

		SyncEvent synevt;
		TTree *tree;

		// Values to do linear regression to find slope and intercept correlating
		// 250MHz clock time to unix timestamp
		double sum_n;
		double sum_x;
		double sum_y;
		double sum_xy;
		double sum_x2;

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_syncskim_

