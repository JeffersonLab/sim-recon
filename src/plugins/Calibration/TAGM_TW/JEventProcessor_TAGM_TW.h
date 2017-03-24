// $Id$
//
//    File: JEventProcessor_TAGM_TW.h
// Created: Thu Aug  6 07:33:06 EDT 2015
// Creator: pooser (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TAGM_TW_
#define _JEventProcessor_TAGM_TW_

#include <JANA/JEventProcessor.h>

#include "TAGGER/DTAGMHit.h"
#include "TAGGER/DTAGMGeometry.h"

class JEventProcessor_TAGM_TW:public jana::JEventProcessor{
	public:
		JEventProcessor_TAGM_TW();
		~JEventProcessor_TAGM_TW();
		const char* className(void){return "JEventProcessor_TAGM_TW";}


	private:
		// For the timewalk
		double tw_c0[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
		double tw_c1[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
		double tw_c2[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
		double thresh[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];
		double P_0[TAGM_MAX_ROW+1][TAGM_MAX_COLUMN+1];

		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TAGM_TW_

