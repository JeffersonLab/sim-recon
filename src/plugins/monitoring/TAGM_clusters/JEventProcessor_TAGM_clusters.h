// $Id$
//
//    File: JEventProcessor_TAGM_clusters.h
// Created: Tue Jul  5 21:19:22 EDT 2016
// Creator: barnes (on Linux gluey.phys.uconn.edu 2.6.32-573.22.1.el6.x86_64 x86_64)
//

#ifndef _JEventProcessor_TAGM_clusters_
#define _JEventProcessor_TAGM_clusters_

#include <JANA/JEventProcessor.h>
#include <TAGGER/DTAGMHit.h>

class JEventProcessor_TAGM_clusters:public jana::JEventProcessor{
	public:
		JEventProcessor_TAGM_clusters();
		~JEventProcessor_TAGM_clusters();
		const char* className(void){return "JEventProcessor_TAGM_clusters";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int32_t runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, uint64_t eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _JEventProcessor_TAGM_clusters_

