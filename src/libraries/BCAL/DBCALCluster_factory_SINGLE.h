// $Id$
//
//    File: DBCALCluster_factory_SINGLE.h
// Created: Fri Sep  7 12:13:07 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#ifndef _DBCALCluster_factory_SINGLE_
#define _DBCALCluster_factory_SINGLE_

#include <JANA/JFactory.h>
#include "DBCALCluster.h"

class DBCALCluster_factory_SINGLE:public jana::JFactory<DBCALCluster>{

	/// This factory will create a single DBCALCluster objects from
	/// all of the DBCALPoint objects. It is intended only for
	/// debugging of simulated data where a single cluster is expected.
	///
	/// If there are no DBCALPoint objects, then the DBCALCluster
	/// object is not created.

	public:
		DBCALCluster_factory_SINGLE(){};
		~DBCALCluster_factory_SINGLE(){};
		const char* Tag(void){return "SINGLE";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DBCALCluster_factory_SINGLE_

