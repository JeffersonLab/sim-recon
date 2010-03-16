// $Id$
//
//    File: DChargedTruthMatch_factory.h
// Created: Sun Jan 31 08:45:38 EST 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DChargedTruthMatch_factory_
#define _DChargedTruthMatch_factory_

#include <JANA/JFactory.h>
#include "DChargedTruthMatch.h"

class DChargedTruthMatch_factory:public jana::JFactory<DChargedTruthMatch>{
	public:
		DChargedTruthMatch_factory(){};
		~DChargedTruthMatch_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DChargedTruthMatch_factory_

