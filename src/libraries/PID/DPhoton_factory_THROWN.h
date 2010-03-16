// $Id$
//
//    File: DPhoton_factory_THROWN.h
// Created: Tue Mar  9 22:34:29 EST 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#ifndef _DPhoton_factory_THROWN_
#define _DPhoton_factory_THROWN_

#include <JANA/JFactory.h>
#include "DPhoton.h"

class DPhoton_factory_THROWN:public jana::JFactory<DPhoton>{
	public:
		DPhoton_factory_THROWN(){};
		~DPhoton_factory_THROWN(){};
		const char* Tag(void){return "THROWN";}

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPhoton_factory_THROWN_

