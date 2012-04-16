// $Id$
//
//    File: DPiPlus_factory.h
// Created: Sat Apr 14 12:13:05 EDT 2012
// Creator: davidl (on Darwin genmacbook.local 11.3.0 i386)
//

#ifndef _DPiPlus_factory_
#define _DPiPlus_factory_

#include <JANA/JFactory.h>
#include "DPiPlus.h"

class DPiPlus_factory:public jana::JFactory<DPiPlus>{
	public:
		DPiPlus_factory(){};
		~DPiPlus_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPiPlus_factory_

