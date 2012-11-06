// $Id$
//
//    File: DPiMinus_factory.h
// Created: Sat Apr 14 12:13:05 EDT 2012
// Creator: davidl (on Darwin genmacbook.local 11.3.0 i386)
//

#ifndef _DPiMinus_factory_
#define _DPiMinus_factory_

#include <JANA/JFactory.h>
#include "DPiMinus.h"

class DPiMinus_factory:public jana::JFactory<DPiMinus>{
	public:
		DPiMinus_factory(){};
		~DPiMinus_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DPiMinus_factory_

