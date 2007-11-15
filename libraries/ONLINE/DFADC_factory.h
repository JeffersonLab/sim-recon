// $Id$
//
//    File: DFADC_factory.h
// Created: Thu Nov 15 09:51:29 EST 2007
// Creator: bellis (on Linux mordor 2.6.22.1 unknown)
//

#ifndef _DFADC_factory_
#define _DFADC_factory_

#include <JANA/JFactory.h>
#include "DFADC.h"

class DFADC_factory:public JFactory<DFADC>{
	public:
		DFADC_factory(){};
		~DFADC_factory(){};
		const string toString(void);


	private:
		//jerror_t init(void);						///< Called once at program start.
		//jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DFADC_factory_

