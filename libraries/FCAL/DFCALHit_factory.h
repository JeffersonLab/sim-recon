// $Id$
//
//    File: DFCALHit_factory.h
// Created: Thu Jun  9 10:29:52 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFCALHit_factory_
#define _DFCALHit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DFCALHit.h"
#include "units.h"

class DFCALHit_factory:public JFactory<DFCALHit>{
	
	public:
	
	    DFCALHit_factory(){};
		~DFCALHit_factory(){};
		const string toString(void);
		
	protected:
		//jerror_t init(void);						///< Called once at program start.
		//jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		//jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DFCALHit_factory_

