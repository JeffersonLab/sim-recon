// $Id$
//
//    File: DVertexTimeBased_factory.h
// Created: Wed Apr  7 10:54:41 EDT 2010
// Creator: davidl (on Darwin harriet.jlab.org 9.8.0 i386)
//

#ifndef _DVertexTimeBased_factory_
#define _DVertexTimeBased_factory_

#include <JANA/JFactory.h>
#include <PID/DVertexTimeBased.h>
#include <PID/DVertexCalculator.h>

class DVertexTimeBased_factory:public jana::JFactory<DVertexTimeBased>, public DVertexCalculator{
	public:
		DVertexTimeBased_factory(){};
		~DVertexTimeBased_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DVertexTimeBased_factory_

