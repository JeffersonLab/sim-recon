// $Id$
//
//    File: DCCALGeometry_factory.h
// Created: Tue Nov 30 15:42:41 EST 2010
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//

#ifndef _DCCALGeometry_factory_
#define _DCCALGeometry_factory_

#include <JANA/JFactory.h>
#include "DCCALGeometry.h"

class DCCALGeometry_factory:public jana::JFactory<DCCALGeometry>{
	public:
		DCCALGeometry_factory(){};
		~DCCALGeometry_factory(){};


	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(jana::JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(jana::JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		jerror_t fini(void);						///< Called after last event of last event source has been processed.
};

#endif // _DCCALGeometry_factory_

