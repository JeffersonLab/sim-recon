// $Id$
//
//    File: DTrack_factory.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrack_factory_
#define _DTrack_factory_

#include "JANA/JFactory.h"
#include "JANA/JGeometry.h"
#include "DMagneticFieldMap.h"
#include "DTrack.h"


class DTrack_factory:public JFactory<DTrack>{
	public:
		DTrack_factory(){};
		~DTrack_factory(){};
		const string toString(void);
	
	private:
		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		const JGeometry *dgeom;
		const DMagneticFieldMap *bfield;
};

#endif // _DTrack_factory_

