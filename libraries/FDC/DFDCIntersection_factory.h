// $Id$
//
//    File: DFDCIntersection_factory.h
// Created: Tue Oct 30 11:24:53 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DFDCIntersection_factory_
#define _DFDCIntersection_factory_

#include <JANA/JFactory.h>
#include "DFDCIntersection.h"

class DFDCIntersection_factory:public JFactory<DFDCIntersection>{
	public:
		DFDCIntersection_factory(){};
		~DFDCIntersection_factory(){};
		const string toString(void);


	private:
		//jerror_t init(void);						///< Called once at program start.
		//jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void MakeIntersectionPoints(vector<vector<const DFDCHit*> >&hits_by_layer);

};

#endif // _DFDCIntersection_factory_

