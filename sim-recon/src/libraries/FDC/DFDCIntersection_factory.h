// $Id$
//
//    File: DFDCIntersection_factory.h
// Created: Tue Oct 30 11:24:53 EDT 2007
// Creator: davidl (on Darwin fwing-dhcp95.jlab.org 8.10.1 i386)
//

#ifndef _DFDCIntersection_factory_
#define _DFDCIntersection_factory_

#include <JANA/JFactory.h>
#include <DVector2.h>

#include "DFDCIntersection.h"
#include "DFDCGeometry.h"

class DFDCIntersection_factory:public JFactory<DFDCIntersection>{
	public:
		DFDCIntersection_factory(){};
		~DFDCIntersection_factory(){};

	private:
		jerror_t init(void);						///< Called once at program start.
		jerror_t brun(JEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Called every event.
		//jerror_t erun(void);						///< Called everytime run number changes, provided brun has been called.
		//jerror_t fini(void);						///< Called after last event of last event source has been processed.

		void MakeIntersectionPoints(vector<vector<const DFDCHit*> >&hits_by_layer);
		void MakeRestrictedIntersectionPoints(vector<vector<const DFDCHit*> >&hits_by_layer);
		
		void FindIntersections(vector<const DFDCHit*> &layer1, vector<const DFDCHit*> &layer2, vector<DFDCIntersection*> &intersections);

		vector<vector<vector<const DFDCHit*> > > fdchits_by_package; ///< fdchits_by_package[package][layer][hit]
		double MAX_DIST2;

	vector<vector<DFDCWire*> >fdcwires;
	bool USE_FDC;

};

#endif // _DFDCIntersection_factory_

