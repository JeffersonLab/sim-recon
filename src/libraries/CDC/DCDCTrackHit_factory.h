// $Id$
//
//    File: DCDCTrackHit_factory.h
// Created: Mon Oct 16 10:20:07 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#ifndef _DCDCTrackHit_factory_
#define _DCDCTrackHit_factory_

#include <JANA/JFactory.h>
#include "DCDCTrackHit.h"
#include "HDGEOMETRY/DGeometry.h"

#define CDC_MAX_STRAWS 222
#define CDC_MAX_RINGS 28

/// Provide CDC track hit objects. Currently, these objects are
/// supplied by the simulated data file. Once real data is 
/// available, this factory will be made to take lower-level
/// uncalibrated hit objects and apply calibrations to generate
/// the DCDCHit objects.

class DCDCTrackHit_factory:public JFactory<DCDCTrackHit>{
	public:
		DCDCTrackHit_factory(){};
		~DCDCTrackHit_factory(){};
		
		static const DCDCWire* GetCDCWire(int ring, int straw);

	private:
		jerror_t init(void);
		jerror_t brun(JEventLoop *loop, int runnumber);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method


		DGeometry *dgeom;
		float Z_MIN, Z_MAX;
};

#endif // _DCDCTrackHit_factory_

