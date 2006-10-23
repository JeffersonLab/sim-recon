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

#define CDC_MAX_STRAWS 228
#define CDC_MAX_RINGS 23

class DCDCTrackHit_factory:public JFactory<DCDCTrackHit>{
	public:
		DCDCTrackHit_factory(){};
		~DCDCTrackHit_factory(){};
		const string toString(void);

	private:
		jerror_t init(void);
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		float Z_MIN, Z_MAX;
};

#endif // _DCDCTrackHit_factory_

