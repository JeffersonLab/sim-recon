// $Id$
//
//    File: DTrackHit_factory.h
// Created: Tue Aug 23 05:00:03 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackHit_factory_
#define _DTrackHit_factory_

#include "JANA/JFactory.h"
#include "DTrackHit.h"

class DTrackHit_factory:public JFactory<DTrackHit>{
	public:
		DTrackHit_factory(){};
		~DTrackHit_factory(){};
		const string toString(void);

	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DTrackHit_factory_

