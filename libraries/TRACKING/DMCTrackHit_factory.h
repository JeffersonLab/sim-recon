// $Id$
//
//    File: DMCTrackHit_factory.h
// Created: Mon Apr  4 08:18:07 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DMCTrackHit_factory_
#define _DMCTrackHit_factory_

#include "JANA/JFactory.h"
#include "HDDM/hddm_s.h"
#include "DMCTrackHit.h"

class DMCTrackHit_factory:public JFactory<DMCTrackHit>{
	public:
		DMCTrackHit_factory(){};
		~DMCTrackHit_factory(){};
		const string toString(void);

	private:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DMCTrackHit_factory_

