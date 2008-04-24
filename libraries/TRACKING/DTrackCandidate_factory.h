// $Id$
//
//    File: DTrackCandidate_factory.h
// Created: Mon Jul 18 15:23:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DTrackCandidate_factory_
#define _DTrackCandidate_factory_


#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>
using namespace jana;

#include "DTrackCandidate.h"

class DMagneticFieldMap;

class DTrackCandidate_factory:public JFactory<DTrackCandidate>{
	public:
		DTrackCandidate_factory(){};
		~DTrackCandidate_factory(){};


	protected:
		virtual jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method
};

#endif // _DTrackCandidate_factory_

