// $Id$
//
//    File: DFCALShower_factory.h
// Created: Tue May 17 11:57:50 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8smp i686)
//

#ifndef _DFCALShower_factory_
#define _DFCALShower_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DFCALShower.h"

class DFCALShower_factory:public JFactory<DFCALShower>{
	public:
		DFCALShower_factory();
		~DFCALShower_factory(){};
		const string toString(void);
	
	private:
		jerror_t evnt(JEventLoop *eventLoop, int eventnumber);	///< Invoked via JEventProcessor virtual method

		float MAX_SHOWER_DIST;
};

#endif // _DFCALShower_factory_

