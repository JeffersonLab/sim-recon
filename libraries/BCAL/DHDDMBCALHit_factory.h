// $Id$
//
//    File: DHDDMBCALHit_factory.h
// Created: Thu Jun  9 10:25:22 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DHDDMBCALHit_factory_
#define _DHDDMBCALHit_factory_

#include "JANA/JFactory.h"
#include "JANA/JEventLoop.h"
#include "DHDDMBCALHit.h"

class DHDDMBCALHit_factory:public JFactory<DHDDMBCALHit>{
	public:
		DHDDMBCALHit_factory(){};
		~DHDDMBCALHit_factory(){};	    
		const string toString(void);

        private:
                jerror_t evnt(JEventLoop *loop, int eventnumber);       ///< Invoked via JEventProcessor virtual method
	
};

#endif // _DHDDMBCALHit_factory_

