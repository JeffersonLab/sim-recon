// $Id$
//
//    File: DFactory_DHDDMBCALHit.h
// Created: Thu Jun  9 10:25:22 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DFactory_DHDDMBCALHit_
#define _DFactory_DHDDMBCALHit_

#include "DFactory.h"
#include "DEventLoop.h"
#include "DHDDMBCALHit.h"

class DFactory_DHDDMBCALHit:public DFactory<DHDDMBCALHit>{
	public:
		DFactory_DHDDMBCALHit(){};
		~DFactory_DHDDMBCALHit(){};	    
		const string toString(void);

        private:
                derror_t evnt(DEventLoop *loop, int eventnumber);       ///< Invoked via DEventProcessor virtual method
	
};

#endif // _DFactory_DHDDMBCALHit_

