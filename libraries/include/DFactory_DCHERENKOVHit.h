// $Id$
//
//    File: DFactory_DCHERENKOVHit.h
// Created: Sun Apr  3 10:45:03 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DCHERENKOVHit_
#define _DFactory_DCHERENKOVHit_

#include "DFactory.h"
#include "DCHERENKOVHit.h"

class DFactory_DCHERENKOVHit:public DFactory<DCHERENKOVHit>{
	public:
		DFactory_DCHERENKOVHit(){};
		~DFactory_DCHERENKOVHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DCHERENKOVHit_

