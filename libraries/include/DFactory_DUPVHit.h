// $Id$
//
//    File: DFactory_DUPVHit.h
// Created: Sun Apr  3 10:35:31 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DUPVHit_
#define _DFactory_DUPVHit_

#include "DFactory.h"
#include "DUPVHit.h"

class DFactory_DUPVHit:public DFactory<DUPVHit>{
	public:
		DFactory_DUPVHit(){};
		~DFactory_DUPVHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DUPVHit_

