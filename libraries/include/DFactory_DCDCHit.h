// $Id$
//
//    File: DFactory_DCDCHit.h
// Created: Sun Apr  3 10:46:28 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DCDCHit_
#define _DFactory_DCDCHit_

#include "DFactory.h"
#include "DCDCHit.h"

class DFactory_DCDCHit:public DFactory<DCDCHit>{
	public:
		DFactory_DCDCHit(){};
		~DFactory_DCDCHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DCDCHit_

