// $Id$
//
//    File: DFactory_DTAGGERHit.h
// Created: Sun Apr  3 10:29:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DTAGGERHit_
#define _DFactory_DTAGGERHit_

#include "DFactory.h"
#include "DTAGGERHit.h"

class DFactory_DTAGGERHit:public DFactory<DTAGGERHit>{
	public:
		DFactory_DTAGGERHit(){};
		~DFactory_DTAGGERHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DTAGGERHit_

