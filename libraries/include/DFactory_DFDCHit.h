// $Id$
//
//    File: DFactory_DFDCHit.h
// Created: Sun Apr  3 10:39:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DFDCHit_
#define _DFactory_DFDCHit_

#include "DFactory.h"
#include "DFDCHit.h"

class DFactory_DFDCHit:public DFactory<DFDCHit>{
	public:
		DFactory_DFDCHit(){};
		~DFactory_DFDCHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DFDCHit_

