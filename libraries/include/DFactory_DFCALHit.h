// $Id$
//
//    File: DFactory_DFCALHit.h
// Created: Sun Apr  3 10:41:55 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DFCALHit_
#define _DFactory_DFCALHit_

#include "DFactory.h"
#include "DFCALHit.h"


class DFactory_DFCALHit:public DFactory<DFCALHit>{
	public:
		DFactory_DFCALHit(){};
		~DFactory_DFCALHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DFCALHit_

