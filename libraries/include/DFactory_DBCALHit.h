// $Id$
//
//    File: DFactory_DBCALHit.h
// Created: Sun Apr  3 10:49:22 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DBCALHit_
#define _DFactory_DBCALHit_

#include "DFactory.h"
#include "DBCALHit.h"

class DFactory_DBCALHit:public DFactory<DBCALHit>{
	public:
		DFactory_DBCALHit(){};
		~DFactory_DBCALHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DBCALHit_

