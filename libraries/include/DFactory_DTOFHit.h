// $Id$
//
//    File: DFactory_DTOFHit.h
// Created: Sun Apr  3 10:31:26 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DTOFHit_
#define _DFactory_DTOFHit_

#include "DFactory.h"
#include "DTOFHit.h"

class DFactory_DTOFHit:public DFactory<DTOFHit>{
	public:
		DFactory_DTOFHit(){};
		~DFactory_DTOFHit(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DTOFHit_

