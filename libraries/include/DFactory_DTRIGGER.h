// $Id$
//
//    File: DFactory_DTRIGGER.h
// Created: Mon Apr  4 21:45:02 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DTRIGGER_
#define _DFactory_DTRIGGER_

#include "DFactory.h"
#include "DTRIGGER.h"

class DFactory_DTRIGGER:public DFactory<DTRIGGER>{
	public:
		DFactory_DTRIGGER(){};
		~DFactory_DTRIGGER(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DTRIGGER_

