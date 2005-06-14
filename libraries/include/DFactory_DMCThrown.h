// $Id$
//
//    File: DFactory_DMCThrown.h
// Created: Sun Apr  3 12:22:09 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DMCThrown_
#define _DFactory_DMCThrown_

#include "DFactory.h"
#include "DMCThrown.h"
#include "DEventLoop.h"

class DFactory_DMCThrown:public DFactory<DMCThrown>{
	public:
		DFactory_DMCThrown(){};
		~DFactory_DMCThrown(){};
		derror_t Extract_HDDM(s_HDDM_t *hddm_s, vector<void*> &v);
		const string toString(void);
	
	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DMCThrown_

