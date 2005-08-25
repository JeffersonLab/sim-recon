// $Id$
//
//    File: DFactory_DTrackHit.h
// Created: Tue Aug 23 05:00:03 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DTrackHit_
#define _DFactory_DTrackHit_

#include "DFactory.h"
#include "DTrackHit.h"

class DFactory_DTrackHit:public DFactory<DTrackHit>{
	public:
		DFactory_DTrackHit(){};
		~DFactory_DTrackHit(){};
		const string toString(void);


	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DTrackHit_

