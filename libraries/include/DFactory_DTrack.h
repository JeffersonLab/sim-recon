// $Id$
//
//    File: DFactory_DMCReconstructed.h
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DMCReconstructed_
#define _DFactory_DMCReconstructed_

#include "DFactory.h"
#include "DMCReconstructed.h"

class DFactory_DMCReconstructed:public DFactory<DMCReconstructed>{
	public:
		DFactory_DMCReconstructed(){};
		~DFactory_DMCReconstructed(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
};

#endif // _DFactory_DMCReconstructed_

