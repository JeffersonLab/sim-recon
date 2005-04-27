// $Id$
//
//    File: DFactory_DMCFitStats.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DMCFitStats_
#define _DFactory_DMCFitStats_

#include <TH1.h>
#include <TH2.h>

#include "DFactory.h"
#include "DMCFitStats.h"

class DFactory_DMCFitStats:public DFactory<DMCFitStats>{

	/// This factory is unusual in that it produces only a
	/// single object. The object (DMCFitStats) is really a
	/// collection of histograms which are maintained across
	/// all events (the factory has its PERSISTENT flag set).
	/// Pretty much all of the work is done in the DMCFitStats
	/// object itself, or rather, it's AddEvent() method.

	public:
		DFactory_DMCFitStats();
		~DFactory_DMCFitStats();
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method
		derror_t fini(void);					///< Invoked via DEventProcessor virtual method
		DMCFitStats *fitstats;

};

#endif // _DFactory_DMCFitStats_

