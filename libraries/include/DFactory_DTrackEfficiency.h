// $Id$
//
//    File: DFactory_DTrackEfficiency.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DTrackEfficiency_
#define _DFactory_DTrackEfficiency_


#include "DFactory.h"
#include "DTrackEfficiency.h"

class DFactory_DTrackEfficiency:public DFactory<DTrackEfficiency>{

	public:
		DFactory_DTrackEfficiency(){};
		~DFactory_DTrackEfficiency(){};
		const string toString(void);
	
	private:
		derror_t evnt(DEventLoop *loop, int eventnumber);	///< Invoked via DEventProcessor virtual method

};

#endif // _DFactory_DTrackEfficiency_

