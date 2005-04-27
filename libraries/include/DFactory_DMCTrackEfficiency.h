// $Id$
//
//    File: DFactory_DMCTrackEfficiency.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DFactory_DMCTrackEfficiency_
#define _DFactory_DMCTrackEfficiency_


#include "DFactory.h"
#include "DMCTrackEfficiency.h"

class DFactory_DMCTrackEfficiency:public DFactory<DMCTrackEfficiency>{

	public:
		DFactory_DMCTrackEfficiency(){};
		~DFactory_DMCTrackEfficiency(){};
		const string toString(void);
	
	private:
		derror_t evnt(int eventnumber);	///< Invoked via DEventProcessor virtual method

};

#endif // _DFactory_DMCTrackEfficiency_

