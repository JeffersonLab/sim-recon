// $Id$
//
//    File: DTrackEfficiency_factory.h
// Created: Sun Apr 24 06:45:21 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#ifndef _DTrackEfficiency_factory_
#define _DTrackEfficiency_factory_


#include "JANA/JFactory.h"
#include "DTrackEfficiency.h"

class DTrackEfficiency_factory:public JFactory<DTrackEfficiency>{

	public:
		DTrackEfficiency_factory(){};
		~DTrackEfficiency_factory(){};
		const string toString(void);
	
	private:
		jerror_t evnt(JEventLoop *loop, int eventnumber);	///< Invoked via JEventProcessor virtual method

};

#endif // _DTrackEfficiency_factory_

