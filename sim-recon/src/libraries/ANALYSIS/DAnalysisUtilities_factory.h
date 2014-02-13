// $Id$
//
//    File: DAnalysisUtilities_factory.h
// Created: Mon Feb 28 14:12:16 EST 2011
// Creator: pmatt (on Linux ifarml1 2.6.18-128.el5 x86_64)
//

#ifndef _DAnalysisUtilities_factory_
#define _DAnalysisUtilities_factory_

#include <JANA/JFactory.h>
#include "DAnalysisUtilities.h"

class DAnalysisUtilities_factory : public jana::JFactory<DAnalysisUtilities>
{
	public:
		DAnalysisUtilities_factory(){};
		~DAnalysisUtilities_factory(){};
  
	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber)
		{
			// Create single DAnalysisUtilities object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DAnalysisUtilities *locAnalysisUtilities = new DAnalysisUtilities(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(locAnalysisUtilities);
			return NOERROR;
		}
};

#endif // _DAnalysisUtilities_factory_

