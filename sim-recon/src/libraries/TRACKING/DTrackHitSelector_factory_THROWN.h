// $Id$
//
//    File: DTrackHitSelector_factory_THROWN.h
// Created: Mon Mar  9 09:00:38 EDT 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackHitSelector_factory_THROWN_
#define _DTrackHitSelector_factory_THROWN_

#include <JANA/JFactory.h>
#include "DTrackHitSelector.h"
#include "DTrackHitSelectorTHROWN.h"

class DTrackHitSelector_factory_THROWN:public jana::JFactory<DTrackHitSelector>{
	public:
		DTrackHitSelector_factory_THROWN(){};
		~DTrackHitSelector_factory_THROWN(){};
		const char* Tag(void){return "THROWN";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackHitSelector object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackHitSelector *selector = new DTrackHitSelectorTHROWN(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(selector);
			
			return NOERROR;
		}
};

#endif // _DTrackHitSelector_factory_THROWN_

