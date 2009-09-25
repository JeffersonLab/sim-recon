// $Id$
//
//    File: DTrackHitSelector_factory_ALT1.h
// Created: Fri Feb  6 08:11:38 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackHitSelector_factory_ALT1_
#define _DTrackHitSelector_factory_ALT1_

#include <JANA/JFactory.h>
#include "DTrackHitSelector.h"
#include "DTrackHitSelectorALT1.h"

class DTrackHitSelector_factory_ALT1:public jana::JFactory<DTrackHitSelector>{
	public:
		DTrackHitSelector_factory_ALT1(){};
		~DTrackHitSelector_factory_ALT1(){};
		const char* Tag(void){return "ALT1";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackHitSelector object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackHitSelector *selector = new DTrackHitSelectorALT1(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(selector);
			
			return NOERROR;
		}
};

#endif // _DTrackHitSelector_factory_ALT1_

