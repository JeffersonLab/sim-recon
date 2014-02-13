// $Id$
//
//    File: DTrackHitSelector_factory_ALT2.h
// Created: Wed Jan 19 08:28:53 EST 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.6.0 i386)
//

#ifndef _DTrackHitSelector_factory_ALT2_
#define _DTrackHitSelector_factory_ALT2_

#include <JANA/JFactory.h>
#include "DTrackHitSelector.h"
#include "DTrackHitSelectorALT2.h"

class DTrackHitSelector_factory_ALT2:public jana::JFactory<DTrackHitSelector>{
	public:
		DTrackHitSelector_factory_ALT2(){};
		~DTrackHitSelector_factory_ALT2(){};
		const char* Tag(void){return "ALT2";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackHitSelector object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackHitSelector *selector = new DTrackHitSelectorALT2(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(selector);
			
			return NOERROR;
		}
};

#endif // _DTrackHitSelector_factory_ALT2_

