// $Id$
//
//    File: DTrackFitter_factory_LSLM.h
// Created: Wed Jan 14 08:59:27 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//

#ifndef _DTrackFitter_factory_LSLM_
#define _DTrackFitter_factory_LSLM_

#include <JANA/JFactory.h>
#include "DTrackLSFitter.h"

class DTrackFitter_factory_LSLM:public jana::JFactory<DTrackFitter>{
	public:
		DTrackFitter_factory_LSLM(){};
		~DTrackFitter_factory_LSLM(){};
		const char* Tag(void){return "LSLM";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackFitter object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackFitter *fitter = new DTrackLSFitter(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(fitter);
			
			return NOERROR;
		}
};

#endif // _DTrackFitter_factory_LSLM_

