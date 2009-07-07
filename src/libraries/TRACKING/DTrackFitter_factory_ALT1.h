// $Id$
//
//    File: DTrackFitter_factory_ALT1.h
// Created: Mon Sep  1 10:29:51 EDT 2008
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DTrackFitter_factory_ALT1_
#define _DTrackFitter_factory_ALT1_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitterALT1.h>

class DTrackFitter_factory_ALT1:public jana::JFactory<DTrackFitter>{
	public:
		DTrackFitter_factory_ALT1(){};
		~DTrackFitter_factory_ALT1(){};
		const char* Tag(void){return "ALT1";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackFitter object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackFitter *fitter = new DTrackFitterALT1(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(fitter);
			
			return NOERROR;
		}
};

#endif // _DTrackFitter_factory_ALT1_

