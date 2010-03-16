// $Id$
//
//    File: DTrackFitter_factory_Kalman.h
// Created: Wed Jul 29 23:16:48 EDT 2009
// Creator: davidl (on Darwin Amelia.local 8.11.1 i386)
//

#ifndef _DTrackFitter_factory_Kalman_
#define _DTrackFitter_factory_Kalman_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitterKalman.h>

class DTrackFitter_factory_Kalman:public jana::JFactory<DTrackFitter>{
	public:
		DTrackFitter_factory_Kalman(){};
		~DTrackFitter_factory_Kalman(){};
		const char* Tag(void){return "Kalman";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackFitter object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackFitter *fitter = new DTrackFitterKalman(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(fitter);
			
			return NOERROR;
		}
};

#endif // _DTrackFitter_factory_Kalman_

