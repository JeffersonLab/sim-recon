// $Id$
//
//    File: DTrackFitter_factory_Kalman_SIMD_ALT1.h

#ifndef _DTrackFitter_factory_KalmanSIMD_ALT1_
#define _DTrackFitter_factory_KalmanSIMD_ALT1_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitterKalmanSIMD_ALT1.h>

class DTrackFitter_factory_KalmanSIMD_ALT1:public jana::JFactory<DTrackFitter>{
	public:
		DTrackFitter_factory_KalmanSIMD_ALT1(){};
		~DTrackFitter_factory_KalmanSIMD_ALT1(){};
		const char* Tag(void){return "KalmanSIMD_ALT1";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackFitter object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackFitter *fitter = new DTrackFitterKalmanSIMD_ALT1(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(fitter);
			
			return NOERROR;
		}
};

#endif // _DTrackFitter_factory_KalmanSIMD_ALT1_

