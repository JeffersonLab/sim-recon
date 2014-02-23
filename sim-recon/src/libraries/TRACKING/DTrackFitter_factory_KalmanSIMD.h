// $Id$
//
//    File: DTrackFitter_factory_Kalman_SIMD.h

#ifndef _DTrackFitter_factory_KalmanSIMD_
#define _DTrackFitter_factory_KalmanSIMD_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitterKalmanSIMD.h>

class DTrackFitter_factory_KalmanSIMD:public jana::JFactory<DTrackFitter>{
	public:
		DTrackFitter_factory_KalmanSIMD(){};
		~DTrackFitter_factory_KalmanSIMD(){};
		const char* Tag(void){return "KalmanSIMD";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackFitter object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackFitter *fitter = new DTrackFitterKalmanSIMD(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(fitter);
			
			return NOERROR;
		}
};

#endif // _DTrackFitter_factory_KalmanSIMD_

