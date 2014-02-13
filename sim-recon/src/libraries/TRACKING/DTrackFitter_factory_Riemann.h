// $Id$
//
//    File: DTrackFitter_factory_Riemann.h

#ifndef _DTrackFitter_factory_Riemann_
#define _DTrackFitter_factory_Riemann_

#include <JANA/JFactory.h>
#include <TRACKING/DTrackFitterRiemann.h>

class DTrackFitter_factory_Riemann:public jana::JFactory<DTrackFitter>{
	public:
		DTrackFitter_factory_Riemann(){};
		~DTrackFitter_factory_Riemann(){};
		const char* Tag(void){return "Riemann";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, int eventnumber){

			// Create single DTrackFitter object and mark the factory as
			// persistent so it doesn't get deleted every event.
			DTrackFitter *fitter = new DTrackFitterRiemann(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(fitter);
			
			return NOERROR;
		}
};

#endif // _DTrackFitter_factory_Riemann_

