// $Id$
//
//    File: Df125EmulatorAlgorithm_factory.h
// Created: Mar 20, 2016
// Creator: mstaib
//

#ifndef _Df125EmulatorAlgorithm_factory_v2_
#define _Df125EmulatorAlgorithm_factory_v2_

#include <JANA/JFactory.h>
#include <DAQ/Df125EmulatorAlgorithm_v2.h>

class Df125EmulatorAlgorithm_factory_v2:public jana::JFactory<Df125EmulatorAlgorithm>{
	public:
		Df125EmulatorAlgorithm_factory_v2(){};
		~Df125EmulatorAlgorithm_factory_v2(){};
		const char* Tag(void){return "v2";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, uint64_t eventnumber){

			// Create single Df125EmulatorAlgorithm object and mark the factory as
			// persistent so it doesn't get deleted every event.
			Df125EmulatorAlgorithm *emulator = new Df125EmulatorAlgorithm_v2(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(emulator);
			
			return NOERROR;
		}
};

#endif // _Df125EmulatorAlgorithm_factory_v2_

