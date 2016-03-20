// $Id$
//
//    File: Df250EmulatorAlgorithm_factory.h
// Created: Mar 20, 2016
// Creator: mstaib
//

#ifndef _Df250EmulatorAlgorithm_factory_v1_
#define _Df250EmulatorAlgorithm_factory_v1_

#include <JANA/JFactory.h>
#include <DAQ/Df250EmulatorAlgorithm_v1.h>

class Df250EmulatorAlgorithm_factory_v1:public jana::JFactory<Df250EmulatorAlgorithm>{
	public:
		Df250EmulatorAlgorithm_factory_v1(){};
		~Df250EmulatorAlgorithm_factory_v1(){};
		const char* Tag(void){return "v1";}

	private:
		jerror_t evnt(jana::JEventLoop *loop, uint64_t eventnumber){

			// Create single Df250EmulatorAlgorithm object and mark the factory as
			// persistent so it doesn't get deleted every event.
			Df250EmulatorAlgorithm *emulator = new Df250EmulatorAlgorithm_v1(loop);
			SetFactoryFlag(PERSISTANT);
			ClearFactoryFlag(WRITE_TO_OUTPUT);
			_data.push_back(emulator);
			
			return NOERROR;
		}
};

#endif // _Df250EmulatorAlgorithm_factory_v1_

