// $Id$
//
//    File: Df125EmulatorAlgorithm_factory.h
// Created: Mar 20, 2016
// Creator: mstaib
//

#ifndef _Df125EmulatorAlgorithm_factory_
#define _Df125EmulatorAlgorithm_factory_

#include <JANA/JFactory.h>
#include <DAQ/Df125EmulatorAlgorithm.h>

class Df125EmulatorAlgorithm_factory:public jana::JFactory<Df125EmulatorAlgorithm>{
	public:
		Df125EmulatorAlgorithm_factory(){};
		~Df125EmulatorAlgorithm_factory(){};

	private:
		jerror_t evnt(jana::JEventLoop *loop, uint64_t eventnumber){

            // This is a trivial class that simply implements the
            // v2 tagged factory as the default. It is here so
            // that the default can be changed easily by simply
            // changing the tag here or on the command line.

            // v1 = ported f250 code (has not been implemented)
            // v2 = firmware using the upsampling technique.

            vector<const Df125EmulatorAlgorithm*> emulators;
            Df125EmulatorAlgorithm_factory *f125EmFac = static_cast<Df125EmulatorAlgorithm_factory*>(loop->GetFactory("Df125EmulatorAlgorithm","v2"));
            if(f125EmFac) f125EmFac->Get(emulators);
            for(unsigned int i=0; i< emulators.size(); i++){
                _data.push_back(const_cast<Df125EmulatorAlgorithm*>(emulators[i]));
            }
            SetFactoryFlag(NOT_OBJECT_OWNER);

            return NOERROR;
        }
};

#endif // _Df125EmulatorAlgorithm_factory_

