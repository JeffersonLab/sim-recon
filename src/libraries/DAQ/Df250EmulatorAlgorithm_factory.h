// $Id$
//
//    File: Df250EmulatorAlgorithm_factory.h
// Created: Mar 20, 2016
// Creator: mstaib
//

#ifndef _Df250EmulatorAlgorithm_factory_
#define _Df250EmulatorAlgorithm_factory_

#include <JANA/JFactory.h>
#include <DAQ/Df250EmulatorAlgorithm.h>

class Df250EmulatorAlgorithm_factory:public jana::JFactory<Df250EmulatorAlgorithm>{
	public:
		Df250EmulatorAlgorithm_factory(){};
		~Df250EmulatorAlgorithm_factory(){};

	private:
		jerror_t evnt(jana::JEventLoop *loop, uint64_t eventnumber){

            // This is a trivial class that simply implements the
            // v1 tagged factory as the default. It is here so
            // that the default can be changed easily by simply
            // changing the tag here or on the command line.

            // v1 = f250 code used until mid 2016
            // v1 = f250 code used starting fall 2016

            vector<const Df250EmulatorAlgorithm*> emulators;
            Df250EmulatorAlgorithm_factory *f250EmFac = static_cast<Df250EmulatorAlgorithm_factory*>(loop->GetFactory("Df250EmulatorAlgorithm","v1"));
            if (f250EmFac) f250EmFac->Get(emulators);
            for(unsigned int i=0; i< emulators.size(); i++){
                _data.push_back(const_cast<Df250EmulatorAlgorithm*>(emulators[i]));
            }

            // add the v2 factory for now as well
            Df250EmulatorAlgorithm_factory *f250EmFacV2 = static_cast<Df250EmulatorAlgorithm_factory*>(loop->GetFactory("Df250EmulatorAlgorithm","v2"));
            if (f250EmFacV2) f250EmFacV2->Get(emulators);
            for(unsigned int i=0; i< emulators.size(); i++){
                _data.push_back(const_cast<Df250EmulatorAlgorithm*>(emulators[i]));
            }
            SetFactoryFlag(NOT_OBJECT_OWNER);

            return NOERROR;
        }
};

#endif // _Df250EmulatorAlgorithm_factory_

