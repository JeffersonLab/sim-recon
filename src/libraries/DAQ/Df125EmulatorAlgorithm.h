#ifndef _Df125EmulatorAlgorithm_
#define _Df125EmulatorAlgorithm_
#include <JANA/JObject.h>
#include <JANA/JEventLoop.h>

#include <stdint.h>
#include <vector>
#include <iostream>
using namespace std;
using namespace jana;
#include <Rtypes.h>

#include <DAQ/Df125WindowRawData.h>
#include <DAQ/Df125CDCPulse.h>
#include <DAQ/Df125FDCPulse.h>
#include <DAQ/Df125Config.h>
#include <DAQ/Df125BORConfig.h>

/////////////////////////////////////////////////////////////////
// This implements the base class for the f125 firmware emulation
// EmulateFirmware needs to be virtually overwritten by the user
////////////////////////////////////////////////////////////////

class Df125EmulatorAlgorithm:public jana::JObject{
    public:
        JOBJECT_PUBLIC(Df125EmulatorAlgorithm);

        Df125EmulatorAlgorithm(){};
        ~Df125EmulatorAlgorithm(){};

        // The main emulation routines are overwritten in the inherited classes
        virtual void EmulateFirmware(const Df125WindowRawData*, Df125CDCPulse*, Df125FDCPulse*) = 0;

    protected:
	//        Df125EmulatorAlgorithm(){};

};

#endif // _Df125EmulatorAlgorithm_factory_
