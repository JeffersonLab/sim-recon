#ifndef _Df250EmulatorAlgorithm_
#define _Df250EmulatorAlgorithm_
#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250PulsePedestal.h>
#include <DAQ/Df250PulseIntegral.h>
#include <DAQ/Df250Config.h>
#include <DAQ/Df250BORConfig.h>

using namespace std;
using namespace jana;

/////////////////////////////////////////////////////////////////
// This implements the base class for the f250 firmware emulation
// EmulateFirmware needs to be virtually overwritten by the user
////////////////////////////////////////////////////////////////

class Df250EmulatorAlgorithm:public jana::JObject{
    public:
        JOBJECT_PUBLIC(Df250EmulatorAlgorithm);
        Df250EmulatorAlgorithm(JEventLoop *loop){};
        ~Df250EmulatorAlgorithm(){};

        // The main emulation routines are overwritten in the inherited classes
        virtual void EmulateFirmware(const Df250WindowRawData* wrd,
                                     std::vector<JObject*> &pt_objs,
                                     std::vector<JObject*> &pp_objs,
                                     std::vector<JObject*> &pi_objs) = 0;

    protected:
        // Suppress default constructor
        Df250EmulatorAlgorithm(){};

};

#endif // _Df250EmulatorAlgorithm_factory_
