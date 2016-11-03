#ifndef _Df250EmulatorAlgorithm_
#define _Df250EmulatorAlgorithm_
#include <JANA/JObject.h>
#include <JANA/JFactory.h>
#include <JANA/JEventLoop.h>

#include <DAQ/Df250WindowRawData.h>
#include <DAQ/Df250PulseTime.h>
#include <DAQ/Df250PulseData.h>
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
        
        // firmware v1 data format
        virtual void EmulateFirmware(const Df250WindowRawData* wrd,
                                     std::vector<Df250PulseTime*> &pt_objs,
                                     std::vector<Df250PulsePedestal*> &pp_objs,
                                     std::vector<Df250PulseIntegral*> &pi_objs)=0;

        virtual void EmulateFirmware(const Df250WindowRawData* rawData,
                                     std::vector<JObject*> &pt_objs,
                                     std::vector<JObject*> &pp_objs,
                                     std::vector<JObject*> &pi_objs)
			{
				std::vector<Df250PulseTime*> mypt_objs;
				std::vector<Df250PulsePedestal*> mypp_objs;
				std::vector<Df250PulseIntegral*> mypi_objs;
				EmulateFirmware(rawData, mypt_objs, mypp_objs, mypi_objs);
				for(auto p : mypt_objs) pt_objs.push_back(p);
				for(auto p : mypp_objs) pp_objs.push_back(p);
				for(auto p : mypi_objs) pi_objs.push_back(p);

			}

        // firmware v2 data format
        virtual void EmulateFirmware(const Df250WindowRawData* wrd,
                                     std::vector<Df250PulseData*> &pdat_objs)=0;

        virtual void EmulateFirmware(const Df250WindowRawData* rawData,
                                     std::vector<JObject*> &pdat_objs)
			{
				std::vector<Df250PulseData*> mypdat_objs;
				EmulateFirmware(rawData, mypdat_objs);
				for(auto p : mypdat_objs) pdat_objs.push_back(p);

			}
    protected:
        // Suppress default constructor
        Df250EmulatorAlgorithm(){};

};

#endif // _Df250EmulatorAlgorithm_factory_
