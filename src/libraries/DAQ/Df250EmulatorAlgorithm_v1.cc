#include <DAQ/Df250EmulatorAlgorithm_v1.h>

Df250EmulatorAlgorithm_v1::Df250EmulatorAlgorithm_v1(JEventLoop *loop){
    ;
}

void Df250EmulatorAlgorithm_v1::EmulateFirmware(const Df250WindowRawData*, Df250PulseTime*, Df250PulsePedestal*, Df250PulseIntegral*){
   jout << " Df250EmulatorAlgorithm_v1::EmulateFirmware Emulating f250 Firmware" << endl;    
   return;
}
