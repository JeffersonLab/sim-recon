#ifndef _Df125EmulatorAlgorithm_v2_
#define _Df125EmulatorAlgorithm_v2_

#include "Df125EmulatorAlgorithm.h"

/////////////////////////////////////////////////////////////////
// This implements the base class for the f125 firmware emulation
// EmulateFirmware needs to be virtually overwritten by the user
////////////////////////////////////////////////////////////////

class Df125EmulatorAlgorithm_v2:public Df125EmulatorAlgorithm{
    public:

        Df125EmulatorAlgorithm_v2(JEventLoop *loop);
        ~Df125EmulatorAlgorithm_v2(){};

        //Only the emulation routines need to be overwritten
        void EmulateFirmware(const Df125WindowRawData*, Df125CDCPulse*, Df125FDCPulse*);

    protected:

        Df125EmulatorAlgorithm_v2(){};
        // Many helper functions from the old fa125algos files
        void cdc_hit(Int_t&, Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);   // look for a hit
        void cdc_time(Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t); // find hit time
        void cdc_integral(Long_t&, Int_t&, Int_t, Int_t[], Int_t, Int_t); // find integral
        void cdc_max(Int_t&, Int_t, Int_t[], Int_t); // find first max amplitude after hit
        void upsamplei(Int_t[], Int_t, Int_t[], Int_t);   // upsample
        void fa125_algos(Int_t&, Int_t&, Int_t&, Long_t&, Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);

};

#endif // _Df125EmulatorAlgorithm_v2_
