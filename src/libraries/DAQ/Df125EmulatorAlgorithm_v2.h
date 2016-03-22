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
        void cdc_hit(Int_t&, Int_t&, Int_t&, const uint16_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);   // look for a hit
        void cdc_time(Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t); // find hit time
        void cdc_integral(Long_t&, Int_t&, Int_t, const uint16_t[], Int_t, Int_t); // find integral
        void cdc_max(Int_t&, Int_t, const uint16_t[], Int_t); // find first max amplitude after hit
        void upsamplei(Int_t[], Int_t, Int_t[], Int_t);   // upsample
        void fa125_algos(Int_t&, Int_t&, Int_t&, Long_t&, Int_t&, Int_t&, const uint16_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);

        // Enables forced use of default values
        int FORCE_DEFAULT_CDC;
        int FORCE_DEFAULT_FDC;

        // Default values for the essential parameters
        Int_t CDC_WS_DEF ;      // hit window start - must be >= F125_CDC_NP
        Int_t CDC_WE_DEF ;      // hit window end - must be at least 20 less than number of samples available
        Int_t CDC_IE_DEF ;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH
        Int_t CDC_NP_DEF ;      // # samples used for pedestal used to find hit. 2**integer
        Int_t CDC_NP2_DEF;      // # samples used for pedestal calculated just before hit. 2**integer
        Int_t CDC_PG_DEF ;      // # samples between hit threshold crossing and local pedestal sample
        Int_t CDC_H_DEF  ;      // 5 sigma hit threshold
        Int_t CDC_TH_DEF ;      // 4 sigma high timing threshold
        Int_t CDC_TL_DEF ;      // 1 sigma low timing threshold

        Int_t FDC_WS_DEF ;      // hit window start - must be >= F125_FDC_NP
        Int_t FDC_WE_DEF ;      // hit window end - must be at least 20 less than number of samples available
        Int_t FDC_IE_DEF ;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH
        Int_t FDC_NP_DEF ;      // # samples used for pedestal used to find hit. 2**integer
        Int_t FDC_NP2_DEF;      // # samples used for pedestal calculated just before hit. 2**integer
        Int_t FDC_PG_DEF ;      // # samples between hit threshold crossing and local pedestal sample
        Int_t FDC_H_DEF  ;      // 5 sigma hit threshold
        Int_t FDC_TH_DEF ;      // 4 sigma high timing threshold
        Int_t FDC_TL_DEF ;      // 1 sigma low timing threshold

        // Set verbosity
        int VERBOSE;
};

#endif // _Df125EmulatorAlgorithm_v2_
