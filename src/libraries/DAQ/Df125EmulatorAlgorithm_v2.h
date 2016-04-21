#ifndef _Df125EmulatorAlgorithm_v2_
#define _Df125EmulatorAlgorithm_v2_

#include "Df125EmulatorAlgorithm.h"

/////////////////////////////////////////////////////////////////
// This implements the base class for the f125 firmware emulation
// EmulateFirmware needs to be virtually overwritten by the user
////////////////////////////////////////////////////////////////

class Df125EmulatorAlgorithm_v2:public Df125EmulatorAlgorithm{
    public:

        Df125EmulatorAlgorithm_v2();
        ~Df125EmulatorAlgorithm_v2(){};

        //Only the emulation routines need to be overwritten
        void EmulateFirmware(const Df125WindowRawData*, Df125CDCPulse*, Df125FDCPulse*);

        // Many helper functions from the old fa125algos files
        void fa125_hit(Int_t&, Int_t&, Int_t&, const uint16_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);   // look for a hit
        void fa125_time(Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t); // find hit time
        void fa125_integral(Long_t&, Int_t&, Int_t, const uint16_t[], Int_t, Int_t); // find integral
        void fa125_max(Int_t&, Int_t&, Int_t, const uint16_t[], Int_t); // find first max amplitude after hit
        void fa125_algos(Int_t&, Int_t&, Int_t&, Long_t&, Int_t&, Int_t&, Int_t&, const uint16_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);


    protected:
        //Df125EmulatorAlgorithm_v2(){};

        void upsamplei(Int_t[], Int_t, Int_t[], Int_t);   // upsample

   private:

        // Enables forced use of default values
        int FORCE_DEFAULT_CDC;
        int FORCE_DEFAULT_FDC;

        // Default values, used if FORCE_DEFAULT_xDC=1 or if BORConfig is not found
        Int_t CDC_WS_DEF ;      // hit window start - must be >= F125_CDC_NP
        Int_t CDC_WE_DEF ;      // hit window end - must be at least 20 less than number of samples available
        Int_t CDC_IE_DEF ;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH
        Int_t CDC_P1_DEF ;      // 2**P1 = # samples used for initial pedestal, used to find hit
        Int_t CDC_P2_DEF ;      // 2**P2 = # samples used for pedestal calculated just before hit
        Int_t CDC_PG_DEF ;      // # samples between hit threshold crossing and local pedestal sample
        Int_t CDC_H_DEF  ;      // 5 sigma hit threshold
        Int_t CDC_TH_DEF ;      // 4 sigma high timing threshold
        Int_t CDC_TL_DEF ;      // 1 sigma low timing threshold
        Int_t CDC_IBIT_DEF ;    // scaling factor for integral
        Int_t CDC_ABIT_DEF ;    // scaling factor for amplitude
        Int_t CDC_PBIT_DEF ;    // scaling factor for pedestal


        Int_t FDC_WS_DEF ;      // hit window start - must be >= F125_FDC_NP
        Int_t FDC_WE_DEF ;      // hit window end - must be at least 20 less than number of samples available
        Int_t FDC_IE_DEF ;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH
        Int_t FDC_P1_DEF ;      // 2**P1 = # samples used for initial pedestal, used to find hit
        Int_t FDC_P2_DEF ;      // 2**P2 = # samples used for pedestal calculated just before hit
        Int_t FDC_PG_DEF ;      // # samples between hit threshold crossing and local pedestal sample
        Int_t FDC_H_DEF  ;      // 5 sigma hit threshold
        Int_t FDC_TH_DEF ;      // 4 sigma high timing threshold
        Int_t FDC_TL_DEF ;      // 1 sigma low timing threshold
        Int_t FDC_IBIT_DEF ;    // scaling factor for integral
        Int_t FDC_ABIT_DEF ;    // scaling factor for amplitude
        Int_t FDC_PBIT_DEF ;    // scaling factor for pedestal


        //Override values, used if given in command line
        Int_t CDC_WS ;      // hit window start - must be >= F125_CDC_NP
        Int_t CDC_WE ;      // hit window end - must be at least 20 less than number of samples available
        Int_t CDC_IE ;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH
        Int_t CDC_P1 ;      // 2**P1 = # samples used for initial pedestal, used to find hit
        Int_t CDC_P2 ;      // 2**P2 = # samples used for pedestal calculated just before hit
        Int_t CDC_PG ;      // # samples between hit threshold crossing and local pedestal sample
        Int_t CDC_H  ;      // 5 sigma hit threshold
        Int_t CDC_TH ;      // 4 sigma high timing threshold
        Int_t CDC_TL ;      // 1 sigma low timing threshold
        Int_t CDC_IBIT ;    // scaling factor for integral
        Int_t CDC_ABIT ;    // scaling factor for amplitude
        Int_t CDC_PBIT ;    // scaling factor for pedestal

        Int_t FDC_WS ;      // hit window start - must be >= F125_FDC_NP
        Int_t FDC_WE ;      // hit window end - must be at least 20 less than number of samples available
        Int_t FDC_IE ;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH
        Int_t FDC_P1 ;      // 2**P1 = # samples used for initial pedestal, used to find hit
        Int_t FDC_P2 ;      // 2**P2 = # samples used for pedestal calculated just before hit
        Int_t FDC_PG ;      // # samples between hit threshold crossing and local pedestal sample
        Int_t FDC_H  ;      // 5 sigma hit threshold
        Int_t FDC_TH ;      // 4 sigma high timing threshold
        Int_t FDC_TL ;      // 1 sigma low timing threshold
        Int_t FDC_IBIT ;    // scaling factor for integral
        Int_t FDC_ABIT ;    // scaling factor for amplitude
        Int_t FDC_PBIT ;    // scaling factor for pedestal

        // Set verbosity
        int VERBOSE;
};

#endif // _Df125EmulatorAlgorithm_v2_



