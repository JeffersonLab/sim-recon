#include "Df125EmulatorAlgorithm_v2.h"

// Some masks that are useful for getting data from BORConfig registers
// from fa125Lib.h
/* 0x1058 FE nw register defintions */
#define FA125_FE_NW_MASK          0x000003FF
/* 0x105C FE pl register defintions */
#define FA125_FE_PL_MASK           0x0000FFFF
/* 0xN070 - 0xN084 threshold register defintions */
#define FA125_FE_THRESHOLD_MASK          0x00000FFF
/* 0xN0A0 FE ped_sf definitions */
#define FA125_FE_PED_SF_P1_MASK       0x000000FF
#define FA125_FE_PED_SF_P2_MASK       0x0000FF00
#define FA125_FE_PED_SF_IBIT_MASK     0x00070000
#define FA125_FE_PED_SF_ABIT_MASK     0x00380000
#define FA125_FE_PED_SF_PBIT_MASK     0x01C00000
#define FA125_FE_PED_SF_PBIT_SIGN     (1<<25)
/* 0xN0A4 FE timing_thres definitions */
#define FA125_FE_TIMING_THRES_HI_MASK(x) (0x1FF<<((x%3)*9))
#define FA125_FE_TIMING_THRES_LO_MASK(x) (0xFF<<(8+((x%2)*16)))
/* 0xN0B0 FE integration_end definitions */
#define FA125_FE_IE_INTEGRATION_END_MASK  0x00000FFF
#define FA125_FE_IE_PEDESTAL_GAP_MASK     0x000FF000

Df125EmulatorAlgorithm_v2::Df125EmulatorAlgorithm_v2(){

  // Normal use - BORConfig params are used if found.  If not, default values are used
  // EMULATION125:FORCE_DEFAULT_CDC=1 - use default values for CDC instead of BORConfig params
  // EMULATION125:FORCE_DEFAULT_FDC=1 - use default values for FDC instead of BORConfig params
  // Params in command line override both default and BORConfig
  //
  // eg EMULATION125:CDC_H=200 will cause the emulation to use CDC hit threshold=200, and either the BORConfig values (if present) 
  // or the default emulation parameter values for everything else.
  //
  // EMULATION125:CDC_H=200 used with EVIO:F125_EMULATION_MODE=1 will replace the firmware quantities with those emulated using CDC_H=200.
  //

    // Enables forced use of default values
    FORCE_DEFAULT_CDC = 0;
    FORCE_DEFAULT_FDC = 0;
    
    // Default values for the essential parameters

    // PG or 'pedestal gap' is the number of samples between the end of the initial pedestal window and start of the hit search
    // and also the number of samples between the local pedestal window and the hit threshold crossing sample

    // 2**P2 samples are summed to give the integrated local pedestal
    // The firmware returns the (integrated local pedestal)>>(P2+PBIT)
    // Firmware 2_00_F + accepts a negative value for PBIT, but the driver ensures that 3 >= PBIT >= -P2
    // Using PBIT=-P2 obtains pedestal = P2 summed samples 
    // Using PBIT=0 obtains pedestal = mean of P2 samples


    CDC_WS_DEF  =  16;      // hit window start - must be >= F125_CDC_NP
    CDC_WE_DEF  = 179;      // hit window end - must be at least 20 less than number of samples available, eg WE=179 for 200 samples
    CDC_IE_DEF  = 200;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH  
    CDC_P1_DEF  =   4;      // 2**P1 = # samples used for initial pedestal, used to find hit
    CDC_P2_DEF  =   4;      // 2**P2 = # samples used for pedestal calculated just before hit, and returned 
    CDC_PG_DEF  =   4;      // # samples between hit threshold crossing and local pedestal sample 'pedestal gap'
    CDC_H_DEF   = 125;      // 5 sigma hit threshold
    CDC_TH_DEF  = 100;      // 4 sigma high timing threshold
    CDC_TL_DEF  =  25;      // 1 sigma low timing threshold

    CDC_IBIT_DEF =  4;      // Scaling factor for integral
    CDC_ABIT_DEF =  3;      // Scaling factor for amplitude
    CDC_PBIT_DEF =  0;      // Scaling factor for pedestal

    FDC_WS_DEF  =  16;      // hit window start - must be >= F125_FDC_NP
    FDC_WE_DEF  =  59;      // hit window end - must be at least 20 less than number of samples available, eg WE=59 for 80 samples
    FDC_IE_DEF  =  16;      // end integration at the earlier of WE, or this many samples after threshold crossing of TH  
    FDC_P1_DEF  =   4;      // 2**P1 = # samples used for pedestal used to find hit
    FDC_P2_DEF  =   4;      // 2**P2 = # samples used for pedestal calculated just before hit, and returned
    FDC_PG_DEF  =   4;      // # samples between hit threshold crossing and local pedestal sample
    FDC_H_DEF   =  50;      // 5 sigma hit threshold
    FDC_TH_DEF  =  40;      // 4 sigma high timing threshold
    FDC_TL_DEF  =  10;      // 1 sigma low timing threshold

    FDC_IBIT_DEF =  4;      // Scaling factor for integral
    FDC_ABIT_DEF =  0;      // Scaling factor for amplitude
    FDC_PBIT_DEF =  0;      // Scaling factor for pedestal

    
    // Default values for parameters which can be set in the command line
    // The parameters are ignored if their value is negative or for IBIT, below -7.
    // If set in the command line, they override both BORConfig and the default parameter values above
    // The default parameters above are used if BORConfig is not found and there are no command line parameter overrides

    CDC_WS = -9;         // hit window start - must be >= F125_CDC_NP
    CDC_WE = -9;         // hit window end - must be at least 20 less than number of samples available, eg WE=179 for 200 samples
    CDC_IE = -9;         // end integration at the earlier of WE, or this many samples after threshold crossing of TH  
    CDC_P1 = -9;         // 2**P1 = # samples used for initial pedestal, used to find hit
    CDC_P2 = -9;         // 2**P2 = # samples used for pedestal calculated just before hit, and returned 
    CDC_PG = -9;         // # samples between hit threshold crossing and local pedestal sample 'pedestal gap'
    CDC_H  = -9;         // 5 sigma hit threshold
    CDC_TH = -9;         // 4 sigma high timing threshold
    CDC_TL = -9;         // 1 sigma low timing threshold

    CDC_IBIT = -9;       // Scaling factor for integral
    CDC_ABIT = -9;       // Scaling factor for amplitude 
    CDC_PBIT = -9;       // Scaling factor for pedestal


    FDC_WS = -9;         // hit window start - must be >= F125_FDC_NP
    FDC_WE = -9;         // hit window end - must be at least 20 less than number of samples available, eg WE=179 for 200 samples
    FDC_IE = -9;         // end integration at the earlier of WE, or this many samples after threshold crossing of TH  
    FDC_P1 = -9;         // 2**P1 = # samples used for initial pedestal, used to find hit
    FDC_P2 = -9;         // 2**P2 = # samples used for pedestal calculated just before hit, and returned 
    FDC_PG = -9;         // # samples between hit threshold crossing and local pedestal sample 'pedestal gap'
    FDC_H  = -9;         // 5 sigma hit threshold
    FDC_TH = -9;         // 4 sigma high timing threshold
    FDC_TL = -9;         // 1 sigma low timing threshold

    FDC_IBIT = -9;       // Scaling factor for integral
    FDC_ABIT = -9;       // Scaling factor for amplitude 
    FDC_PBIT = -9;       // Scaling factor for pedestal


    
    // Set verbosity
    VERBOSE = 0;

    if(gPARMS){
        gPARMS->SetDefaultParameter("EMULATION125:FORCE_DEFAULT_CDC", FORCE_DEFAULT_CDC,"Set to >0 to force use of CDC default values");
        gPARMS->SetDefaultParameter("EMULATION125:FORCE_DEFAULT_FDC", FORCE_DEFAULT_FDC,"Set to >0 to force use of FDC default values");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_WS",  CDC_WS,  "Set CDC_WS for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_WE",  CDC_WE,  "Set CDC_WE for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_IE",  CDC_IE,  "Set CDC_IE for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_P1",  CDC_P1,  "Set CDC_P1 for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_P2",  CDC_P2, "Set CDC_P2 for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_PG",  CDC_PG,  "Set CDC_PG for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_H",   CDC_H,   "Set CDC_H (hit threshold) for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_TH",  CDC_TH,  "Set CDC_TH for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_TL",  CDC_TL,  "Set CDC_TL for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_IBIT",CDC_IBIT,"Set CDC_IBIT for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_ABIT",CDC_ABIT,"Set CDC_ABIT for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:CDC_PBIT",CDC_PBIT,"Set CDC_PBIT for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_WS",  FDC_WS,  "Set FDC_WS for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_WE",  FDC_WE,  "Set FDC_WE for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_IE",  FDC_IE,  "Set FDC_IE for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_P1",  FDC_P1,  "Set FDC_P1 for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_P2",  FDC_P2,  "Set FDC_P2 for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_PG",  FDC_PG,  "Set FDC_PG for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_H",   FDC_H,   "Set FDC_H (hit threshold) for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_TH",  FDC_TH,  "Set FDC_TH for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_TL",  FDC_TL,  "Set FDC_TL for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_IBIT",FDC_IBIT,"Set FDC_IBIT for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_ABIT",FDC_ABIT,"Set FDC_ABIT for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:FDC_PBIT",FDC_PBIT,"Set FDC_PBIT for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION125:VERBOSE", VERBOSE,"Set verbosity for f125 emulation");
    }
}

void Df125EmulatorAlgorithm_v2::EmulateFirmware(const Df125WindowRawData *rawData, Df125CDCPulse *cdcPulse, Df125FDCPulse *fdcPulse){

    if (VERBOSE > 0) {
        jout << "=== Entering f125 Firmware Emulation === ROCID: " << rawData->rocid << " SLOT: " << rawData->slot << " CHANNEL: " << rawData->channel <<endl;
    }

    // This is the main routine called by JEventSource_EVIO::GetObjects() and serves as the entry point for the code.
    bool isCDC = cdcPulse != NULL ? true : false;
    bool isFDC = fdcPulse != NULL ? true : false;

    if (isCDC && isFDC){
        jout << " Df125EmulatorAlgorithm_v2::EmulateFirmware Both FDC and CDC words present??? " << endl;
        return;
    } else if (!isCDC && !isFDC){ 
        jout << " Df125EmulatorAlgorithm_v2::EmulateFirmware Neither FDC or CDC words present??? " << endl;
        return;
    }

    // channel is needed for the config lookup
    uint32_t channel = rawData->channel;

    // The following are the essential values needed for the emulation 
    // (will use ROOT types since that is what the existing f125_algos code uses)
    Int_t WS=0,WE=0,IE=0,P1=0,P2=0,PG=0,H=0,TH=0,TL=0;
    Int_t IBIT=0, ABIT=0, PBIT=0;

    Int_t NE = 20; // This is a hardcoded constant in the firmware WE = NW - NE - 1  

    // Field max (saturation) values

    Int_t CDC_IMAX = 16383; //field max for integral
    Int_t CDC_AMAX = 511; //field max for max amp
    Int_t CDC_PMAX = 255; //field max for pedestal
    Int_t CDC_OMAX = 7; // field max for overflows

    Int_t FDC_IMAX = 4095; //field max for integral
    Int_t FDC_AMAX = 4095; //field max for max amp
    Int_t FDC_PMAX = 2047; //field max for pedestal
    Int_t FDC_OMAX = 7; // field max for overflows


    // Now try to get the configuration form the BOR record, if this does not exist,
    // or is forced, use the default values.
    const Df125BORConfig *BORConfig = NULL;
    //rawData->GetSingle(BORConfig);  // DEBUG

    if(BORConfig != NULL){
        if (VERBOSE > 0) jout << "--Found BORConfig" << endl; 
        Int_t NW = BORConfig->fe[0].nw;
        P1 = BORConfig->fe[0].ped_sf & FA125_FE_PED_SF_P1_MASK;
        P2 = (BORConfig->fe[0].ped_sf & FA125_FE_PED_SF_P2_MASK) >> 8;
        if (NW != Int_t(rawData->samples.size())) jout << "WARNING Df125EmulatorAlgorithm_v2::EmulateFirmware NW != rawData->samples.size()" << endl;
        WE  = NW - NE - 1;
        IE  = BORConfig->fe[0].ie & FA125_FE_IE_INTEGRATION_END_MASK;
        PG  = (BORConfig->fe[0].ie & FA125_FE_IE_PEDESTAL_GAP_MASK) >> 12;
        WS  = 1<<P1;
        H   = BORConfig->fe[channel/6].threshold[channel%6] & FA125_FE_THRESHOLD_MASK;
        TH  = (BORConfig->fe[channel/6].timing_thres_hi[(channel/3)%2] & FA125_FE_TIMING_THRES_HI_MASK(channel)) >> ((channel%3) * 9);
        TL  = (BORConfig->fe[channel/6].timing_thres_lo[(channel/2)%3] & FA125_FE_TIMING_THRES_LO_MASK(channel)) >> ((channel%2) * 16 + 8);
        IBIT = (BORConfig->fe[0].ped_sf & FA125_FE_PED_SF_IBIT_MASK)>>16;
        ABIT = (BORConfig->fe[0].ped_sf & FA125_FE_PED_SF_ABIT_MASK)>>19;
        PBIT = (BORConfig->fe[0].ped_sf & FA125_FE_PED_SF_PBIT_MASK)>>22;
        if (BORConfig->fe[0].ped_sf & FA125_FE_PED_SF_PBIT_MASK) PBIT = -1*PBIT;
    }

    // Set defaults if BORConfig missing or if forced
    if( isCDC && ((BORConfig == NULL) || FORCE_DEFAULT_CDC) ){
        if (VERBOSE > 0) jout << "WARNING Df125EmulatorAlgorithm_v2::EmulateFirmware Using CDC Default values" << endl;
        WS  = CDC_WS_DEF;
        WE  = CDC_WE_DEF;
        IE  = CDC_IE_DEF;
        P1  = CDC_P1_DEF;
        P2  = CDC_P2_DEF;
        PG  = CDC_PG_DEF;
        H   = CDC_H_DEF;
        TH  = CDC_TH_DEF;
        TL  = CDC_TL_DEF;
        IBIT= CDC_IBIT_DEF;
        ABIT= CDC_ABIT_DEF;
        PBIT= CDC_PBIT_DEF;
    }



    //Set override params if given in command line  PBIT can be negative
    if ( isCDC ) {
      if (VERBOSE > 0) jout << "WARNING Df125EmulatorAlgorithm_v2::EmulateFirmware Using CDC override values" << endl;
      if (CDC_WS >= 0) WS  = CDC_WS;
      if (CDC_WE >= 0) WE  = CDC_WE;
      if (CDC_IE >= 0) IE  = CDC_IE;
      if (CDC_P1 >= 0) P1  = CDC_P1;
      if (CDC_P2 >= 0) P2  = CDC_P2;
      if (CDC_PG >= 0) PG  = CDC_PG;
      if (CDC_H >= 0)  H   = CDC_H;
      if (CDC_TH >= 0) TH  = CDC_TH;
      if (CDC_TL >= 0) TL  = CDC_TL;
      if (CDC_IBIT >= 0) IBIT = CDC_IBIT;
      if (CDC_ABIT >= 0) ABIT = CDC_ABIT;
      if (CDC_PBIT >= -7) PBIT = CDC_PBIT;
    }


    if( isFDC && ((BORConfig == NULL) || FORCE_DEFAULT_FDC) ){
        if (VERBOSE > 0) jout << "WARNING Df125EmulatorAlgorithm_v2::EmulateFirmware Using FDC Default values" << endl;
        WS  = FDC_WS_DEF;
        WE  = FDC_WE_DEF;
        IE  = FDC_IE_DEF;
        P1  = FDC_P1_DEF;
        P2  = FDC_P2_DEF;
        PG  = FDC_PG_DEF;
        H   = FDC_H_DEF;
        TH  = FDC_TH_DEF;
        TL  = FDC_TL_DEF;
        IBIT= FDC_IBIT_DEF;
        ABIT= FDC_ABIT_DEF;
        PBIT= FDC_PBIT_DEF;
    }


    //Set override params if given in command line  PBIT can be -ve
    if ( isFDC ) {
      if (VERBOSE > 0) jout << "WARNING Df125EmulatorAlgorithm_v2::EmulateFirmware Using FDC override values" << endl;
      if (FDC_WS >= 0) WS  = FDC_WS;
      if (FDC_WE >= 0) WE  = FDC_WE;
      if (FDC_IE >= 0) IE  = FDC_IE;
      if (FDC_P1 >= 0) P1  = FDC_P1;
      if (FDC_P2 >= 0) P2  = FDC_P2;
      if (FDC_PG >= 0) PG  = FDC_PG;
      if (FDC_H >= 0)  H   = FDC_H;
      if (FDC_TH >= 0) TH  = FDC_TH;
      if (FDC_TL >= 0) TL  = FDC_TL;
      if (FDC_IBIT >= 0) IBIT = FDC_IBIT;
      if (FDC_ABIT >= 0) ABIT = FDC_ABIT;
      if (FDC_PBIT >= -7) PBIT = FDC_PBIT;
    }



    if (VERBOSE > 0) {
        jout << "============ Parameters used for emulation ================" << endl;
        jout << "WS: " << WS << " WE: " << WE << " IE: " << IE << " P1: " << P1 << " P2: " << P2 << endl;
        jout << "PG: " << PG << " H: " << H << " TH: " << TH << " TL: " << TL << endl;
        jout << "IBIT: " << IBIT << " ABIT: " << ABIT << " PBIT: " << PBIT << endl;
    }

    // The calculated quantities are passed by reference
    Int_t time=0, q_code=0, pedestal=0, overflows=0, maxamp=0, pktime=0;
    Long_t integral=0;

   

    // Perform the emulation
    fa125_algos(time, q_code, pedestal, integral, overflows, maxamp, pktime, &rawData->samples[0], WS, WE, IE, P1, P2, PG, H, TH, TL);

    // Scale down
    
    integral = integral >> IBIT;
    maxamp = maxamp >> ABIT;
    pedestal = pedestal >> (P2 + PBIT);


    // Put the emulated values in the objects
    if (isCDC){
        cdcPulse->le_time_emulated = time;
        cdcPulse->time_quality_bit_emulated = q_code;
        cdcPulse->overflow_count_emulated = ( overflows > CDC_OMAX) ? CDC_OMAX : overflows;
        cdcPulse->pedestal_emulated = ( pedestal > CDC_PMAX ) ? CDC_PMAX : pedestal;
        cdcPulse->integral_emulated = ( integral > CDC_IMAX ) ? CDC_IMAX : integral;
        cdcPulse->first_max_amp_emulated = ( maxamp > CDC_AMAX ) ? CDC_AMAX : maxamp;
    }
    else{
        fdcPulse->le_time_emulated = time;
        fdcPulse->time_quality_bit_emulated = q_code;
        fdcPulse->overflow_count_emulated =  ( overflows > FDC_OMAX) ? FDC_OMAX : overflows;
        fdcPulse->pedestal_emulated = ( pedestal > FDC_PMAX ) ? FDC_PMAX : pedestal;
        fdcPulse->integral_emulated = ( integral > FDC_IMAX ) ? FDC_IMAX : integral;
        fdcPulse->peak_amp_emulated = ( maxamp > FDC_AMAX ) ? FDC_AMAX : maxamp;
        fdcPulse->peak_time_emulated = pktime;
    }

    // Copy the emulated values to the main values if needed
    if (isCDC && cdcPulse->emulated){
        cdcPulse->le_time = cdcPulse->le_time_emulated;
        cdcPulse->time_quality_bit = cdcPulse->time_quality_bit_emulated;
        cdcPulse->overflow_count = cdcPulse->overflow_count_emulated;
        cdcPulse->pedestal = cdcPulse->pedestal_emulated;
        cdcPulse->integral = cdcPulse->integral_emulated;
        cdcPulse->first_max_amp = cdcPulse->first_max_amp_emulated;
    }
    else if (isFDC && fdcPulse->emulated){
        fdcPulse->le_time = fdcPulse->le_time_emulated;
        fdcPulse->time_quality_bit = fdcPulse->time_quality_bit_emulated;
        fdcPulse->overflow_count = fdcPulse->overflow_count_emulated;
        fdcPulse->pedestal = fdcPulse->pedestal_emulated;
        fdcPulse->integral = fdcPulse->integral_emulated;
        fdcPulse->peak_amp = fdcPulse->peak_amp_emulated;
        fdcPulse->peak_time = fdcPulse->peak_time_emulated;
    }

    if (VERBOSE > 0) jout << "=== Exiting f125 Firmware Emulation === " << endl;

    return;
}

void Df125EmulatorAlgorithm_v2::fa125_algos(Int_t &time, Int_t &q_code, Int_t &pedestal, Long_t &integral, Int_t &overflows, Int_t &maxamp, Int_t &pktime, const uint16_t adc[], Int_t WINDOW_START, Int_t WINDOW_END, Int_t INT_END, Int_t P1, Int_t P2, Int_t PG, Int_t HIT_THRES, Int_t HIGH_THRESHOLD, Int_t LOW_THRESHOLD) {

    const Int_t NU = 20;  //number of samples sent to time algo
    const Int_t PED = 5;  //sample to be used as pedestal for timing is in place 5

    const Int_t XTHR_SAMPLE = PED + PG;
    Int_t adc_subset[NU];

    Int_t hitfound=0; //hit found or not (1=found,0=not)
    Int_t hitsample=-1;  // if hit found, sample number of threshold crossing
    Int_t timesample=0;

    Int_t i=0;

    time=0;       // hit time in 0.1xsamples since start of buffer passed to fa125_time
    q_code=-1;    // quality code, 0=good, 1=returned rough estimate
    pedestal=0;   // pedestal just before hit
    integral=0;   // signal integral, total
    overflows=0;  // count of samples with overflow bit set (need raw data, not possible from my root files)
    maxamp=0;     // signal amplitude at first max after hit
    pktime=0;     // sample number containing maxamp

    // look for hit using mean pedestal of NPED samples before trigger
    fa125_hit(hitfound, hitsample, pedestal, adc, WINDOW_START, WINDOW_END, HIT_THRES, P1, P2, PG);

    if (hitfound==1) {
        for (i=0; i<NU; i++) {
            adc_subset[i] = adc[hitsample+i-XTHR_SAMPLE];
        }

        fa125_time(time, q_code, adc_subset, NU, PG, HIGH_THRESHOLD, LOW_THRESHOLD);

        timesample = hitsample-XTHR_SAMPLE + (Int_t)(0.1*time);  //sample number containing leading edge sample

        fa125_integral(integral, overflows, timesample, adc, WINDOW_END, INT_END);

        fa125_max(maxamp, pktime, timesample, adc, WINDOW_END);

        time = 10*(hitsample-XTHR_SAMPLE) + time;   // integer number * 0.1 samples
    }
}

void Df125EmulatorAlgorithm_v2::fa125_hit(Int_t &hitfound, Int_t &hitsample, Int_t &pedestal, const uint16_t adc[], Int_t WINDOW_START, Int_t WINDOW_END, Int_t HIT_THRES, Int_t P1, Int_t P2, Int_t PG) {

    pedestal=0;  //pedestal
    Int_t threshold=0;
    Int_t i=0;

    Int_t NPED = 1<<P1;
    Int_t NPED2 = 1<<P2;

    // calc pedestal as mean of NPED samples before trigger
    for (i=0; i<NPED; i++) {
        pedestal += adc[WINDOW_START-NPED+i];
    }

    pedestal = pedestal>>P1; 
    threshold = pedestal + HIT_THRES;



    // look for threshold crossing
    i = WINDOW_START - 1 + PG;
    hitfound = 0;

    while ((hitfound==0) && (i<WINDOW_END-1)) {
        i++;
        if (adc[i] >= threshold) {
            if (adc[i+1] >= threshold) {
                hitfound = 1;
                hitsample = i;
            }
        }
    }


    if (hitfound == 1) {
        //calculate INTEGRATED pedestal ending just before the hit (this is rightshifted by P2+PBIT later on)
        pedestal = 0;
        for (i=0; i<NPED2; i++) {
            pedestal += adc[hitsample-PG-i];
        }
    }


}

void Df125EmulatorAlgorithm_v2::fa125_integral(Long_t& integral, Int_t& overflows, Int_t timesample, const uint16_t adc[], Int_t WINDOW_END, Int_t INT_END) {

    Int_t i=0;

    integral = 0;
    overflows = 0;

    if (timesample <= WINDOW_END) {

        Int_t lastsample = timesample + INT_END - 1;

        if (lastsample > WINDOW_END) lastsample = WINDOW_END;

        for (i = timesample; i <= lastsample; i++ ) {

            integral += (Long_t)adc[i];
            if (adc[i]==(Long_t)4095) overflows++;   

        }

    }

}

void Df125EmulatorAlgorithm_v2::fa125_max(Int_t& maxamp, Int_t &maxsample, Int_t hitsample, const uint16_t adc[], Int_t WINDOW_END) {

    int i;
    int ndec = 0;  //number of decreasing samples

    //make sure we are on an up-slope

    while ((adc[hitsample] <= adc[hitsample-1]) && (hitsample <= WINDOW_END)) hitsample++;

    maxamp = adc[hitsample];
    maxsample = hitsample;

    for (i=hitsample; i<=WINDOW_END; i++) {

        if (adc[i] > adc[i-1]) {
            maxamp = adc[i];
            maxsample = i;
            ndec = 0;
        }

        if (adc[i] <= adc[i-1]) ndec++;
        if (ndec==2) break;
    }

    if (hitsample >= WINDOW_END) {
      maxamp = adc[WINDOW_END];
      maxsample = WINDOW_END;
    }

}

void Df125EmulatorAlgorithm_v2::fa125_time(Int_t &le_time, Int_t &q_code, Int_t adc[], Int_t NU, Int_t PG, Int_t THRES_HIGH, Int_t THRES_LOW) {

    // adc[NU]     array of samples
    // NU=20       size of array
    // PG          pedestal gap - hit threshold crossing is PG samples after adc[PED] - the single sample to be used as pedestal here
    // THRES_HIGH  high timing threshold (eg 4 x pedestal-width )
    // THRES_LOW   high timing threshold (eg 1 x pedestal-width )
    //
    // le_time     leading edge time as 0.1x number of samples since first sample supplied
    // q_code      quality code, 0=good, >0=rough estimate (firmware returns 0=good, 1=not so good)
    //
    // q_code  Time returned       Condition
    // 0       Leading edge time   Good
    // 1       X*10 - 29           Sample value of 0 found
    // 1       X*10 - 28           Sample value greater than PED_MAX found in adc[0 to PED]
    // 1       X*10 - 27           Sample values lie below the high timing threshold
    // 1       TCL*10 + 4          Low timing threshold crossing sample TCL occurs too late to upsample
    // 1       TCL*10 + 5          One or more upsampled values are negative
    // 1       TCL*10 + 9          The upsampled values are too low
    //
    const Int_t NUPSAMPLED = 6;       // number of upsampled values to calculate, minimum is 6
    const Int_t SET_ADC_MIN = 20;     // adjust adc values so that the min is at 20
    const Int_t LIMIT_PED_MAX = 511;  // max acceptable value in adc[0 to PED]
    const Int_t PED = 5;              // take local pedestal to be adc[PED]

    const Int_t START_SEARCH = PED+1; // -- start looking for hi threshold xing with this sample

    const Int_t X = PED + PG;         // hit threshold crossing sample is adc[X]
    const Int_t ROUGH_TIME = (X*10)-30; // -- add onto this to return rough time estimates

    Int_t iubuf[NUPSAMPLED] = {0};  // array of upsampled values; iubuf[0] maps to low thres xing sample

    Int_t adc_thres_hi = 0; // high threshold
    Int_t adc_thres_lo = 0; // low threshold

    //    -- contributions to hit time, these are summed together eventually, units of sample/10
    Int_t itime1 = 0; // which sample
    Int_t itime2 = 0; // which minisample
    Int_t itime3 = 0; // correction from interpolation

    //    -- search vars
    Int_t adc_sample_hi = 0; // integer range 0 to NU := 0;  --sample number for adc val at or above hi thres
    Int_t adc_sample_lo = 0; // integer range 0 to NU := 0;  -- sample num for adc val at or below lo thres
    Int_t adc_sample_lo2 = 0; // integer range 0 to 12:= 0;  -- minisample num for adc val at or below lo thres

    Bool_t over_threshold = kFALSE;
    Bool_t below_threshold = kFALSE;

    // upsampling checks
    Int_t ups_adjust = 0;
    Int_t i = 0;


    //check all samples are >0
    //check all samples from 0 to pedestal are <= LIMIT_PED_MAX
    //adc=zero and pedestal limit checks are in same fadc clock cycle

    while (i<NU) {
      if (adc[i] == 0) {
        le_time = ROUGH_TIME + 1;
        q_code = 1;
      }

      if ((i < PED+1) && (adc[i] > LIMIT_PED_MAX)) {
        le_time = ROUGH_TIME + 2;
        q_code = 2;
      }
      i++;
    }

    if (q_code>0) return; 



    //  add offset to move min val in subset equal to SET_ADC_MIN
    //  this is to move samples away from 0 to avoid upsampled pts going -ve (on a curve betw 2 samples)
    Int_t adcmin = 4095;
    i=0;
    while (i<NU) {
        if (adc[i] < adcmin) {
            adcmin = adc[i];
        }
        i++;
    }

    Int_t adcoffset = SET_ADC_MIN - adcmin;
    i=0;
    while (i<NU) {
        adc[i] = adc[i] + adcoffset;
        if (adc[i] > 4095) adc[i] = 4095;  //saturate
        i++;
    }

    // eg if adcmin is 100, setmin is 30, adcoffset = 30 - 100 = -70, move adc down by 70

    //////////////////////////////

    // calc thresholds
    adc_thres_hi = adc[PED] + THRES_HIGH;
    adc_thres_lo = adc[PED] + THRES_LOW;

    // search for high threshold crossing
    over_threshold = kFALSE;
    i = START_SEARCH;

    while ((!over_threshold)&&(i<NU)) {
        if (adc[i] >= adc_thres_hi) {
            adc_sample_hi = i;
            over_threshold = kTRUE;
        }
        i++;
    }

    if (!over_threshold) {
        le_time = ROUGH_TIME + 3;
        q_code = 3;
        return;
    }

    // search for low threshold crossing
    below_threshold = kFALSE;
    i = adc_sample_hi-1;

    while ((!below_threshold) && (i>=PED)) {
        if (adc[i] <= adc_thres_lo) {
            adc_sample_lo = i;
            itime1 = i*10;
            below_threshold = kTRUE;
        }
        i--;
    }

    if (adc[adc_sample_lo] == adc_thres_lo) {   // no need to upsample
        le_time = itime1;
        q_code = 0;
        return;
    }

    if (adc_sample_lo > NU-7) {   // too late to upsample
        le_time = itime1 + 4;
        q_code = 4;
        return;
    }

    //upsample values from adc_sample_lo to adc_sample_lo + 1
    upsamplei(adc, adc_sample_lo, iubuf, NUPSAMPLED);

    //check upsampled values are >0
    Bool_t negups = kFALSE;
    i=0;
    while ((!negups)&&(i<NUPSAMPLED)) {
        if (iubuf[i] < 0 ) {
            negups = kTRUE;
        }
        i++;
    }

    if (negups) {
        le_time = itime1 + 5;
        q_code = 5;
        return;
    }

    // correct errors
    // iubuf[0] should be equal to adc[adc_sample_lo] and iubuf[5] should equal adc[adc_sample_lo+1]
    // very steep pulse gradients cause errors in upsampling with larger errors in the later values
    // match iubuf[0] to adc[adc_sample_lo] so that the threshold crossing must be at or after iubuf[0]
    ups_adjust = iubuf[0] - adc[adc_sample_lo];

    // move threshold correspondingly instead of correcting upsampled values
    adc_thres_lo = adc_thres_lo + ups_adjust;

    // check that threshold crossing lies within the range of iubuf[0 to 5]
    if (iubuf[NUPSAMPLED-1]<= adc_thres_lo) { //bad upsampling
        le_time = itime1 + 9;   //midway
        q_code = 6;
        return;
    }

    // search through upsampled array
    below_threshold = kFALSE;
    i = NUPSAMPLED-2;

    while ((!below_threshold) && (i>=0)) {
        if (iubuf[i] <= adc_thres_lo) {
            adc_sample_lo2 = i;
            below_threshold = kTRUE;
        }
        i--;
    }

    if (!below_threshold) { //upsampled points did not go below thres
        printf("upsampled points did not go below threshold - should be impossible\n");
        le_time = 0;
        q_code = 9;
        return;
    }

    itime2 = adc_sample_lo2*2;  //  convert from sample/5 to sample/10

    //interpolate
    itime3 = 0;
    if (iubuf[adc_sample_lo2] != adc_thres_lo) {
        if (2*adc_thres_lo >= iubuf[adc_sample_lo2] + iubuf[adc_sample_lo2+1]) itime3 = 1;
    }
    le_time = itime1 + itime2 + itime3;  //   -- this is time from first sample point, in 1/10ths of samples
    q_code = 0;
}


void Df125EmulatorAlgorithm_v2::upsamplei(Int_t x[], Int_t startpos, Int_t z[], const Int_t NUPSAMPLED) {

    // x is array of samples
    // z is array of upsampled data
    // startpos is where to start upsampling in array x, only need to upsample a small region
    const Int_t nz = NUPSAMPLED;

    Int_t k,j,dk;
    const Int_t Kscale = 16384;
    const Int_t K[43] = {-4, -9, -13, -10, 5, 37, 82, 124, 139, 102, -1, -161, -336, -455, -436, -212, 241, 886, 1623, 2309, 2795, 2971, 2795, 2309, 1623, 886, 241, -212, -436, -455, -336, -161, -1, 102, 139, 124, 82, 37, 5, -10, -13, -9, -4};

    //don't need to calculate whole range of k possible
    //earliest value k=42 corresponds to sample 4.2
    //               k=43                sample 4.4
    //               k=46                sample 5.0

    // sample 4 (if possible) would be at k=41
    // sample 4.2                         k=42
    // sample 5                           k=46

    // sample x                           k=41 + (x-4)*5
    // sample x-0.2                       k=40 + (x-4)*5

    Int_t firstk = 41 + (startpos-4)*5;
    for (k=firstk; k<firstk+nz; k++) {
        dk = k - firstk;
        z[dk]=0.0;
        for (j=k%5;j<43;j+=5) {
            z[dk] += x[(k-j)/5]*K[j];
        }
        //    printf("dk %i z %i 5z %i  5z/scale %i\n",dk,z[dk],5.0*z[dk],5.0*z[dk]/Kscale);
        z[dk] = (Int_t)(5*z[dk])/Kscale;
    }
}
