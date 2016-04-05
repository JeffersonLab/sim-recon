#include <DAQ/Df250EmulatorAlgorithm_v1.h>

Df250EmulatorAlgorithm_v1::Df250EmulatorAlgorithm_v1(JEventLoop *loop){
    // Enables forced use of default values 
    FORCE_DEFAULT = 0;
    // Default values for the essential parameters
    NSA_DEF = 20;
    NSB_DEF = 5;
    THR_DEF = 120;
    // Set verbosity
    VERBOSE = 0;

    if(gPARMS){
        gPARMS->SetDefaultParameter("EMULATION250:FORCE_DEFAULT", FORCE_DEFAULT,"Set to >0 to force use of default values");
        gPARMS->SetDefaultParameter("EMULATION250:NSA", NSA_DEF,"Set NSA for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:NSB", NSB_DEF,"Set NSB for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:THR", THR_DEF,"Set threshold for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:VERBOSE", VERBOSE,"Set verbosity for f250 emulation");
    }
}

void Df250EmulatorAlgorithm_v1::EmulateFirmware(const Df250WindowRawData* rawData,
                                                std::vector<JObject*> &pt_objs,
                                                std::vector<JObject*> &pp_objs,
                                                std::vector<JObject*> &pi_objs)
{
    // This is the main routine called by JEventSource_EVIO::GetObjects() and serves as the entry point for the code.
    if (VERBOSE > 0) {
        jout << " Df250EmulatorAlgorithm_v1::EmulateFirmware ==> Starting emulation <==" << endl;
        jout << "rocid : " << rawData->rocid << " slot: " << rawData->slot << " channel: " << rawData->channel << endl;
    }

    // First check that we have window raw data available
    if (rawData == NULL) {
        jout << " ERROR: Df250EmulatorAlgorithm_v1::EmulateFirmware - raw sample data is missing" << endl;
        jout << " Contact mstaib@jlab.org" << endl;
    } 

    // We need the channel number to get the threshold
    uint32_t channel = rawData->channel;

    // First grab the config objects from the raw data and get the quantities we need from them
    // The only things we need for this version of the f250 firmware are NSB, NSA, and the threshold.
    // These are all stored in the BOR config. We can grab this from the raw data since that was already associated in JEventSource_EVIO::GetObjects. 
    const Df250BORConfig *f250BORConfig = NULL;
    rawData->GetSingle(f250BORConfig);

    uint32_t NSA, NSB;
    uint16_t THR;
    //If this does not exist, or we force it, use the default values
    if (f250BORConfig == NULL || FORCE_DEFAULT){
        static int counter = 0;
        NSA = NSA_DEF;
        NSB = NSA_DEF;
        THR = THR_DEF;
        if (counter < 10){
            counter++;
            if (counter == 10) jout << " WARNING Df250EmulatorAlgorithm_v1::EmulateFirmware No Df250BORConfig == Using default values == LAST WARNING" << endl;
            else jout << " WARNING Df250EmulatorAlgorithm_v1::EmulateFirmware No Df250BORConfig == Using default values" << endl;
        }
    }
    else{
        NSA = f250BORConfig->adc_nsa & 0x7F;
        NSB = f250BORConfig->adc_nsb & 0x7F;
        THR = f250BORConfig->adc_thres[channel];
        if (VERBOSE > 0) jout << "Df250EmulatorAlgorithm_v1::EmulateFirmware NSA: " << NSA << " NSB: " << NSB << " THR: " << THR << endl; 
    }

    std::vector<const Df250PulseTime*> pulseTimes;
    std::vector<const Df250PulsePedestal*> pulsePedestals;
    std::vector<const Df250PulseIntegral*> pulseIntegrals;
    rawData->Get(pulseTimes);
    rawData->Get(pulsePedestals);
    rawData->Get(pulseIntegrals);
 
    // Now we can start to loop over the raw data
    // This requires a few passes due to some features in the way the quantities are calculated...
    // The first step is to scan the samples for TC (threshold crossing sample) and compute the
    // integrals of all pulses found.

    vector<uint16_t> samples = rawData->samples; 
    uint16_t NW = samples.size();
    int npulses = 0;
    const int max_pulses = 3;
    uint32_t TC[max_pulses] = {};
    uint32_t TMIN[max_pulses] = {3};
    uint32_t pulse_integral[max_pulses] = {};

    for (unsigned int i=0; i < NW; i++) {
        if (VERBOSE > 5) jout << "Df250EmulatorAlgorithm_v1::EmulateFirmware samples[" << i << "]: " << samples[i] << endl;
    }

    for (unsigned int i=0; i < NW; i++) {
        if ((samples[i] & 0xfff) > THR) {
            TC[npulses] = i+1;
            unsigned int ibegin = i > NSB ? (i - NSB) : 0; // Set to beginning of window if too early
            unsigned int iend = (i + NSA) < uint32_t(NW) ? (i + NSA) : NW; // Set to last sample if too late
            for (i = ibegin; i < iend; ++i)
                pulse_integral[npulses] += (samples[i] & 0xfff);
            for (; i < NW && (samples[i] & 0xfff) > THR; ++i) {}
            if (++npulses == max_pulses)
               break;
            TMIN[npulses] = i;
        }
    }

    // That concludes the first pass over the data.
    // Now we can head into the fine timing pass over the data.

    uint32_t VPEAK[max_pulses] = {};
    uint16_t TMID[max_pulses] = {};
    uint16_t VMID[max_pulses] = {};
    uint16_t TFINE[max_pulses] = {};
    uint32_t pulse_time[max_pulses] = {};
    bool reportTC[max_pulses] = {};

    // The first 4 samples are used for the pedestal calculation
    uint32_t VMIN = 0;
    for (unsigned int i=0; i < 4; i++) {
        VMIN += (samples[i] & 0xfff);
        // There is an error if any of the first four samples are above threshold
        // "A problem in the computation of the high resolution time will occur when VMIN is greater than VPEAK.  
        // In the current implementation of the algorithm, the simplest way to protect against this situation is 
        // to require that all 4 samples that determine the VMIN must be at or below threshold for the high resolution 
        // timing algorithm to be used.  If this condition is not satisfied, the reported pulse time is TC, and 
        // both VMIN and VPEAK are reported as 0 in the pulse parameter data word (type 10) to identify the condition." 
        if ((samples[i] & 0xfff) > THR){
            if (VERBOSE > 1) jout << "WARNING Df250EmulatorAlgorithm_v1::EmulateFirmware sample above threshold in pedestal calculation " << endl;
            VMIN = 99999 << 2;
            break;
        }
    }
    VMIN = VMIN >> 2;

    for (int p=0; p < npulses; ++p) {
        while (true) {
            if (VMIN == 99999) {
                VPEAK[p] = 0;
                reportTC[p] = true;
                pulse_time[p] = (TC[p] << 6);
                break;
            }
            int ipeak;
            for (ipeak = TC[p]; ipeak < NW; ++ipeak) {
                if ((samples[ipeak] & 0xfff) < (samples[ipeak-1] & 0xfff)) {
                    VPEAK[p] = (samples[ipeak-1] & 0xfff);
                    break;
                }
            }
        // There is an error condition if the TC sample is within 5 from the end of the window
        // (and a typo in the document, the firmware is looking for (NW - TC) <= 5).
        // "In the current implementation of the algorithm, a technical difficulty arises when TC is near the 
        // end of the trigger window.  If (NW - TC) < 5, the reported pulse time is TC, and the pulse parameter 
        // data word (type 10) reports VPEAK = 0 and VMIN as measured. 
        //
        // Note by RTJ:
        // I believe that the algorithmic glitch is associated with (NW - TPEAK) < 5,
        // which Mike may have found often corresponds to (NW - TC) <= 5, but not always.
            if (NW - ipeak < 5) {
                VPEAK[p] = 0;
                pulse_time[p] = ((TC[p] - 1) << 6);
                reportTC[p] = true;
                break;
            }
        // If the peak search failed, there is another special error condition
        // "A problem with the algorithm occurs if VPEAK is not found within the trigger window. In this case, the reported 
        // pulse time is TC. To identify this condition, the pulse parameter data word (type 10) reports VPEAK = 0 and VMIN as measured." 
            if (VPEAK[p] == 0) {
                pulse_time[p] = (TC[p] << 6);
                reportTC[p] = true;
                break;
            }

            if (VERBOSE > 1) {
                jout << " pulse " << p << ": VMIN: " << VMIN << " reportTC: " << int(reportTC[p]) 
                     << " TC: " << TC[p] << " VPEAK: " << VPEAK[p] << endl;  
            }

            VMID[p] = (VMIN + VPEAK[p]) >> 1;
            for (unsigned int i = TMIN[p] + 1; i < (uint32_t)ipeak; ++i) {
                if ((samples[i] & 0xfff) > VMID[p]) {
                    TMID[p] = i;
                    break;
                }
            }
            if (TMID[p] == 0) {
#if EMULATION250_MODE_8
                if (p == 0) {
#else
                if (false) {
#endif
                    TMID[p] = 1;
                    TFINE[p] = 0; // empirical constant
                }
                else {
                    TMID[p] = TC[p];
                    TFINE[p] = 0;
                }
            }
            else {
                int Vnext = (samples[TMID[p]] & 0xfff);
                int Vlast = (samples[TMID[p]-1] & 0xfff);
                if (Vnext > Vlast && VMID[p] >= Vlast)
                    TFINE[p] = 64 * (VMID[p] - Vlast) / (Vnext - Vlast);
                else
                    TFINE[p] = 62;
            }
            pulse_time[p] = (TMID[p] << 6) + TFINE[p];
            break;
        }
        VMIN = (VMIN < 99999)? VMIN : 0;

        if (VERBOSE > 1) {
            jout << " pulse " << p << ": VMID: " << VMID[p] << " TMID: " << TMID[p] 
                 << " TFINE: " << TFINE[p] << " time: " << pulse_time[p]
                 << " integral: " << pulse_integral[p] << endl;
        }

        Df250PulseTime* f250PulseTime = new Df250PulseTime;
        f250PulseTime->rocid = rawData->rocid;
        f250PulseTime->slot =  rawData->slot;
        f250PulseTime->channel = rawData->channel;
        f250PulseTime->itrigger = rawData->itrigger;
        f250PulseTime->pulse_number = p;
        f250PulseTime->quality_factor = reportTC[p];
        f250PulseTime->time = pulse_time[p];
        f250PulseTime->emulated = true;
        f250PulseTime->time_emulated = pulse_time[p];
        f250PulseTime->quality_factor_emulated = reportTC[p];
        f250PulseTime->AddAssociatedObject(rawData);
        pt_objs.push_back(f250PulseTime);

        Df250PulsePedestal* f250PulsePedestal = new Df250PulsePedestal;
        f250PulsePedestal->rocid = rawData->rocid;
        f250PulsePedestal->slot = rawData->slot;
        f250PulsePedestal->channel = rawData->channel;
        f250PulsePedestal->itrigger = rawData->itrigger;
        f250PulsePedestal->pulse_number = p;
        f250PulsePedestal->pedestal = VMIN;
        f250PulsePedestal->pulse_peak = VPEAK[p];
        f250PulsePedestal->emulated = true;
        f250PulsePedestal->pedestal_emulated = VMIN;
        f250PulsePedestal->pulse_peak_emulated = VPEAK[p];
        f250PulsePedestal->AddAssociatedObject(rawData);
        pp_objs.push_back(f250PulsePedestal);

        Df250PulseIntegral* f250PulseIntegral = new Df250PulseIntegral;
        f250PulseIntegral->rocid = rawData->rocid;
        f250PulseIntegral->slot = rawData->slot;
        f250PulseIntegral->channel = rawData->channel;
        f250PulseIntegral->itrigger = rawData->itrigger;
        f250PulseIntegral->pulse_number = p;
        f250PulseIntegral->quality_factor = reportTC[p];
        f250PulseIntegral->integral = pulse_integral[p];
        f250PulseIntegral->pedestal = VMIN;
        f250PulseIntegral->nsamples_integral = NSA + NSB;
        f250PulseIntegral->nsamples_pedestal = 4;
        f250PulseIntegral->emulated = true;
        f250PulseIntegral->integral_emulated = pulse_integral[p];
        f250PulseIntegral->pedestal_emulated = VMIN;
        f250PulseIntegral->AddAssociatedObject(rawData);
        pi_objs.push_back(f250PulseIntegral);
    }

    if (VERBOSE > 0) jout << " Df250EmulatorAlgorithm_v1::EmulateFirmware ==> Emulation complete <==" << endl;    
    return;
}
