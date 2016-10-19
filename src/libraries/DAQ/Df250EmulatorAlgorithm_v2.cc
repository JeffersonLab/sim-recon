#include <DAQ/Df250EmulatorAlgorithm_v2.h>

Df250EmulatorAlgorithm_v2::Df250EmulatorAlgorithm_v2(JEventLoop *loop){
    // Enables forced use of default values 
    FORCE_DEFAULT = 0;
    // Default values for the essential parameters
    NSA_DEF = 20;
    NSB_DEF = 5;
    THR_DEF = 120;
    NPED_DEF = 4;
    MAXPED_DEF = 512;

    // Set verbosity
    VERBOSE = 0;

    if(gPARMS){
        gPARMS->SetDefaultParameter("EMULATION250:FORCE_DEFAULT", FORCE_DEFAULT,"Set to >0 to force use of default values");
        gPARMS->SetDefaultParameter("EMULATION250:NSA", NSA_DEF,"Set NSA for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:NSB", NSB_DEF,"Set NSB for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:THR", THR_DEF,"Set threshold for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:NPED", NPED_DEF,"Set NPED for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:MAXPED", MAXPED_DEF,"Set MAXPED for firmware emulation, will be overwritten by BORConfig if present");
        gPARMS->SetDefaultParameter("EMULATION250:VERBOSE", VERBOSE,"Set verbosity for f250 emulation");
    }
}

void Df250EmulatorAlgorithm_v2::EmulateFirmware(const Df250WindowRawData* rawData,
                                                std::vector<Df250PulseData*> &pdat_objs)
{
    // This is the main routine called by JEventSource_EVIO::GetObjects() and serves as the entry point for the code.
    if (VERBOSE > 0) {
        jout << " Df250EmulatorAlgorithm_v2::EmulateFirmware ==> Starting emulation <==" << endl;
        jout << "rocid : " << rawData->rocid << " slot: " << rawData->slot << " channel: " << rawData->channel << endl;
    }

    // First check that we have window raw data available
    if (rawData == NULL) {
        jout << " ERROR: Df250EmulatorAlgorithm_v2::EmulateFirmware - raw sample data is missing" << endl;
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
    uint32_t NPED, MAXPED;
    uint16_t THR;
    //If this does not exist, or we force it, use the default values
    if (f250BORConfig == NULL || FORCE_DEFAULT){
        static int counter = 0;
        NSA = NSA_DEF;
        NSB = NSA_DEF;
        THR = THR_DEF;
        NPED = NPED_DEF;
        MAXPED = MAXPED_DEF;
        if (counter < 10){
            counter++;
            if (counter == 10) jout << " WARNING Df250EmulatorAlgorithm_v2::EmulateFirmware No Df250BORConfig == Using default values == LAST WARNING" << endl;
            else jout << " WARNING Df250EmulatorAlgorithm_v2::EmulateFirmware No Df250BORConfig == Using default values" << endl;
        }
    }
    else{
        NSA = f250BORConfig->adc_nsa & 0x7F;
        NSB = f250BORConfig->adc_nsb & 0x7F;
        THR = f250BORConfig->adc_thres[channel];
        // set more parameters once the BOR is updated
        NPED = f250BORConfig->nped;
        MAXPED = MAXPED_DEF;
        if (VERBOSE > 0) jout << "Df250EmulatorAlgorithm_v2::EmulateFirmware NSA: " << NSA << " NSB: " << NSB << " THR: " << THR << endl; 
    }

    // quality bits
    bool bad_pedestal = false;
    bool bad_timing_pedestal = false;
    bool no_timing_calculation = false;

    // Now we can start to loop over the raw data
    // This requires a few passes due to some features in the way the quantities are calculated...
    // The first step is to scan the samples for TC (threshold crossing sample) and compute the
    // integrals of all pulses found.

    vector<uint16_t> samples = rawData->samples; 
    uint16_t NW = samples.size();
    uint32_t npulses = 0;
    const int max_pulses = 3;
    uint32_t TC[max_pulses] = {};
    uint32_t TMIN[max_pulses] = {3};
    uint32_t pulse_integral[max_pulses] = {};
    bool has_overflow_samples[max_pulses] = {false};
    bool has_underflow_samples[max_pulses] = {false};
    uint32_t number_samples_above_threshold[max_pulses] = {0};
    bool NSA_beyond_PTW[max_pulses] = {false};
    bool vpeak_beyond_NSA[max_pulses] = {false};
    bool vpeak_not_found[max_pulses] = {false};

    for (unsigned int i=0; i < NW; i++) {
        if (VERBOSE > 5) jout << "Df250EmulatorAlgorithm_v2::EmulateFirmware samples[" << i << "]: " << samples[i] << endl;
    }

    for (unsigned int i=0; i < NW; i++) {
        if ((samples[i] & 0xfff) > THR) {
            if (VERBOSE > 1) {
                jout << "threshold crossing at " << i << endl;
            }
            TC[npulses] = i+1;
            unsigned int ibegin = i > NSB ? (i - NSB) : 0; // Set to beginning of window if too early
            unsigned int iend = (i + NSA) < uint32_t(NW) ? (i + NSA) : NW; // Set to last sample if too late
            // check to see if NSA extends beyond the end of the window
            NSA_beyond_PTW[npulses] = (i + NSA) >= uint32_t(NW);
            for (i = ibegin; i < iend; ++i) {
                pulse_integral[npulses] += (samples[i] & 0xfff);
                // quality monitoring
                if(samples[i] == 0x1fff)
                    has_overflow_samples[npulses] = true;
                if(samples[i] == 0x1000)
                    has_underflow_samples[npulses] = true;
                // count number of samples within NSA that are above thresholds
                if( (i+1>=TC[npulses]) && ((samples[i] & 0xfff) > THR) )
                    number_samples_above_threshold[npulses]++;
            }
            for (; i < NW && (samples[i] & 0xfff) > THR; ++i) {}
            if (++npulses == max_pulses)
               break;
            TMIN[npulses] = i;
        }
    }

    // That concludes the first pass over the data.
    // Now we can head into the fine timing pass over the data.

    uint32_t VPEAK[max_pulses] = {};
    uint32_t TPEAK[max_pulses] = {};
    uint16_t TMID[max_pulses] = {};
    uint16_t VMID[max_pulses] = {};
    uint16_t TFINE[max_pulses] = {};
    uint32_t pulse_time[max_pulses] = {};

    // The pulse pedestal is the sum of NPED (4-15) samples at the beginning of the window
    uint32_t pedestal = 0;
    uint32_t VMIN = 0;   // VMIN is just the average of the first 4 samples, needed for timing algorithm
    for (unsigned int i=0; i < NPED; i++) {
        pedestal += (samples[i] & 0xfff);
        if(i<4)
            VMIN += (samples[i] & 0xfff);
        // error condition
        if ((samples[i] & 0xfff) > MAXPED) {
            bad_pedestal = true;
        }
    }
    VMIN /= 4;

    // error conditions for timing algorithm
    for (unsigned int i=0; i < 5; i++) {
        // "If any of the first 5 samples is greater than MaxPed but less than TET, the TDC algorithm will proceed
        // and Time Quality bit 0 will be set to 1"
        if ( ((samples[i] & 0xfff) > MAXPED) && ((samples[i] & 0xfff) < THR) ) {
            bad_timing_pedestal = true;
        }
        // "If any of the first 5 samples is greater than TET or underflow the TDC will NOT proceed..."
        // Waiit for iiit...
        if( ((samples[i] & 0xfff) > THR) || (samples[i] == 0x1000) )
            no_timing_calculation = true;
    }

    for (unsigned int p=0; p < npulses; ++p) {
        // "If any of the first 5 samples is greater than TET or underflow the TDC will NOT proceed
        //   1. pulse time is set to TC
        //   2. pulse peak is set to zero
        //   3. Time quality bits 0 and 1 are set to 1"
        if(no_timing_calculation) {
            // this will need to be changed when the timing algorithm is fixed
            VMID[p] = TC[p];
            TFINE[p] = 0;
            VPEAK[p] = 0;
            bad_timing_pedestal = true;
            vpeak_not_found[npulses] = true;
            //continue;
        }

        while (!no_timing_calculation && true) {
            //if (VMIN == 99999) {
            //    VPEAK[p] = 0;
            //    reportTC[p] = true;
            //    pulse_time[p] = (TC[p] << 6);
            //    break;
            // }
            unsigned int ipeak;
            for (ipeak = TC[p]; ipeak < NW; ++ipeak) {
                if ((samples[ipeak] & 0xfff) < (samples[ipeak-1] & 0xfff)) {
                    VPEAK[p] = (samples[ipeak-1] & 0xfff);
                    TPEAK[p] = ipeak-1;
                    break;
                }
            }
            // check to see if the peak is beyond the NSA
            if(ipeak > TC[p]+NSA)
                vpeak_beyond_NSA[p] = true;

            if (VERBOSE > 1) {
                jout << " pulse " << p << ": VMIN: " << VMIN 
                     << " TC: " << TC[p] << " VPEAK: " << VPEAK[p] << endl;  
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
            // NOTE: this should only be in the previous version
            // though something like this should set vpeak_beyond_NSA
            //if (NW - ipeak < 5) {
            //    VPEAK[p] = 0;
            //    pulse_time[p] = ((TC[p] - 1) << 6);
            //    reportTC[p] = true;
            //    break;
            //}
            // If the peak search failed, there is another special error condition
            // "A problem with the algorithm occurs if VPEAK is not found within the trigger window. In this case,
            //  the reported parameters are as follows:
            //   1. Pulse time is set to TC
            //   2. Pulse peak is set to 0
            //   3. Time quality bit 1 is set to 1
            if (VPEAK[p] == 0) { 
                TMID[p] = TC[p];
                TFINE[p] = 0;
                VPEAK[p] = 0;
                bad_timing_pedestal = true;
                vpeak_not_found[npulses] = true;
                break;
            }

            VMID[p] = (VMIN + VPEAK[p]) >> 1;
            
            for (unsigned int i = TMIN[p] + 1; i < (uint32_t)ipeak; ++i) {
                if ((samples[i] & 0xfff) > VMID[p]) {
                    TMID[p] = i;
                    break;
                }
            }
            /*
            for (unsigned int i = TPEAK[p]; i > 0; --i) {
                if ((samples[i] & 0xfff) < VMID[p]) {
                    TMID[p] = i+1;
                    break;
                }
            }
            */
            if (TMID[p] == 0) {  // redundant?
                TFINE[p] = 0;
            }
            else {
                int Vnext = (samples[TMID[p]] & 0xfff); 
                int Vlast = (samples[TMID[p]-1] & 0xfff);
                //int Vnext = (samples[TMID[p]-1] & 0xfff);  // should TMID be 1 smaller?
                //int Vlast = (samples[TMID[p]-2] & 0xfff);
                if (VERBOSE > 2) {
                    jout << "   TMIN = " << TMIN[p] << "  TMID  = " << TMID[p] << "  TPEAK = " << TPEAK[p] << endl
                         << "   VMID = " << VMID[p] << "  Vnext = " << Vnext   << "  Vlast = " << Vlast    << endl;
                }
                if (Vnext > Vlast && VMID[p] >= Vlast)
                    TFINE[p] = 64 * (VMID[p] - Vlast) / (Vnext - Vlast);
                else
                    TFINE[p] = 62;
            }
            pulse_time[p] = ((TMID[p]-1) << 6) + TFINE[p];
            break;
        }
        VMIN = (VMIN < 99999)? VMIN : 0;

        if (VERBOSE > 1) {
            jout << " pulse " << p << ": VMID: " << VMID[p] << " TMID: " << TMID[p] 
                 << " TFINE: " << TFINE[p] << " time: " << pulse_time[p]
                 << " integral: " << pulse_integral[p] << endl;
        }

        /*
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
        */

        Df250PulseData* f250PulseData;
        if( p < pdat_objs.size() ) {
            f250PulseData = pdat_objs[p];
        } else {
            // make a fresh object if one does not exist
            f250PulseData = new Df250PulseData;

            f250PulseData->rocid = rawData->rocid;
            f250PulseData->slot = rawData->slot;
            f250PulseData->channel = rawData->channel;
            f250PulseData->itrigger = rawData->itrigger;
            // word 1
            f250PulseData->event_within_block = 1;
            f250PulseData->QF_pedestal = bad_pedestal;  // is this right?
            f250PulseData->pedestal = pedestal;
            // word 2
            f250PulseData->integral = pulse_integral[p];
            f250PulseData->QF_NSA_beyond_PTW = NSA_beyond_PTW[npulses];  // is this right?
            f250PulseData->QF_overflow = has_overflow_samples[npulses];  // is this right?
            f250PulseData->QF_underflow = has_underflow_samples[npulses];  // is this right?
            f250PulseData->nsamples_over_threshold = number_samples_above_threshold[npulses];  // is this right?
            // word 3
            //f250PulseTime->time = pulse_time[p];
            f250PulseData->course_time = TMID[p]; // ?????
            f250PulseData->fine_time = TFINE[p];  // ???????
            //f250PulseData->time = pulse_time[p];
            f250PulseData->QF_vpeak_beyond_NSA = vpeak_beyond_NSA[npulses];  // is this right?
            f250PulseData->QF_vpeak_not_found = vpeak_not_found[npulses];  // is this right?
            f250PulseData->QF_bad_pedestal = bad_timing_pedestal;  // is this right?
            // other information
            f250PulseData->pulse_number = p;
            //f250PulseData->quality_factor = reportTC[p];
            f250PulseData->nsamples_integral = NSA + NSB;
            f250PulseData->nsamples_pedestal = NPED;
            f250PulseData->emulated = true;

            f250PulseData->AddAssociatedObject(rawData);
            pdat_objs.push_back(f250PulseData);
        }

        // copy over emulated values
        f250PulseData->integral_emulated = pulse_integral[p];
        f250PulseData->pedestal_emulated = pedestal;
        //f250PulseData->time_emulated = pulse_time[p]; 
        f250PulseData->pulse_peak_emulated = VPEAK[p]; 
        //f250PulseData->time_emulated = pulse_time[p]; 
        f250PulseData->course_time_emulated = TMID[p];
        f250PulseData->fine_time_emulated = TFINE[p];

        /*
        // if we are using the emulated values, copy them
        if( f250PulseData->emulated ) {
            f250PulseData->integral    = f250PulseData->integral_emulated;
            f250PulseData->pedestal    = f250PulseData->pedestal_emulated;
            //f250PulseData->time       = f250PulseData->time_emulated;
            f250PulseData->pulse_peak  = f250PulseData->pulse_peak_emulated;
            f250PulseData->course_time = f250PulseData->course_time_emulated;
            f250PulseData->fine_time   = f250PulseData->fine_time_emulated;
        }
        */
    }

    if (VERBOSE > 0) jout << " Df250EmulatorAlgorithm_v2::EmulateFirmware ==> Emulation complete <==" << endl;    
    return;
}
