#include <DAQ/Df250EmulatorAlgorithm_v1.h>

Df250EmulatorAlgorithm_v1::Df250EmulatorAlgorithm_v1(JEventLoop *loop){
    // Enables forced use of default values 
    FORCE_DEFAULT = 0;
    // Default values for the essential parameters
    NSA_DEF = 20;
    NSB_DEF = 5;
    THR_DEF = 120;
    // Set verbosity
    VERBOSE = 10;

    if(gPARMS){
        //gPARMS->
    }
}

void Df250EmulatorAlgorithm_v1::EmulateFirmware(const Df250WindowRawData* rawData, Df250PulseTime* pulseTime, Df250PulsePedestal* pulsePedestal, Df250PulseIntegral* pulseIntegral){

    if (VERBOSE > 0) {
        jout << " Df250EmulatorAlgorithm_v1::EmulateFirmware ==> Starting emulation <==" << endl;
        jout << "rocid : " << rawData->rocid << " slot: " << rawData->slot << " channel: " << rawData->channel << endl;
    }

    //First check that we have all of the words available
    if(rawData == NULL || pulseTime == NULL || pulsePedestal == NULL || pulseIntegral == NULL){
        jout << " ERROR: Df250EmulatorAlgorithm_v1::EmulateFirmware Some of the data objects are NULL" << endl;
        jout << " WRD: " << rawData << " PT: " << pulseTime << " PP: " << pulsePedestal << " PI: " << pulseIntegral << endl;
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

    // Now we can start to loop over the raw data
    // This requires a few passes due to some features in the way the quantities are calculated...
    // The first step is to find the VMIN (pedestal), TC (threshold crossing sample) and VPEAK (pulse height).
    vector<uint16_t> samples = rawData->samples; 
    uint16_t NW = samples.size();
    uint32_t VMIN = 0, VPEAK = 0, TC = 0;
    bool reportTC = false, peakSearch = false, foundPeak = false;
    for(unsigned int i=0; i< NW; i++){
        if (VERBOSE > 5) jout << "Df250EmulatorAlgorithm_v1::EmulateFirmware samples[" << i << "]: " << samples[i] << endl;
        // The first 4 samples are used for the pedestal calculation
        if (i < 4) {
            VMIN += samples[i];
            // There is an error if any of the first four samples are above threshold
            // "A problem in the computation of the high resolution time will occur when VMIN is greater than VPEAK.  
            // In the current implementation of the algorithm, the simplest way to protect against this situation is 
            // to require that all 4 samples that determine the VMIN must be at or below threshold for the high resolution 
            // timing algorithm to be used.  If this condition is not satisfied, the reported pulse time is TC, and 
            // both VMIN and VPEAK are reported as 0 in the pulse parameter data word (type 10) to identify the condition." 
            if (samples[i] > THR){
                if (VERBOSE > 1) jout << "WARNING Df250EmulatorAlgorithm_v1::EmulateFirmware sample above threshold in pedestal calculation " << endl;
                TC = i;
                VMIN = 0; VPEAK = 0; reportTC = true;
                break;
            }
        }
        else if (i == 4) VMIN = VMIN >> 2; // We made it through the pedestal calculation without anything above threshold
        if(samples[i] > THR && peakSearch == false) {
            TC = i; // Threshold crossing sample determined

            // There is an error condition if the TC sample is within 5 from the end of the window (and a typo in the document)
            // The firmware is looking for (NW - TC) <= 5
            // "In the current implementation of the algorithm, a technical difficulty arises when TC is near the 
            // end of the trigger window.  If (NW - TC) < 5, the reported pulse time is TC, and the pulse parameter 
            // data word (type 10) reports VPEAK = 0 and VMIN as measured. 
            if((NW - TC) <= 5){
                VPEAK = 0; reportTC = true;
                break;
            }
            VPEAK = samples[i];
            peakSearch = true; // begin the peak search
        } 
        else if (peakSearch == true){
            if (samples[i] >= samples[i-1]) VPEAK = samples[i];
            else {
                foundPeak = true;
                break;
            }
        }
    }

    // That concludes the first pass over the data, if the peak search failed, there is another special error condition
    // "A problem with the algorithm occurs if VPEAK is not found within the trigger window. In this case, the reported 
    // pulse time is TC.  To identify this condition, the pulse parameter data word (type 10) reports VPEAK = 0 and VMIN as measured." 

    if (!foundPeak){
        reportTC = true; VPEAK = 0;
    }

    // Now we can head into the fine timing and integral pass over the data.
    uint32_t VMID=0;
    if (foundPeak){
        VMID = (VMIN + VPEAK) >> 1;
    }

    // Need to determine the starting and finishing sample for the integral
    int16_t sample_integral_start, sample_integral_end;
    sample_integral_start = (TC - NSB) > 0 ? (TC - NSB) : 0; // Set to beginning of window if too early
    sample_integral_end = (int16_t(TC + NSA) - 1) < (int16_t(NW) - 1) ? (TC + NSA - 1) : (NW - 1); // Set to last sample if too late
    
    // Now run through the data a second time to calculate the integral and determine the fine time.
    bool foundTime = false;
    uint16_t TF = 0; // fine time
    uint32_t integral = 0;
    for(unsigned int i=0; i< NW; i++){
        if (int16_t(i) >= sample_integral_start && int16_t(i) <= sample_integral_end){
            integral += samples[i];
        }
        if (!reportTC && foundPeak && !foundTime){
            //if ( i > 3){
            if ( i > 3){ // only search in the region after the pedestal samples
                if(samples[i] > VMID){ 
                    // The line below is a bug in the current firmware...
                    // Should look for a 4ns difference between mode 7 and mode 8...
                    //TC = i-1;
                    TC = i;
                    if((samples[i] - samples[i-1]) != 0) TF = 64 * (VMID - samples[i-1]) / (samples[i] - samples[i-1]);
                    foundTime = true;
                }
            }
        }
    }


    // Now we have all of the quantities that we wanted. Set them as the emulated values in the data words
    pulseIntegral->integral_emulated = integral;
    pulseIntegral->pedestal_emulated = VMIN;

    uint32_t time = TC << 6 | TF;
    pulseTime->time_emulated = time;
    pulseTime->quality_factor_emulated = reportTC;

    pulsePedestal->pedestal_emulated = VMIN;
    pulsePedestal->pulse_peak_emulated = VPEAK;

    // if emulated, copy from the emulated quantities to the real values.

    if (VERBOSE > 0) jout << " Df250EmulatorAlgorithm_v1::EmulateFirmware ==> Emulation complete <==" << endl;    
    return;
}
