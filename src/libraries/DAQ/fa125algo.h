//
// $Id$
//

// This file and the fa125algo.cc file were made by splitting the file
// fa125algo.inc.C sent to me by Naomi via e-mail on 3/18/2015. Here is
// here message from that e-mail:
//
//==============================================================================
//	Hi David,
//
//	I have attached a new version of the fa125 new algorithm emulation code
//	which has one set of constants for the CDC, another lovely set of
//	constants for the FDC, and it chooses which set to use after looking at
//	the rocid   :-)  but before calling the algorithm function (as sending
//	rocid to the algorithm code is taking the ugliness too far).
//	I put in my best guess at the constants.
//
//	Also I cleaned out some of the diagnostic-type-parameters which nobody
//	else would be interested in.
//
//	The constants which would be altered the most often would be the hit
//	window start and end, and the hit and timing thresholds.
//	eg for the CDC
//
//	#define CDC_WS 46       //hit window start
//	#define CDC_WE 150      //hit window end
//
//	#define CDC_H 125       // 5 sigma hit threshold
//	#define CDC_TH 100      // 4 sigma threshold
//	#define CDC_TL 25       // 1 sigma threshold
//
//	If you wish to make any accessible as hd_root command line options, those
//	(and the equivalent for the FDC) would be the best candidates.
//
//	The attached code can run on the DAQTree plugin output & makes a root file
//	(with a horrible name) containing its calculations.  The DAQTree itrigger
//	values are still crazy.
//
//	Please let me know if there is anything I can do to this to make it easier
//	to implement.
//
//	Best,
//
//	Naomi.
//==============================================================================


#include <stdint.h>
#include <vector>
#include <iostream>
using namespace std;

#include <Rtypes.h>  // ROOT type definitions (e.g. Int_t)

// The following class is used to help wrap the fa125_algos call, making it easier
// to access from the JEventSource_EVIO class and to make that code more compact.
//#define NUPSAMPLED 8       // number of upsampled values to calculate, minimum is 8
class fa125_algos_data_t{
	public:

		Int_t time;
		Int_t q_code;
		Int_t pedestal;
		Long_t integral;
		Int_t overflows;
		Int_t maxamp;
		// Int_t adc[]; // Commented out since this is input
		Int_t NSAMPLES;
		Int_t WINDOW_START;
		Int_t WINDOW_END;
		Int_t HIT_THRES;
		Int_t NPED;
		Int_t NPED2;
		Int_t HIGH_THRESHOLD;
		Int_t LOW_THRESHOLD;
		Int_t ROUGH_DT;
		Int_t INT_SAMPLE;
		Int_t INT_END;
		Int_t LIMIT_PED_MAX;
		Int_t LIMIT_ADC_MAX;
		Int_t XTHR_SAMPLE;
		Int_t PED_SAMPLE;
		Int_t SET_ADC_MIN;
		Int_t LIMIT_UPS_ERR;
};



// FA125 emulation functions cdc_hit, cdc_time etc


  // cdc_time q_code values:
  //
  //   0: Good
  //   1: ADC data did not go over threshold adc_thres_hi 
  //   2: Leading edge time is outside the upsampled region (cross adc_thres_lo too late in the buffer subset ) 
  //   3: Last upsampled point is <= low timing threshold
  //   4: Upsampled points did not go below low timing threshold
  //   5: ADC sample value of 0 found
  //   6: ADC sample value > ADC_LIMIT found
  //   7: Pedestal ADC[PED_SAMPLE] value > LIMIT_PED_MAX found
  //   8: Upsampled point is below 0
  //   9: Upsampling error is > SET_UPS_TOL


  // calc pedestal as mean of NPED samples before trigger and again as mean of NPED2 samples before hit
  // hit search is from samples WINDOW_START to WINDOW_END
  // pedestal is not subtracted from max amplitude



void cdc_hit(Int_t&, Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);   // look for a hit
void cdc_time(Int_t&, Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);  // find hit time
void cdc_integral(Long_t&, Int_t&, Int_t, Int_t[], Int_t, Int_t); // find integral
void cdc_max(Int_t&, Int_t, Int_t[], Int_t); // find first max amplitude after hit
void upsamplei(Int_t[], Int_t, Int_t[], const Int_t);   // upsample

void fa125_algos(Int_t&, Int_t&, Int_t&, Long_t&, Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);


// Alternate form that wraps original defined above. (See comment just above
// fa125_algos_data_t class definition)
void fa125_algos(int rocid, vector<uint16_t> samples, fa125_algos_data_t &fa125_algos_data);

