//
// $Id$
//

// This file and the fa125algo.cc file were made by splitting the file
// fa125algo.inc.C sent to me by Naomi via e-mail on 3/18/2015. 
//
// Updated 10.27.15 NSJ

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
		Int_t WINDOW_START;
		Int_t WINDOW_END;
		Int_t INT_END;

		Int_t NPED;
		Int_t NPED2;
        Int_t PG;

		Int_t HIT_THRES;
		Int_t HIGH_THRESHOLD;
		Int_t LOW_THRESHOLD;
};



// FA125 emulation functions cdc_hit, cdc_time etc

// Peak amplitude and integral are returned without pedestal subtraction

  // cdc_time q_code values:

  // q_code  Time returned       Condition
  // 0       Leading edge time   Good 
  // 1       X*10 - 29           Sample value of 0 found 
  // 1       X*10 - 28           Sample value greater than PED_MAX found in adc[0 to PED]
  // 1       X*10 - 27           Sample values lie below the high timing threshold 
  // 1       TCL*10 + 4          Low timing threshold crossing sample TCL occurs too late to upsample
  // 1       TCL*10 + 5          One or more upsampled values are negative 
  // 1       TCL*10 + 9          The upsampled values are too low


void cdc_hit(Int_t&, Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);   // look for a hit
void cdc_time(Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t); // find hit time
void cdc_integral(Long_t&, Int_t&, Int_t, Int_t[], Int_t, Int_t); // find integral
void cdc_max(Int_t&, Int_t, Int_t[], Int_t); // find first max amplitude after hit
void upsamplei(Int_t[], Int_t, Int_t[], Int_t);   // upsample

void fa125_algos(Int_t&, Int_t&, Int_t&, Long_t&, Int_t&, Int_t&, Int_t[], Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t);

//	fa125_algos(d.time, d.q_code, d.pedestal, d.integral, d.overflows, d.maxamp, adc, d.WINDOW_START, d.WINDOW_END, d.INT_END, d.NPED, d.NPED2, d.PG, d.HIT_THRES, d.HIGH_THRESHOLD, d.LOW_THRESHOLD);

// Alternate form that wraps original defined above. (See comment just above
// fa125_algos_data_t class definition)
//void fa125_algos(int rocid, vector<uint16_t> samples, fa125_algos_data_t &fa125_algos_data);

void fa125_algos(int rocid, vector<uint16_t> samples, fa125_algos_data_t &fa125_algos_data, uint32_t CDC_WS, uint32_t CDC_WE, uint32_t CDC_IE, uint32_t CDC_NP, uint32_t CDC_NP2, uint32_t CDC_PG, uint32_t CDC_H, uint32_t CDC_TH, uint32_t CDC_TL, uint32_t FDC_WS, uint32_t FDC_WE, uint32_t FDC_IE, uint32_t FDC_NP, uint32_t FDC_NP2, uint32_t FDC_PG, uint32_t FDC_H, uint32_t FDC_TH, uint32_t FDC_TL);
