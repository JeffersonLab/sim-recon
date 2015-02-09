
//
// $Id:$
//

//
// This file was made by combining two files: cdcalgo16.inc.C and checkalgo.C
// sent by Naomi J. on 1/30/2015 to a select number of people. Below is the 
// message she sent with the e-mail:
//
//==============================================================================
//
//Hi,
//
//I have attached the C++ code which goes with the description here
//https://halldweb1.jlab.org/wiki/index.php/CDC_algo
//
//The time modules are in the inc file.
//
//The .C file is a wrapper for the modules, if you make a root file using
//the DAQTree CDC plugin (or the DAQTree plugin) then it sends the raw data
//in Df125WindowRawData to the modules for analysis and writes a root file
//with the results.  It also includes some extra branches for diagnostics
//that you might, or, more likely, might not, be interested in.
//eg root -l run1280/daqtree.root "../C/checkalgo.C+(0,0)"
//(0,0) means start w event 0 and process them all.  (0,4) means process
//events 0-4.
//
//
//The .C file also contains all the #define constants which the modules in
//the inc file need.
//I've put in suggested values for runs 1602 (fairly horrible) and 2016
//(much better but still some v noisy channels).
//
//
//There are lots of comments, and the variety of quality codes is listed in
//the inc file.
//
//This code matches the VHDL code which I sent to Cody a few days ago.  I
//hope he will be able to implement that code as is - if not, I can trim
//back the VHDL and will have to trim the C++ accordingly.
//
//
//David, if you could move this into the fa125 emulations, that would be
//great.  Thank you!
//
//
//Naomi.
//==============================================================================


#include <stdint.h>
#include <vector>
#include <iostream>
using namespace std;

#include <Rtypes.h>  // ROOT type definitions (e.g. Int_t)

// The following class is used to help wrap the cdc_algos2 call, making it easier
// to access from the JEventSource_EVIO class and to make that code more compact.
#define NUPSAMPLED 8       // number of upsampled values to calculate, minimum is 8
class cdc_algos_data_t{
	public:
		Int_t time;
		Int_t q_code;
		Int_t pedestal;
		Long_t integral;
		Int_t overflows;
		Int_t maxamp;
		Int_t initped;
		// Int_t adc[]; // Commented out since this is input
		Int_t iubuf[NUPSAMPLED+1];
		Int_t ups_err1;
		Int_t ups_err2;
		Int_t subsetmax;
		Int_t subsetmin;
		Int_t subsetstart;
};

// functions cdc_hit, cdc_time etc


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



void cdc_hit(Int_t&, Int_t&, Int_t&, Int_t[], Int_t[]);   // look for a hit
void cdc_time(Int_t&, Int_t&, Int_t&, Int_t[], Int_t[], Int_t&, Int_t&);           // find hit time
void cdc_integral(Long_t&, Int_t&, Int_t, Int_t[]); // find integral
void cdc_max(Int_t&, Int_t, Int_t[]); // find first max amplitude after hit
void upsamplei(Int_t[], Int_t, Int_t[]);   // upsample

void cdc_algos2(Int_t&, Int_t&, Int_t&, Long_t&, Int_t&, Int_t&, Int_t&, Int_t[], Int_t[], Int_t&, Int_t&, Int_t&, Int_t&, Int_t&);

// Alternate form that wraps original defined above. (See comment just above
// cdc_algos_data_t class definition)
void cdc_algos2(vector<uint16_t> samples, cdc_algos_data_t &cdc_algos_data);

