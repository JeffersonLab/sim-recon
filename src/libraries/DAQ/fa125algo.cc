//
// $Id$
//
//
// (See comments in fa125algo16.h)
//
// The #defines below might normally be placed in the header file, but
// these have some generic names and are intended only for this CDC
// timing algorithm. They may also be replaced with JANA configuration
// parameters and/or CCDB constants.
//
// The main entry point here is the fa125_algos routine. The samples
// are passed in along with references to many different return values.
// 

#include <stdlib.h>

#include <fa125algo.h>

#ifndef _DBG_
#define _DBG_ cout<<__FILE__<<":"<<__LINE__<<" "
#define _DBG__ cout<<__FILE__<<":"<<__LINE__<<endl
#endif

// FA125 emulation functions cdc_hit, cdc_time etc for NPK=1 (only the first peak is identified)

// The ADC data buffer has NW samples in total, comprising
// [unused samples][NPED samples for pedestal][samples WINDOW_START to WINDOW_END][20 or more extra samples] 
//
// WINDOW_START is the sample immediately after the NPED samples in the pedestal window.
// The hit search starts with sample WINDOW_START+PG and ends including sample WINDOW_END. This is for ease of use with older data.
// In the new firmware, there are no unused samples at the start, and 20 samples at the end.
// The first sample is numbered sample 0.
// The pedestal returned is the mean of the NPED2 samples ending PED samples before the sample containing the hit threshold crossing
//
// This requires that: 
//                     WINDOW_START >= NPED 
//                     WINDOW_END+20 < NW (number of samples supplied) 
//                     NPED2 >= NPED
//
//
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




// Wrapper for routine below
//void fa125_algos(int rocid, vector<uint16_t> samples, fa125_algos_data_t &fa125_algos_data)

void fa125_algos(int rocid, vector<uint16_t> samples, fa125_algos_data_t &fa125_algos_data, uint32_t CDC_WS, uint32_t CDC_WE, uint32_t CDC_IE, uint32_t CDC_NP, uint32_t CDC_NP2, uint32_t CDC_PG, uint32_t CDC_H, uint32_t CDC_TH, uint32_t CDC_TL, uint32_t FDC_WS, uint32_t FDC_WE, uint32_t FDC_IE, uint32_t FDC_NP, uint32_t FDC_NP2, uint32_t FDC_PG, uint32_t FDC_H, uint32_t FDC_TH, uint32_t FDC_TL)
{
	// Since we would have to write a lot of "cdc_algos_data."'s, make another
	// reference that's smaller so we just have to write "d." instead
	fa125_algos_data_t &d = fa125_algos_data;

	bool cdchit = kFALSE;
	if ((rocid > 24)&&(rocid < 29)) cdchit = kTRUE;   //CDC has rocid 25 to 28

	if (cdchit) {

		d.WINDOW_START = CDC_WS;
		d.WINDOW_END = CDC_WE;
		d.INT_END = CDC_IE;

		d.NPED = CDC_NP;
		d.NPED2 = CDC_NP2;
		d.PG = CDC_PG;

		d.HIT_THRES = CDC_H;
		d.HIGH_THRESHOLD = CDC_TH;
		d.LOW_THRESHOLD = CDC_TL;

	} else {

		d.WINDOW_START = FDC_WS;
		d.WINDOW_END = FDC_WE;
		d.INT_END = FDC_IE;

		d.NPED = FDC_NP;
		d.NPED2 = FDC_NP2;
		d.PG = FDC_PG;

		d.HIT_THRES = FDC_H;
		d.HIGH_THRESHOLD = FDC_TH;
		d.LOW_THRESHOLD = FDC_TL;

	}


	if (samples.size()<=(uint32_t)d.WINDOW_END + 20) {
		static int Nwarn=0;
		if(Nwarn<10){
        	cout << "The number of samples passed into the fa125_algos routine (" << samples.size() << ") is less than the" << endl;
			cout << "minimum required by the parameters in use (" << d.WINDOW_END+21 << "). " << endl;
			cout << "Parameter WE (" << d.WINDOW_END << ") should be decreased to " << samples.size()-21 << " or less." << endl;
			if(++Nwarn==10) cout <<" --- LAST WARNING! ---" << endl;
		}
		return;
	}

	if (d.NPED2 > d.NPED) {
		static int Nwarn=0;
		if(Nwarn<10){
			cout << "Parameter NPED is too small or NPED2 is too large. " << endl;
			cout << "NPED (" << d.NPED << ") should be increased or NPED2 (" << d.NPED2 << ") decreased until NPED >= NPED2." << endl;
			if(++Nwarn==10) cout <<" --- LAST WARNING! ---" << endl;
		}
		return;
	}

    if (d.WINDOW_START < d.NPED) {
		static int Nwarn=0;
		if(Nwarn<10){
			cout << "Parameter WS (" << d.WINDOW_START << ") is too small or NPED (" << d.NPED << ") too large." << endl;
			cout << "WS should be >= NPED." << endl;
			if(++Nwarn==10) cout <<" --- LAST WARNING! ---" << endl;
		}
		return;
    }
	
	// Copy uint16_t samples into Int_t type array so we can pass it into the cdc_algos2
	// routine that does the actual work
	Int_t adc[d.WINDOW_END+21];
	for(uint32_t i=0; i<=(uint32_t)d.WINDOW_END+20; i++) adc[i] = (Int_t)samples[i];
	
	// Call the actual routine that does the heavy lifting
	fa125_algos(d.time, d.q_code, d.pedestal, d.integral, d.overflows, d.maxamp, adc, d.WINDOW_START, d.WINDOW_END, d.INT_END, d.NPED, d.NPED2, d.PG, d.HIT_THRES, d.HIGH_THRESHOLD, d.LOW_THRESHOLD);

}


void fa125_algos(Int_t &time, Int_t &q_code, Int_t &pedestal, Long_t &integral, Int_t &overflows, Int_t &maxamp, Int_t adc[], Int_t WINDOW_START, Int_t WINDOW_END, Int_t INT_END, Int_t NPED, Int_t NPED2, Int_t PG, Int_t HIT_THRES, Int_t HIGH_THRESHOLD, Int_t LOW_THRESHOLD) {


  const Int_t NU = 20;  //number of samples sent to time algo
  const Int_t PED = 5;  //sample to be used as pedestal for timing is in place 5

  const Int_t XTHR_SAMPLE = PED + PG; 


  Int_t adc_subset[NU]; 

  Int_t hitfound=0; //hit found or not (1=found,0=not)
  Int_t hitsample=-1;  // if hit found, sample number of threshold crossing
  Int_t timesample=0;
 
  Int_t i=0;

  time=0;       // hit time in 0.1xsamples since start of buffer passed to cdc_time
  q_code=-1;    // quality code, 0=good, 1=returned rough estimate
  pedestal=0;   // pedestal just before hit
  integral=0;   // signal integral, total
  overflows=0;  // count of samples with overflow bit set (need raw data, not possible from my root files)
  maxamp=0;     // signal amplitude at first max after hit


  // look for hit using mean pedestal of NPED samples before trigger 
  cdc_hit(hitfound, hitsample, pedestal, adc, WINDOW_START, WINDOW_END, HIT_THRES, NPED, NPED2, PG);


  if (hitfound==1) {

    for (i=0; i<NU; i++) {
      adc_subset[i] = adc[hitsample+i-XTHR_SAMPLE];
    }

    cdc_time(time, q_code, adc_subset, NU, PG, HIGH_THRESHOLD, LOW_THRESHOLD);

    timesample = hitsample-XTHR_SAMPLE + (Int_t)(0.1*time);  //sample number containing leading edge sample

    cdc_integral(integral, overflows, timesample, adc, WINDOW_END, INT_END);

    cdc_max(maxamp, hitsample, adc, WINDOW_END);

    time = 10*(hitsample-XTHR_SAMPLE) + time;   // integer number * 0.1 samples

  }

}




void cdc_hit(Int_t &hitfound, Int_t &hitsample, Int_t &pedestal, Int_t adc[], Int_t WINDOW_START, Int_t WINDOW_END, Int_t HIT_THRES, Int_t NPED, Int_t NPED2, Int_t PG) {


  pedestal=0;  //pedestal
  Int_t threshold=0;

  Int_t i=0;

  // calc pedestal as mean of NPED samples before trigger
  for (i=0; i<NPED; i++) {
    pedestal += adc[WINDOW_START-NPED+i];
  }

  pedestal = ( NPED==0 ? 0:(pedestal/NPED) );   // Integer div is ok as fpga will do 2 rightshifts

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

    //calculate new pedestal ending just before the hit

    pedestal = 0;

    for (i=0; i<NPED2; i++) {
      pedestal += adc[hitsample-PG-i];
    }

    pedestal = ( NPED2==0 ? 0:(pedestal/NPED2) );
  }


}




void cdc_integral(Long_t& integral, Int_t& overflows, Int_t timesample, Int_t adc[], Int_t WINDOW_END, Int_t INT_END) {

  Int_t i=0;

  integral = 0;
  overflows = 0;

  Int_t lastsample = timesample + INT_END - 1;

  if (lastsample > WINDOW_END) lastsample = WINDOW_END;

  for (i = timesample; i <= lastsample; i++ ) {

    integral += (Long_t)adc[i];
    if (adc[i]==(Long_t)4095) overflows++;    // only a placeholder at present. need to test on sample's overflow bit

  }


}




void cdc_max(Int_t& maxamp, Int_t hitsample, Int_t adc[], Int_t WINDOW_END) {

  int i;
  int ndec = 0;  //number of decreasing samples

  maxamp = adc[hitsample];

  for (i=hitsample; i<=WINDOW_END; i++) {

    if (adc[i] > adc[i-1]) {
      maxamp = adc[i];
      ndec = 0;
    }

    if (adc[i] <= adc[i-1]) ndec++;
    if (ndec==2) break;
  }


}




void cdc_time(Int_t &le_time, Int_t &q_code, Int_t adc[], Int_t NU, Int_t PG, Int_t THRES_HIGH, Int_t THRES_LOW) {


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

  Bool_t adczero = kFALSE;

  i = 0;
  while ((!adczero)&&(i<NU)) {

    if (adc[i] == 0) {
      adczero = kTRUE;
    }

    i++;
  }


  if (adczero) {

    le_time = ROUGH_TIME + 1;
    q_code = 1;
    return;
   
  }




  //check all samples from 0 to pedestal are <= LIMIT_PED_MAX

  Bool_t pedlimit = kFALSE;

  i = 0;
  while ((!pedlimit)&&(i<PED+1)) {

    if (adc[i] > LIMIT_PED_MAX) {
      pedlimit = kTRUE;
    }

    i++;
  }



  if (pedlimit) {

    le_time = ROUGH_TIME + 2;
    q_code = 2;
    return;
   
  }

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


void upsamplei(Int_t x[], Int_t startpos, Int_t z[], const Int_t NUPSAMPLED) {

  // x is array of samples
  // z is array of upsampled data
  // startpos is where to start upsampling in array x, only need to upsample a small region


  const Int_t nz = NUPSAMPLED;

  //  const Int_t Kscale = 32768;
  //  const Int_t K[43]={-8,-18,-27,-21,10,75,165,249,279,205,-2,-323,-673,-911,-873,-425,482,1773,3247,4618,5591,5943,5591,4618,3247,1773,482,-425,-873,-911,-673,-323,-2,205,279,249,165,75,10,-21,-27,-18,-8}; //32768

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

