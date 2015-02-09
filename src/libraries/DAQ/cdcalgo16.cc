
//
// $Id:$
//
//
// (See comments in cdcalgo16.h)
//
// The #defines below might normally be placed in the header file, but
// these have some generic names and are intended only for this CDC
// timing algorithm. They may also be replaced with JANA configuration
// parameters and/or CCDB constants.
//
// The main entry point here is the cdc_algos2 routine. The samples
// are passed in along with references to many different return values.
// 

#include <cdcalgo16.h>

// data buffer
#define NADCBUFFER 180   //number of samples in buffer from adc
// nadcbuffer should be >= WINDOW_END-XTHR_SAMPLE+NSAMPLES
// to ensure that there will be enough data to upsample near WINDOW_END

// hit search 
// run 1602 use WINDOW_START 46, WINDOW_END ~130, HIT_THRES ~300, LIMIT_PED_MAX ~200
// run 2016 use WINDOW_START 46, WINDOW_END ~130, HIT_THRES 100, LIMIT_PED_MAX 511
// run 2016 has a few channels firing often with value~300

#define WINDOW_START 46  //earliest sample for start of hit
#define WINDOW_END 99   //last sample for hit arrival, and last sample included in integral

#define HIT_THRES 300   // 5 sigma hit threshold
#define NPED 16         // number of samples used for pedestal used to find hit. must be 2**integer
#define NPED2 16        // number of samples used for pedestal calculated just before hit (returned, not used here). must be 2**integer


// cdc_time
#define HIGH_THRESHOLD 80   // 4 sigma threshold
#define LOW_THRESHOLD 20    // 1 sigma threshold
#define ROUGH_DT 24     // if pulse fails QA, return this many tenth-samples before threshold xing
#define INT_SAMPLE 6    // if pulse fails QA, start integration with this sample 

//cdc_time good pulse criteria
#define LIMIT_PED_MAX 511    //return rough time if any sample in 0 to PED_SAMPLE exceeds this
#define LIMIT_ADC_MAX 4095  // return rough time if any sample in PED_SAMPLE+1 to NSAMPLES exceeds this

// cdc_time upsampling
#define NSAMPLES 15     //number of samples in subset of array to pass to cdc_time
#define XTHR_SAMPLE 9   // the hit_thres xing sample is sample[9] passed into cdc_time, starting with sample[0]
#define PED_SAMPLE 5    // take local ped as sample[5] passed into cdc_time
//#define NUPSAMPLED 8       // number of upsampled values to calculate, minimum is 8 (moved to header so it can set size of iubuf array DL)
#define SET_ADC_MIN 20     // add an offset to the adc values to set the min value equal to SET_ADC_MIN
#define LIMIT_UPS_ERR 30   // upsampling error tolerance, return midpoint time if error is greater than this


//diagnostic output
#define VERBOSE_T 0   // set to 1/0 to switch cdc_time text output on/off
#define VERBOSE_U 0   // set to 1/0 to switch upsampling text output on/off
#define VERBOSE_S 0   // set to 1/0 to switch hit search output on/off
#define VERBOSE 0     // set to 1/0 to switch other text output on/off


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


// Wrapper for routine below
void cdc_algos2(vector<uint16_t> samples, cdc_algos_data_t &cdc_algos_data)
{
	if(samples.size()<=WINDOW_END){
		cout << "The number of samples passed into the cdc_algos2 routine is less than the" << endl;
		cout << "minimum (" << samples.size() << " <= " << WINDOW_END << "). The code is" << endl;
		cout << "currently not capable of handling this. " << endl;
		exit(-1);
	}
	
	// Copy uint16_t samples into Int_t type array so we can pass it into the cdc_algos2
	// routine that does the actual work
	Int_t adc[WINDOW_END+1];
	for(uint32_t i=0; i<=WINDOW_END; i++) adc[i] = (Int_t)samples[i];
	
	// Since we would have to write a lot of "cdc_algos_data."'s, make another
	// reference that's smaller so we just have to write "d." instead
	cdc_algos_data_t &d = cdc_algos_data;
	cdc_algos2(d.time, d.q_code, d.pedestal, d.integral, d.overflows, d.maxamp, d.initped, adc, d.iubuf, d.ups_err1, d.ups_err2, d.subsetmax, d.subsetmin, d.subsetstart);

}



void cdc_algos2(Int_t &time, Int_t &q_code, Int_t &pedestal, Long_t &integral, Int_t &overflows, Int_t &maxamp, Int_t &initped, Int_t adc[], Int_t iubuf[], Int_t &ups_err1, Int_t &ups_err2, Int_t &subsetmax, Int_t &subsetmin, Int_t &subsetstart) {

 
  Int_t adc_subset[NSAMPLES] = {0}; 

  Int_t hitfound=0; //hit found or not (1=found,0=not)
  Int_t hitsample=-1;  // if hit found, sample number of threshold crossing

  Float_t time_samples=0; // hit time in samples since start of adc buffer
  Float_t time_ns=0; // hit time in ns since start of adc buffer

  //  Int_t time=0;    // hit time in 0.1xsamples since start of buffer passed to cdc_time
  //  Int_t q_code=-1; //quality code, 0=good, 1=returned rough estimate (high threshold not met)

  //  Int_t pedestal=0;   // pedestal just before hit

  Int_t integral1=0;   // signal integral, from le time to threshold crossing
  Long_t integral2=0;   // signal integral, from threshold crossing to end of signal

  //  Int_t integral=0;   // signal integral, total
  //  Int_t overflows=0;  // count of samples with overflow bit set (need raw data, not possible from my root files)
  //  Int_t maxamp=0;     // signal amplitude at first max after hit

  Int_t i=0;

  Int_t pulselength=0; //duration of pulse, from sample following le_time to WINDOW_END

  cdc_hit(hitfound,hitsample,pedestal,adc,adc_subset);   // look for hit using mean pedestal of 16 samples before trigger 

  subsetmax = 0;
  subsetmin = 4095; 
  subsetstart = 0;

  if (hitfound==1) {

    for (i=0; i<NSAMPLES; i++) {
      adc_subset[i] = adc[hitsample+i-XTHR_SAMPLE];
      if (VERBOSE) printf("adc_subset[%i] %i\n",i,adc_subset[i]);
      if (adc_subset[i]<subsetmin) subsetmin = adc_subset[i];
    }

    subsetstart = hitsample-XTHR_SAMPLE;
    subsetmax = adc_subset[NSAMPLES-1];

    cdc_time(time, q_code, integral1, adc_subset, iubuf, ups_err1, ups_err2);

    cdc_integral(integral2, overflows, hitsample, adc);

    cdc_max(maxamp, hitsample, adc);

    time = 10*(hitsample-XTHR_SAMPLE) + time;   // integer number * 0.1 samples

    integral = (Long_t)integral1 + integral2;

    //repeat initial pedestal calc here for diagnostics as it is not supposed to be returned by cdc_hit
    for (i=0; i<NPED; i++) initped += adc[WINDOW_START-NPED+i];
    initped = initped/NPED;   // Integer div is ok as fpga will do 2 rightshifts



    if (VERBOSE) {
      time_samples = 0.1*time;
      time_ns = time_samples*8.0;

      pulselength = WINDOW_END - (Int_t)time_samples;

      printf ("Integral part 1 %i\n",integral1);
      printf ("Integral part 2 %li\n",integral2);
      printf ("Integral starts at sample %i, includes %i, ends at %i\n",(Int_t)time_samples,hitsample,WINDOW_END);

      printf("\n   Hit found...\n");
      printf("\n   Pedestal %i  (initial pedestal %i)\n",pedestal, initped);
      printf("\n   Max amplitude %i\n",maxamp);
      printf("\n   Integral %li - samples %i x pedestal %i = %li\n",integral,pulselength,pedestal,integral-(Long_t)(pulselength*pedestal));
      printf("\n   Overflow count %i\n",overflows);
      printf("\n   Time %5.1f samples (%4.1f ns) since start of ADC buffer\n",time_samples,time_ns);
      printf("\n   Q-code is %i \n\n",q_code);
      if (q_code>0) printf("\n*** This is a rough estimate ***\n\n");
    }
  } else { 

    //    if (VERBOSE && !VERBOSE_S) printf("\n   ADC values did not go over threshold \n\n");

  }

}




void cdc_hit(Int_t &hitfound, Int_t &hitsample, Int_t &ped, Int_t adc[], Int_t adc_subset[]){

  ped=0;  //pedestal
  Int_t threshold=0;

  Int_t i=0;

  // calc pedestal as mean of NPED samples before trigger

  for (i=0; i<NPED; i++) ped += adc[WINDOW_START-NPED+i];

  ped = ped/NPED;   // Integer div is ok as fpga will do 2 rightshifts

  threshold = ped + HIT_THRES;

  // look for threshold crossing

  if (VERBOSE) printf("Searching for ADC value >= %i\n",threshold);

  i = WINDOW_START - 1;
  hitfound = 0;

  while ((hitfound==0) && (i<WINDOW_END)) {

    i++;

    if (VERBOSE_S) printf("  ADC[%i] %i\n",i,adc[i]);

    if (adc[i] >= threshold) {

      hitfound = 1;
      hitsample = i;

    }
  }

  if (hitfound == 1) {

    //copy values to small array of samples to send to time module

    for (i=0;i<NSAMPLES;i++) {
      adc_subset[i] = adc[hitsample+i];
    }


    //calculate new pedestal ending just before the hit

    ped = 0;

    for (i=0; i<NPED2; i++) {
      ped += adc[hitsample-XTHR_SAMPLE+PED_SAMPLE-i];
    }

    ped = ped/NPED2;

    if (VERBOSE) printf("\n   ADC value %i at sample %i over threshold %i, pedestal %i\n\n",adc[hitsample],hitsample,threshold,ped);

  }


  if ((VERBOSE) && (hitfound==0)) {
    printf("\n   ADC values did not go over threshold %i \n\n",threshold);
  }


}




void cdc_integral(Long_t& integral, Int_t& overflows, Int_t hitsample, Int_t adc[]) {

  Int_t i=0;

  integral = 0;
  overflows = 0;


  for (i = hitsample; i <= WINDOW_END; i++ ) {

    integral += (Long_t)adc[i];
    if (adc[i]>4095) overflows++;    // only a placeholder at present. need to test on sample's overflow bit

  }


}




void cdc_max(Int_t& maxamp, Int_t hitsample, Int_t adc[]) {

  maxamp = 0;

  Int_t maxbin=0;

  maxbin = hitsample;

  while ( (adc[maxbin]<=adc[maxbin+1]) && (maxbin <= WINDOW_END ) ){
    maxbin++;
  }

  maxamp = adc[maxbin];

}






void cdc_time(Int_t &le_time, Int_t &q_code, Int_t &integral, Int_t adc[], Int_t iubuf[], Int_t &ups_err1, Int_t &ups_err2){


  // returned quantities:

  //  Int_t le_time = 0;  // leading edge time as 0.1 samples since time of first sample supplied
  //  Int_t q_code = 0;   // quality code, 0=good, 1=rough estimate as data did not exceed threshold

  // q_code list:
  //   0: Good
  //   1: ADC data did not go over threshold adc_thres_hi 
  //   2: Leading edge time is outside the upsampled region (cross adc_thres_lo too late in the buffer subset ) 
  //   3: Last upsampled point is <= low timing threshold
  //   4: Upsampled points did not go below low timing threshold
  //   5: ADC sample value of 0 found
  //   6: ADC sample value > LIMIT_ADC_MAX found
  //   7: Pedestal ADC[PED_SAMPLE] value > LIMIT_PED_MAX found
  //   8: Upsampled point is below 0
  //   9: Difference between upsampled and sampled values > LIMIT_UPS_ERR

  // Input:
  //  Int_t adc[NSAMPLES] = {65,62,56,46,41,56,85,109,120,122,127,150,181,197};
  //  defined as Int_t to save type casts later on


  // config constants, defined above as globals

  //#define NSAMPLES 15;    //number of samples to pass to cdc_time
  //#define XTHR_SAMPLE 9;  // the 5 sigma thres xing is sample[8] passed into cdc_time, starting with sample[0]
  //#define PED_SAMPLE 5;   // take local ped as sample[4] passed into cdc_time

  //#define HIGH_THRESHOLD 64;    // 4 sigma
  //#define LOW_THRESHOLD 16;    // 1 sigma
  //#define ROUGH_DT 20;    // if algo fails, return this many tenth-samples before threshold xing
  //#define NUPSAMPLED 8;   // number of upsampled values to calculate
  //#define INT_SAMPLE 6; // if algo fails start integration with this sample 


  // internal constants


  const Int_t START_SEARCH = PED_SAMPLE+1; // -- start looking for hi threshold xing with this sample

  const Int_t ROUGH_TIME = (XTHR_SAMPLE*10)-ROUGH_DT; // --return this for time if the algo fails

  //  Int_t iubuf[NUPSAMPLED] = {0};  // array of upsampled values
   
  //	-- iubuf 0 corresponds to 0.2 before low thres xing sample
  //	-- iubuf 1 maps to low thres xing sample


  Int_t adc_thres_hi = 0; // high threshold
  Int_t adc_thres_lo = 0; // low threshold

  //    -- contributions to hit time, these are summed together eventually, units of sample/10

  Int_t itime1 = 0; // which sample
  Int_t itime2 = 0; // which minisample
  Int_t itime3 = 0; // correction from interpolation
    

  //    -- search vars

  Int_t adc_sample_hi = 0; // integer range 0 to NSAMPLES := 0;  --sample number for adc val at or above hi thres
  Int_t adc_sample_lo = 0; // integer range 0 to NSAMPLES := 0;  -- sample num for adc val at or below lo thres
  Int_t adc_sample_lo2 = 0; // integer range 0 to 12:= 0;  -- minisample num for adc val at or below lo thres

  Bool_t over_threshold = kFALSE;
  Bool_t below_threshold = kFALSE;


  // interpolation vars

  Int_t denom = 0; 
  Int_t limit = 0;
  Int_t sum = 0;
  Int_t ifrac = 0;
  

  Int_t i = 0;

  ups_err1 = 0;
  ups_err2 = 0;


  //check all samples are >0

  Bool_t adczero = kFALSE;

  i = 0;
  while ((!adczero)&&(i<NSAMPLES)) {

    if (adc[i] == 0) {
      adczero = kTRUE;
    }

    i++;
  }


  if (adczero) {

    if (VERBOSE) printf ("ADC value 0 found\n");

    le_time = ROUGH_TIME;
    q_code = 5;

    integral = 0;    

    for (i=INT_SAMPLE; i<XTHR_SAMPLE; i++) integral += adc[i];

    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 1 - this is a rough estimate\n");

    return;
   
  }


  //check all samples are <= LIMIT_ADC_MAX

  Bool_t adclimit = kFALSE;

  i = 0;
  while ((!adclimit)&&(i<NSAMPLES)) {

    if (adc[i] > LIMIT_ADC_MAX) {
      adclimit = kTRUE;
    }

    i++;
  }


  if (adclimit) {

    if (VERBOSE) printf ("ADC value >= %i found\n",LIMIT_ADC_MAX);

    le_time = ROUGH_TIME;
    q_code = 6;

    integral = 0;    

    for (i=INT_SAMPLE; i<XTHR_SAMPLE; i++) integral += adc[i];

    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 1 - this is a rough estimate\n");

    return;
   
  }


  //check all samples from 0 to pedestal are <= LIMIT_PED_MAX

  Bool_t pedlimit = kFALSE;

  i = 0;
  while ((!pedlimit)&&(i<PED_SAMPLE+1)) {

    if (adc[i] > LIMIT_PED_MAX) {
      pedlimit = kTRUE;
    }

    i++;
  }



  if (pedlimit) {

    if (VERBOSE) printf ("Pedestal value >= %i found\n",LIMIT_PED_MAX);

    le_time = ROUGH_TIME;
    q_code = 7;

    integral = 0;    

    for (i=INT_SAMPLE; i<XTHR_SAMPLE; i++) integral += adc[i];

    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 1 - this is a rough estimate\n");

    return;
   
  }

  //  add offset to move min val in subset equal to SET_ADC_MIN
  //  this is to move samples away from 0 to avoid upsampled pts going -ve (on a curve betw 2 samples)

  Int_t adcmin = 4095; 

  i=0; 

  while (i<NSAMPLES) {

    if (adc[i] < adcmin) {
      adcmin = adc[i];
    }

    i++;
  }

  Int_t adcoffset = SET_ADC_MIN - adcmin;  

  i=0; 

  while (i<NSAMPLES) {
    adc[i] = adc[i] + adcoffset;
    i++;
  }

  if (VERBOSE_T) printf("adjusted adc values by %i\n",adcoffset);

  // eg if adcmin is 100, setmin is 30, adcoffset = 30 - 100 = -70, move adc down by 70
 

  //////////////////////////////

  // calc thresholds

  adc_thres_hi = adc[PED_SAMPLE] + HIGH_THRESHOLD;
  adc_thres_lo = adc[PED_SAMPLE] + LOW_THRESHOLD;

  if (VERBOSE_T) printf("thresholds: hi %i lo %i\n",adc_thres_hi,adc_thres_lo);

  // search for high threshold crossing

  over_threshold = kFALSE;
  i = START_SEARCH;

  while ((!over_threshold)&&(i<NSAMPLES)) {

    if (adc[i] >= adc_thres_hi) {
      adc_sample_hi = i;
      over_threshold = kTRUE;
    }

    i++;
  }


  if (!over_threshold) {

    if (VERBOSE) printf ("ADC data did not go over threshold %i\n",adc_thres_hi);

    le_time = ROUGH_TIME;
    q_code = 1;

    integral = 0;    

    for (i=INT_SAMPLE; i<XTHR_SAMPLE; i++) integral += adc[i];

    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 1 - this is a rough estimate\n");

    return;
   
  }

  if (VERBOSE_T) printf("threshold %i met or exceeded by sample %i value %i\n",adc_thres_hi,adc_sample_hi,adc[adc_sample_hi]);


  // search for low threshold crossing

  below_threshold = kFALSE;
  i = adc_sample_hi-1;

  while ((!below_threshold) && (i>=PED_SAMPLE)) {   //************changed START_SEARCH to PED_SAMPLE**********

    if (VERBOSE_T) printf("adc[%i] %i\n",i,adc[i]);

    if (adc[i] <= adc_thres_lo) {
      adc_sample_lo = i;
      itime1 = i*10;
      below_threshold = kTRUE;
    }

    i--;
  }

  if (adc_sample_lo > NSAMPLES-7) {

    if (VERBOSE) printf ("Leading edge time is outside the upsampled region\n");

    le_time = ROUGH_TIME;
    q_code = 2;

    integral = 0;    

    for (i=INT_SAMPLE; i<XTHR_SAMPLE; i++) integral += adc[i];

    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 1 - this is a rough estimate\n");

    return;
   
  }



  if (VERBOSE_T) printf("threshold %i met or preceded by sample %i value %i\n",adc_thres_lo,adc_sample_lo,adc[adc_sample_lo]);

  if (VERBOSE_T) printf("itime1 %i x 0.1 samples\n",itime1);
           

  //upsample values from adc_sample_lo - 0.2 to adc_sample_lo + 1.2 

  upsamplei(adc,adc_sample_lo,iubuf);



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

    if (VERBOSE) printf ("Negative upsampled value returned\n");

    le_time = itime1 + 5;   // better bail-out //midway between adc_sample_lo and _lo + 1

    q_code = 8;

 
    integral = 0;    

          
    for (i=PED_SAMPLE-1; i<XTHR_SAMPLE; i++) {
      if (i>adc_sample_lo) integral += adc[i];
    }

    //    for (i=INT_SAMPLE; i<XTHR_SAMPLE; i++) integral += adc[i];

    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 8 - this is a rough estimate\n");

    return;
   
  }




  // correct errors

  ups_err1 = iubuf[1] - adc[adc_sample_lo];
  ups_err2 = iubuf[6] - adc[adc_sample_lo+1];

  Int_t ups_err_sum = ups_err1 + ups_err2;

  Int_t ups_adjust = (Int_t)(0.5*ups_err_sum);

  // if this is more than set limit, bail out

  if (ups_err_sum > LIMIT_UPS_ERR) { //upsampled 

    if (VERBOSE) printf ("Upsampling error sum is over limit %i\n",LIMIT_UPS_ERR);

    le_time = itime1 + 5;   // better bail-out //midway between adc_sample_lo and _lo + 1

    q_code = 9;

    integral = 0;    
       
    for (i=PED_SAMPLE-1; i<XTHR_SAMPLE; i++) {
      if (i>adc_sample_lo) integral += adc[i];
    }


    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 9 - this is a rough estimate\n");

    return;


  }






  // move threshold correspondingly instead of correcting upsampled values

  adc_thres_lo = adc_thres_lo + ups_adjust;

  //iubuf(0) is at adc_sample_lo - 0.2
  if (VERBOSE) printf ("iubuf[1] %i err %i iubuf[6] %i err %i Correction is %i Moved threshold to %i\n",iubuf[1],ups_err1,iubuf[6],ups_err2,ups_adjust,adc_thres_lo);


  //search through upsampled array


  if (iubuf[NUPSAMPLED-1]<= adc_thres_lo) { //bad upsampling

    if (VERBOSE) printf ("Last upsampled point is <= low timing threshold\n");

    le_time = itime1 + 5;   //midway

    q_code = 3;

    integral = 0;    
       
    for (i=PED_SAMPLE-1; i<XTHR_SAMPLE; i++) {
      if (i>adc_sample_lo) integral += adc[i];
    }

    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 3 - this is a rough estimate\n");

    return;

  }





  below_threshold = kFALSE;
  i = NUPSAMPLED-2;

  while ((!below_threshold) && (i>=0)) {

    if (iubuf[i] <= adc_thres_lo) {
      adc_sample_lo2 = i;
      below_threshold = kTRUE;
    }

    i--;
  }


  /************ bug fix ************/
  /*** do not trust upsampling completely **/

  if (!below_threshold) { //upsampled points did not go below thres

    if (VERBOSE) printf ("Upsampled points did not go below low timing threshold\n");


    le_time = itime1 + 5;   //midway

    q_code = 4;

    integral = 0;    
       
    for (i=PED_SAMPLE-1; i<XTHR_SAMPLE; i++) {
      if (i>adc_sample_lo) integral += adc[i];
    }


    if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
    if (VERBOSE_T) printf("Quality code is 1 - this is a rough estimate\n");

    return;


  }





  if (VERBOSE_T) printf("threshold %i met or preceded by upsampled point %i value %i\n",adc_thres_lo,adc_sample_lo2,iubuf[adc_sample_lo2]);

  itime2 = adc_sample_lo2*2;  //  convert from sample/5 to sample/10

  if (VERBOSE_T) printf("itime2 %i x 0.1 samples\n",itime2);


  //interpolate

  itime3 = 0;
            
  if (iubuf[adc_sample_lo2] != adc_thres_lo) {
                        
    denom = iubuf[adc_sample_lo2+1] - iubuf[adc_sample_lo2];
    limit = (adc_thres_lo - iubuf[adc_sample_lo2])*2;

    if (VERBOSE_T) printf("denom %i limit %i\n",denom,limit);
  
    sum = 0;
    ifrac = 0;
                          
    while (sum<limit) {

      sum = sum + denom;
      ifrac = ifrac + 1;           

    }

    if (VERBOSE_T) printf("sum %i ifrac %i\n",sum,ifrac);

    if (2*(sum-limit) > denom) { //  --round
      itime3 = ifrac - 1;        
    } else {
      itime3 = ifrac;
    }

  }

  if (VERBOSE_T) printf("itime3 %i x 0.1 samples\n",itime3);

  //--              need to subtract 1/5 sample from itime2 because upsampling starts
  //--              1 minisample before adc_sample_lo
     
  le_time = itime1 - 2 + itime2 + itime3;  //   -- this is time from first sample point, in 1/10ths of samples
  q_code = 0;
  integral = 0;
          
  for (i=PED_SAMPLE-1; i<XTHR_SAMPLE; i++) {
    if (i>adc_sample_lo) integral += adc[i];
  }

  if (VERBOSE_T) printf("Leading edge time is %i x 0.1 samples since the time of the first sample supplied \n",le_time);
  if (VERBOSE_T) printf("Quality code is 0 - good\n");

}


void upsamplei(Int_t x[], Int_t startpos, Int_t z[]) {

  // x is array of samples
  // z is array of upsampled data
  // startpos is where to start upsampling in array x, only need to upsample a small region


  const Int_t nz = NUPSAMPLED;

  const Int_t K[43]={-8,-18,-27,-21,10,75,165,249,279,205,-2,-323,-673,-911,-873,-425,482,1773,3247,4618,5591,5943,5591,4618,3247,1773,482,-425,-873,-911,-673,-323,-2,205,279,249,165,75,10,-21,-27,-18,-8}; //32768

  Int_t k,j,dk;

  const Int_t Kscale = 32768;

  //don't need to calculate whole range possible
  //earliest value k=42 corresponds to sample 4.2
  //               k=43                sample 4.4
  //               k=46                sample 5.0
  
  // I need my first value to be for sample startpos - 0.2

  // sample 4 (if possible) would be at k=41
  // sample 4.2                         k=42
  // sample 5                           k=46

  // sample x                           k=41 + (x-4)*5
  // sample x-0.2                       k=40 + (x-4)*5


  Int_t firstk = 40 + (startpos-4)*5;


  for (k=firstk; k<firstk+nz; k++) {

    dk = k - firstk;    

    z[dk]=0.0;

    for (j=k%5;j<43;j+=5) {

      z[dk] += x[(k-j)/5]*K[j]; 

    }

    //    printf("dk %i z %i 5z %i  5z/scale %i\n",dk,z[dk],5.0*z[dk],5.0*z[dk]/Kscale);

    z[dk] = (Int_t)(5*z[dk])/Kscale;

    if (VERBOSE_U) printf("upsampled %i %i \n",dk,z[dk]);
    
  }

}
