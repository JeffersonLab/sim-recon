// $Id$
//
//    File: JEventProcessor_CDC_expert.cc
// Created: Wed Oct 22
// Creator: Naomi Jarvis


#include <stdint.h>
#include <vector>

#include <TMath.h>


#include "JEventProcessor_CDC_expert.h"
#include <JANA/JApplication.h>


using namespace std;
using namespace jana;


#include "CDC/DCDCHit.h"
#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125PulseIntegral.h"
#include "DAQ/Df125PulsePedestal.h"
#include "DAQ/Df125WindowRawData.h"
#include "DAQ/Df125CDCPulse.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>


// root hist pointers



static TH1D *cdc_e = NULL; 
static TH2D *cdc_e_vs_n = NULL; 

static TH1D *cdc_t = NULL; 
static TH2D *cdc_t_vs_n = NULL; 


static TH2D *cdc_e_ring[29];
static TH2D *cdc_t_ring[29];


static TH2D *cdc_e_vs_t; 
static TH2D *cdc_e_vs_t_ring[29];

static TH2I *cdc_raw_int_vs_t; 
static TH2I *cdc_raw_int_vs_t_ring[29];


static TH2I *cdc_o_badt;   
static TH2I *cdc_o_overflow;

static TH2I *cdc_ped_ring[29];  
static TH1I *cdc_ped_badt;  
static TH1I *cdc_ped_overflow;

static TH2I *cdc_windata_ped_ring[29];  

static TH2I *cdc_windata_ped_roc25;  
static TH2I *cdc_windata_ped_roc26;  
static TH2I *cdc_windata_ped_roc27;  
static TH2I *cdc_windata_ped_roc28;  

static TH2I *cdc_raw_t_ring[29];
static TH1I *cdc_raw_t_badt;
static TH1I *cdc_raw_t_overflow;

static TH2I *cdc_raw_amp_ring[29];  
static TH1I *cdc_raw_amp_badt;  
static TH1I *cdc_raw_amp_overflow;  



static TH2I *cdc_raw_intpp_ring[29];  
static TH1I *cdc_raw_intpp_badt;  
static TH1I *cdc_raw_intpp_overflow;  



static TH2I *cdc_raw_int_ring[29];  
static TH1I *cdc_raw_int_badt;  
static TH1I *cdc_raw_int_overflow;  




//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_CDC_expert());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_expert::JEventProcessor_CDC_expert() {
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_expert::~JEventProcessor_CDC_expert() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_CDC_expert::init(void) {

  // I moved all the histogram setup into the brun so that I can use different
  // scales for the later runs using the new firmware.  NSJ.

  //  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
  //  japp->RootUnLock(); //RELEASE ROOT LOCK!!


  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes



  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  // max values for histogram scales, modified fa250-format readout

  //  Int_t IMAX = 524288;  //max for raw integral, fa250-format, 19 bits
  Int_t IMAX = 400000;  //max for raw integral
  Int_t PMAX = 512;     //max for pedestal, fa250-format max is 512
  Int_t AMAX = 4096;    //max for amplitude, fa250-format, 12 bits
  Int_t RTMAX = 12000;  
  Int_t RTVSNMAX = 8192;  //raw time vs straw histogram range ends at this value

  Char_t rtunits[8] = "0.125ns";  //raw time is in units of sample/64 = ns/8

  if (runnumber > 3675) { //new fa125 format firmware, from 11 Sept 2015

    // raw quantities for read out (125 format) are
    //   time                    field max 2047   scaled x 1, units 0.8ns
    //   time qf                 field max 1 
    //   overflow count          field max 7
    //   pedestal                field max 255    scaled x 1/4 initially
    //   max amplitude 9 bits,   field max 511    scaled x 1/8
    //   integral                field max 16383  scaled x 1/14


    // max values for histogram scales, fa125-format readout

    IMAX = 16384; //max for raw integral, fa125-format, 14 bits
    PMAX = 256;   //max for pedestal, fa125-format, 8 bits
    AMAX = 512;    //max for amplitude, fa125-format, 9 bits
    RTMAX = 2048;  //max for raw time, fa125-format, 11 bits
    RTVSNMAX = 1024;  //raw time vs straw histogram range ends at this value

    sprintf(rtunits,"0.8ns");  //raw time is in units of sample/10 = 0.8ns

  }
  

  const Int_t EMAX = 21000;  //max for E histograms, fC
  //  const Int_t EMAX = 21000000;  //max for E histograms, fC
  // E histograms filled with a_scale*gains*(integration-pedestal)

  const Int_t TMAX = 2000;    //max for t histograms, ns
  // t histograms filled with t_scale*(raw-t - offset) + tmin


  const Int_t NSTRAWS = 3522;
  const Float_t HALF = 0.5;
  const Float_t NSTRAWSPH = 3522.5;


  // create root folder for cdc and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("CDC_expert")->cd();

 

  // book histograms

  //number of straws in each ring, starts with 0 so that straws[1] is the number of straws in ring 1
  const Int_t straws[29] = {0,42,42,54,54,66,66,80,80,93,93,106,106,123,123,135,135,146,146,158,158,170,170,182,182,197,197,209,209};


  cdc_e = new TH1D("cdc_e","CDC charge (fC);charge (fC)",200,0,EMAX);
  cdc_e_vs_n = new TH2D("cdc_e_vs_n","CDC charge (fC) vs straw number;straw;charge (fC)",NSTRAWS,HALF,NSTRAWSPH,100,0,EMAX);

  cdc_e_vs_t = new TH2D("cdc_e_vs_t","CDC charge (fC) vs time (ns);time (ns);charge (fC)",150,-250,TMAX,100,0,EMAX);


  cdc_t = new TH1D("cdc_t","CDC time (ns);time (ns)",300,-250,TMAX);
  cdc_t_vs_n = new TH2D("cdc_t_vs_n","CDC time (ns) vs straw number;straw;time (ns)",NSTRAWS,HALF,NSTRAWSPH,150,-250,TMAX);




  cdc_raw_int_vs_t   = new TH2I("cdc_raw_int_vs_t",Form("CDC integral (ADC units), pedestal subtracted, vs raw time (units of %s);time (%s);integral, pedestal subtracted (ADC units)",rtunits,rtunits),(Int_t)256,0,RTVSNMAX,100,0,IMAX);

  

  cdc_windata_ped_roc25   = new TH2I("cdc_windata_ped_roc25","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 25;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_windata_ped_roc26   = new TH2I("cdc_windata_ped_roc26","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 26;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_windata_ped_roc27   = new TH2I("cdc_windata_ped_roc27","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 27;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_windata_ped_roc28   = new TH2I("cdc_windata_ped_roc28","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 28;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);



  if (runnumber > 3675) { //new fa125 format firmware, from 11 Sept 2015

    gDirectory->mkdir("bad_t","CDC Bad time flagged")->cd();

    cdc_o_badt     = new TH2I("cdc_o_badt","CDC occupancy by straw,ring, events with bad time flagged;straw;ring",209,0.5,209.5,28,0.5,28.5);
    cdc_ped_badt   = new TH1I("cdc_ped_badt","CDC pedestal, events with bad time flagged;straw;pedestal",256,0,PMAX);
    cdc_raw_t_badt = new TH1I("cdc_raw_t_badt",Form("CDC raw time (units of %s), events with bad time flagged;straw;raw time (%s)",rtunits,rtunits),256,0,RTMAX);
    cdc_raw_amp_badt   = new TH1I("cdc_raw_amp_badt","CDC amplitude (ADC units), events with bad time flagged;ADC units",256,0,AMAX);
    cdc_raw_int_badt   = new TH1I("cdc_raw_intpp_badt","CDC integral (ADC units), pedestal subtracted, events with bad time flagged;ADC units",100,0,IMAX);
    cdc_raw_intpp_badt   = new TH1I("cdc_raw_intpp_badt","CDC integral (ADC units), including pedestal, events with bad time flagged;ADC units",100,0,IMAX);

    gDirectory->cd("../");
    gDirectory->mkdir("overflows","CDC overflow flagged")->cd();

    cdc_o_overflow = new TH2I("cdc_o_overflow","CDC overflow occupancy by straw,ring;straw;ring",209,0.5,209.5,28,0.5,28.5);
    cdc_ped_overflow  = new TH1I("cdc_ped_overflow","CDC pedestal, events with ADC overflow;pedestal",256,0,PMAX);
    cdc_raw_t_overflow = new TH1I("cdc_raw_t_overflow",Form("CDC raw time (units of %s), events with ADC overflow;raw time (%s)",rtunits,rtunits),256,0,RTMAX);
    cdc_raw_amp_overflow   = new TH1I("cdc_raw_amp_overflow","CDC amplitude (ADC units), events with ADC overflow;ADC units",256,0,AMAX);
    cdc_raw_int_overflow   = new TH1I("cdc_raw_intpp_overflow","CDC integral (ADC units), pedestal subtracted, events with ADC overflow;ADC units",100,0,IMAX);
    cdc_raw_intpp_overflow   = new TH1I("cdc_raw_intpp_overflow","CDC integral (ADC units), including pedestal, events with ADC overflow;ADC units",100,0,IMAX);

    gDirectory->cd("../");

  }  

  Int_t i;


  gDirectory->mkdir("rings_e_vs_t","CDC rings: charge vs time")->cd();
  
  for (i=1; i<29; i++) {
    cdc_e_vs_t_ring[i] = new TH2D(Form("cdc_e_vs_t_ring[%i]",i),"CDC charge (fC) vs time (ns);time (ns);charge (fC)",150,0,TMAX,100,0,EMAX);
  }

  gDirectory->cd("../");


  gDirectory->mkdir("rings_int_vs_raw_t","CDC rings: integral vs raw time (pedestal subtracted)")->cd();
  
  for (i=1; i<29; i++) {
    cdc_raw_int_vs_t_ring[i] = new TH2I(Form("cdc_raw_int_vs_t_ring[%i]",i),Form("CDC integral (ADC units), pedestal subtracted, vs raw time (%s);raw time (%s);integral, pedestal subtracted (ADC units)",rtunits,rtunits),256,0,RTVSNMAX,100,0,IMAX);

  }

  gDirectory->cd("../");



  gDirectory->mkdir("rings_e","CDC rings: charge vs straw")->cd();
  
  for (i=1; i<29; i++) {
    cdc_e_ring[i] = new TH2D(Form("cdc_e_ring[%i]",i),Form("CDC charge (fC), ring %i;straw;charge (fC)",i),straws[i],HALF,straws[i]+HALF,100,0,EMAX);
  }

  gDirectory->cd("../");
  gDirectory->mkdir("rings__t","CDC rings: time vs straw")->cd();
  
  for (i=1; i<29; i++) {
    cdc_t_ring[i] = new TH2D(Form("cdc_t_ring[%i]",i),Form("CDC time (ns), ring %i;straw;time (ns)",i),straws[i],HALF,straws[i]+HALF,150,0,TMAX);
  }


  gDirectory->cd("../");
  gDirectory->mkdir("rings_pedestal","CDC rings: pedestal vs straw")->cd();
  
  for (i=1; i<29; i++) {
    cdc_ped_ring[i]   = new TH2I(Form("cdc_ped_ring[%i]",i),Form("CDC pedestal (ADC units), ring %i;straw;pedestal",i),straws[i],HALF,straws[i]+HALF,(Int_t)PMAX/2,0,PMAX);
  }

  gDirectory->cd("../");
  gDirectory->mkdir("rings_windata_pedestal","CDC rings: pedestal from raw window data vs straw")->cd();
  
  for (i=1; i<29; i++) {
    cdc_windata_ped_ring[i]   = new TH2I(Form("cdc_windata_ped_ring[%i]",i),Form("CDC pedestal (ADC units) from raw window data, ring %i;straw;pedestal",i),straws[i],HALF,straws[i]+HALF,(Int_t)PMAX/2,0,PMAX);
  }





  gDirectory->cd("../");
  gDirectory->mkdir("rings_raw_t","CDC rings: raw time vs straw")->cd();

  for (i=1; i<29; i++) {
    cdc_raw_t_ring[i] = new TH2I(Form("cdc_raw_t_ring[%i]",i),Form("CDC raw time (units of %s), ring %i;straw;raw time (%s)",rtunits,i,rtunits),straws[i],HALF,straws[i]+HALF,256,0,RTVSNMAX);
  }

    
  gDirectory->cd("../");
  gDirectory->mkdir("rings_raw_amp","CDC rings: amplitude")->cd();

  for (i=1; i<29; i++) {
    cdc_raw_amp_ring[i]   = new TH2I(Form("cdc_raw_amp_ring[%i]",i),Form("CDC amplitude (ADC units), ring %i",i),straws[i],HALF,straws[i]+HALF,256,0,AMAX);
  }
  


  gDirectory->cd("../");
  gDirectory->mkdir("rings_raw_integral","CDC rings: integral vs straw (pedestal subtracted)")->cd();

  for (i=1; i<29; i++) {
    cdc_raw_int_ring[i]   = new TH2I(Form("cdc_raw_int_ring[%i]",i),Form("CDC integral (ADC units), pedestal subtracted, ring %i",i),straws[i],HALF,straws[i]+HALF,100,0,IMAX);
  }


  gDirectory->cd("../");
  gDirectory->mkdir("rings_raw_integral_incl_ped","CDC rings: integral vs straw (including pedestal)")->cd();

  for (i=1; i<29; i++) {
    cdc_raw_intpp_ring[i]   = new TH2I(Form("cdc_raw_intpp_ring[%i]",i),Form("CDC integral (ADC units), including pedestal, ring %i",i),straws[i],HALF,straws[i]+HALF,100,0,IMAX);
  }



  // back to main dir
  main->cd();

  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  float q,t;         // dcdchits quantities charge, time

  uint32_t qf,ocount; // time quality factor and overflow count from new firmware
  uint32_t tr,p,a; // dcdcdigihits raw quantities: time, pedestal, amplitude, quality factor, overflow count
  uint32_t integral; // dcdcdigihits integral, includes pedestal
  uint32_t integ;    // dcdcdigihits integral minus pedestal

  uint16_t ring,straw; // ring and straw numbers from either dcdchits or dcdcdigihits
  uint16_t n;         // straw number, 1 to 3522

  Bool_t PED_SUB;  // if this is false, integration window info is missing, so don't plot integrals

  uint32_t total_ped;  //total pedestal during integration period
  uint32_t nsamples_integral;    ///< number of samples used in integral 
  uint32_t nsamples_pedestal;    ///< number of samples used in pedestal

  uint32_t rocid;
  uint32_t slot;
  uint32_t channel;


  const uint16_t NPEDSAMPLES=16;

  //add extra 0 at front to use offset[1] for ring 1
  int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};


  //first set of histograms is for dcdchits, these are t and q after calibration
  //second set is for dcdcdigihits, these are the raw quantities

  // get hit data for cdc
  vector<const DCDCHit*> hits;
  eventLoop->Get(hits);
 
  // get raw data for cdc
  vector<const DCDCDigiHit*> digihits;
  eventLoop->Get(digihits);

  //get WRD data for new format (until it is linked to CDCPulse)
  vector<const Df125WindowRawData*> wrdvector;
  eventLoop->Get(wrdvector);

  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!


  for(uint32_t i=0; i<hits.size(); i++) {

    const DCDCHit *hit = hits[i];  // avoids having to use the uglier “cdcdigihits[0]->” syntax
    
    if(hit->q>0.0) {

      q      = hit->q;     // in fC
      t      = hit->t;      // in nanoseconds
      ring   = hit->ring;
      straw  = hit->straw;
      
      n = straw_offset[ring] + straw;
    
      if (q > 0.0) {
        cdc_e->Fill(q);
        cdc_e_vs_n->Fill(n,q);    
      }
    
      //if (t > 0.0) {
        cdc_t->Fill(t);
        cdc_t_vs_n->Fill(n,t);    
      //}

      cdc_e_vs_t->Fill(t,q);
      cdc_e_vs_t_ring[ring]->Fill(t,q);

      cdc_e_ring[ring]->Fill(straw,q);
      cdc_t_ring[ring]->Fill(straw,t);
    }
  }


  for(uint32_t i=0; i<digihits.size(); i++) {

    const DCDCDigiHit *digihit = digihits[i];  // avoids having to use the uglier “cdcdigihits[0]->” syntax

    // Get pointers to the underlying objects of interest
    const Df125PulseIntegral *pi = NULL;
    const Df125PulsePedestal *pp = NULL;
    const Df125WindowRawData *windat = NULL;
    const Df125CDCPulse      *cp = NULL;


    vector<uint16_t> samples;
    uint32_t winped=0;
    
    PED_SUB = kFALSE;  //set this to true when we find the config params
    total_ped = 0;
    
    rocid = 0;
    slot = 0;
    channel = 0;
    qf = 0;
    ocount = 0;
    a = 0;

    //old firmware uses Df125PulseIntegral and Df125PulsePedestal
    digihit->GetSingle(pi);
    if (pi) {
      rocid = pi->rocid;
      slot = pi->slot;
      channel = pi->channel;
      pi->GetSingle(windat);
    } else if (i < (uint32_t)wrdvector.size()) { 
      windat = wrdvector[i];
    }

    nsamples_integral = pi ? pi->nsamples_integral : 0;
    nsamples_pedestal = pi ? pi->nsamples_pedestal : 0;

    if ((nsamples_integral > 0) && (nsamples_pedestal > 0)) PED_SUB = kTRUE;

    digihit->GetSingle(pp);
    if(pp) a = pp->pulse_peak;

    //new firmware uses Df125CDCPulseData
    digihit->GetSingle(cp);
    if (cp) {
      rocid = cp->rocid;
      slot = cp->slot;
      channel = cp->channel;
      a = cp->first_max_amp;
      qf = cp->time_quality_bit;
      ocount = cp->overflow_count;
    }



    if (windat) {
 
      if (windat->samples.size()>=NPEDSAMPLES) {

        winped = 0;

        for (uint16_t j=0; j<NPEDSAMPLES; j++) winped += (uint32_t)windat->samples[j]; 
          
        winped = (uint32_t)winped/16.0;

        if (winped > 0) {

          if (rocid == 25) cdc_windata_ped_roc25->Fill(100*slot + channel,winped);
          if (rocid == 26) cdc_windata_ped_roc26->Fill(100*slot + channel,winped);
          if (rocid == 27) cdc_windata_ped_roc27->Fill(100*slot + channel,winped);
          if (rocid == 28) cdc_windata_ped_roc28->Fill(100*slot + channel,winped);

        }

      }//sample size
    } //windat




    if((digihit->pulse_integral>0)||(digihit->pulse_time>0)) {

      ring     = digihit->ring;
      straw    = digihit->straw;

      p        = digihit->pedestal;   
      tr       = digihit->pulse_time;    // raw time in 0.8 ns units
      integral = digihit->pulse_integral; // pulse integral in fadc units, pedestal not subtracted

      integ = 0;

      //ok to use p for pedestal subtraction here because if fa250 algo fails with p=0, integral=0 and amplitude=0 also

      if (PED_SUB) {
        total_ped = p*nsamples_integral/nsamples_pedestal;
        integ = integral - total_ped;
      }
             
      n = straw_offset[ring] + straw;
       
      if (PED_SUB) cdc_raw_int_vs_t->Fill(tr,integ);
      if (PED_SUB) cdc_raw_int_vs_t_ring[ring]->Fill(tr,integ);
      
      cdc_ped_ring[ring]->Fill(straw,p);
      cdc_raw_t_ring[ring]->Fill(straw,tr);
      cdc_raw_amp_ring[ring]->Fill(straw,a);   //no ped subtraction in case scaling factors differ
      if (PED_SUB) cdc_raw_int_ring[ring]->Fill(straw,integ);
      cdc_raw_intpp_ring[ring]->Fill(straw,integral);
      if (winped) cdc_windata_ped_ring[ring]->Fill(straw,winped);
      

      if (cp) {
	if (qf==1) {  // rough time flag is set
	  cdc_o_badt->Fill(straw,ring);
	  cdc_ped_badt->Fill(p);
	  cdc_raw_t_badt->Fill(tr);
	  cdc_raw_amp_badt->Fill(a); 
	  if (PED_SUB) cdc_raw_int_badt->Fill(integ); 
	  cdc_raw_intpp_badt->Fill(integral); 
	} 

      
	if (ocount>0) {  // overflow samples present
	  cdc_o_overflow->Fill(straw,ring);
	  cdc_ped_overflow->Fill(p);
	  cdc_raw_t_overflow->Fill(tr); 
	  cdc_raw_amp_overflow->Fill(a); 
	  if (PED_SUB) cdc_raw_int_overflow->Fill(integ);   
	  cdc_raw_intpp_overflow->Fill(integral);   
	} 
      }      


    }

  }


  japp->RootUnLock(); //RELEASE ROOT LOCK!!


  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
