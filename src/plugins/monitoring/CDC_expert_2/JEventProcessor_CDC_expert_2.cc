// $Id$
//
//    File: JEventProcessor_CDC_expert_2.cc
// Created: 26 Feb 2016
// Creator: Naomi Jarvis


#include <stdint.h>
#include <vector>

#include <TMath.h>


#include "JEventProcessor_CDC_expert_2.h"
#include <JANA/JApplication.h>


using namespace std;
using namespace jana;


#include "CDC/DCDCHit.h"
#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125WindowRawData.h"
#include "DAQ/Df125CDCPulse.h"
#include "DAQ/Df125Config.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>


// root hist pointers



static TH1D *cdc_e = NULL; 
static TH2D *cdc_e_vs_n = NULL; 

static TH1D *cdc_t = NULL; 
static TH2D *cdc_t_vs_n = NULL; 

static TH1D *cdc_rt = NULL; 
static TH2D *cdc_rt_vs_n = NULL; 

static TH1D *cdc_amp = NULL; 
static TH2D *cdc_amp_vs_n = NULL; 

static TH1D *cdc_rt_qf0 = NULL;

static TH1D *cdc_qf = NULL;
static TH2D *cdc_qf_vs_n = NULL;
static TH2D *cdc_qf_vs_a = NULL;
static TH2D *cdc_qf_vs_rt = NULL;


//static TH2D *cdc_e_ring[29];
static TH2D *cdc_t_ring[29];

static TH2D *cdc_e_vs_t; 
static TH2D *cdc_e_vs_t_ring[29];

static TH2I *cdc_int_vs_raw_t; 
static TH2I *cdc_int_vs_raw_t_ring[29];

static TH2I *cdc_o_overflow;
static TH1I *cdc_ped_overflow;
static TH1I *cdc_raw_t_overflow;

static TH2I *cdc_o_badt;   

static TH2I *cdc_ped_ring[29];  
static TH1I *cdc_ped_badt;  

static TH2I *cdc_raw_t_ring[29];
static TH1I *cdc_raw_t_badt;

static TH2I *cdc_amp_ring[29];  
static TH1I *cdc_amp_badt;  


static TH2I *cdc_intpp_ring[29];  
//static TH2I *cdc_int_ring[29];  

static TH2I *cdc_initped_ring[29];  

static TH2I *cdc_initped_roc25;  
static TH2I *cdc_initped_roc26;  
static TH2I *cdc_initped_roc27;  
static TH2I *cdc_initped_roc28;  

static TH2I *cdc_ped_roc25;  
static TH2I *cdc_ped_roc26;  
static TH2I *cdc_ped_roc27;  
static TH2I *cdc_ped_roc28;  




//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_CDC_expert_2());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_expert_2::JEventProcessor_CDC_expert_2() {
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_expert_2::~JEventProcessor_CDC_expert_2() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_CDC_expert_2::init(void) {


  // raw quantities for read out (fa125 new format) are
  //   time                    field max 2047   scaled x 1, units 0.8ns
  //   time qf                 field max 1 
  //   overflow count          field max 7
  //   pedestal                field max 255    scaled x 1/1
  //   max amplitude 9 bits    field max 511    scaled x 1/8
  //   integral                field max 16383  scaled x 1/16


  // max values for histogram scales

  const Int_t IMAX = 100000;   //max for raw integral
  const Int_t IPPMAX = 150000; //max for raw integral + pedestal
  const Int_t PMAX = 256;      //max for pedestal, fa125-format, 8 bits
  const Int_t AMAX = 4096;    //max for amplitude, fa125-format, 9 bits  * scale factor
  const Int_t RTMAX = 1800;    //max for raw time
  const Int_t RTMIN = 160;
  const Int_t RTBINS = 164;       //bins

  const Int_t RTVSNMAX = 1024;  //raw time vs straw histogram range ends at this value
 
  const Int_t EMAX = 21000;  //max for E histograms, fC
  //  const Int_t EMAX = 21000000;  //max for E histograms, fC
  // E histograms filled with a_scale*gains*(integration-pedestal)

  const Int_t TMAX = 1250;    //max for t histograms, ns
  // t histograms filled with t_scale*(raw-t - offset) + tmin
  const Int_t TMIN = -250;
  const Int_t TBINS = 250;


  const Int_t NSTRAWS = 3522;
  const Float_t HALF = 0.5;
  const Float_t NSTRAWSPH = 3522.5;

  //dead straws: K39 (row 11) and W38 (ring 23)

  Char_t deadstraws[32] = "(#709 and #2384 disconnected)";  //dead
  Char_t deadrow11[30] = "(#39 disconnected)";  //dead
  Char_t deadrow23[30] = "(#38 disconnected)";  //dead


  // create root folder for cdc and cd to it, store main dir
  TDirectory *main = gDirectory;

  gDirectory->mkdir("CDC_expert_2")->cd();
  TDirectory *xd = gDirectory;
 

  // book histograms

  //number of straws in each ring, starts with 0 so that straws[1] is the number of straws in ring 1
  const Int_t straws[29] = {0,42,42,54,54,66,66,80,80,93,93,106,106,123,123,135,135,146,146,158,158,170,170,182,182,197,197,209,209};


  cdc_e = new TH1D("cdc_e","CDC charge (fC);charge (fC)",200,0,EMAX);
  cdc_e_vs_n = new TH2D("cdc_e_vs_n",Form("CDC charge (fC) vs straw number;straw %s;charge (fC)",deadstraws),NSTRAWS,HALF,NSTRAWSPH,100,0,EMAX);

  cdc_e_vs_t = new TH2D("cdc_e_vs_t","CDC charge (fC) vs time (ns);time (ns);charge (fC)",TBINS,TMIN,TMAX,100,0,EMAX);


  cdc_t = new TH1D("cdc_t","CDC time (ns);time (ns)",TBINS*2,TMIN,TMAX);
  cdc_t_vs_n = new TH2D("cdc_t_vs_n",Form("CDC time (ns) vs straw number;straw %s;time (ns)",deadstraws),NSTRAWS,HALF,NSTRAWSPH,TBINS,TMIN,TMAX);


  cdc_rt = new TH1D("cdc_rt","CDC raw time (0.8ns);raw time (0.8ns)",RTMAX-RTMIN,RTMIN,RTMAX);
  cdc_rt_qf0 = new TH1D("cdc_rt_qf0","CDC raw time with qf=0 (0.8ns);raw time (0.8ns)",RTMAX-RTMIN,RTMIN,RTMAX);
  cdc_rt_vs_n = new TH2D("cdc_rt_vs_n",Form("CDC raw time (0.8ns) vs straw number;straw %s;raw time (0.8ns)",deadstraws),NSTRAWS,HALF,NSTRAWSPH,RTBINS,RTMIN,RTMAX);


  cdc_amp = new TH1D("cdc_amp","CDC amplitude;amplitude",256,0,AMAX);
  cdc_amp_vs_n = new TH2D("cdc_amp_vs_n",Form("CDC time (ns) vs straw number;straw %s;time (ns)",deadstraws),NSTRAWS,HALF,NSTRAWSPH,256,0,AMAX);






  cdc_int_vs_raw_t   = new TH2I("cdc_int_vs_raw_t",Form("CDC integral (ADC units), pedestal subtracted, vs raw time (0.8ns);raw time (0.8ns);integral, pedestal subtracted (ADC units)"),(Int_t)RTBINS,RTMIN,RTMAX,100,0,IMAX);
 

  gDirectory->mkdir("pedestals_by_roc","CDC Pedestals for each ROC")->cd();

  cdc_initped_roc25   = new TH2I("cdc_initped_roc25","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 25;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_initped_roc26   = new TH2I("cdc_initped_roc26","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 26;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_initped_roc27   = new TH2I("cdc_initped_roc27","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 27;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_initped_roc28   = new TH2I("cdc_initped_roc28","CDC pedestal (ADC units) from raw window data vs slot*100+channel, ROC 28;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);

  cdc_ped_roc25   = new TH2I("cdc_ped_roc25","CDC pedestal (ADC units) vs slot*100+channel, ROC 25;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_ped_roc26   = new TH2I("cdc_ped_roc26","CDC pedestal (ADC units) vs slot*100+channel, ROC 26;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_ped_roc27   = new TH2I("cdc_ped_roc27","CDC pedestal (ADC units) vs slot*100+channel, ROC 27;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);
  cdc_ped_roc28   = new TH2I("cdc_ped_roc28","CDC pedestal (ADC units) vs slot*100+channel, ROC 28;slot*100 + channel;pedestal",1600,200+HALF,1800+HALF,(Int_t)PMAX/4,0,PMAX);

  xd->cd();


  gDirectory->mkdir("bad_t","CDC Bad time flagged")->cd();

  cdc_o_badt     = new TH2I("cdc_o_badt","CDC occupancy by straw,ring, events with bad time flagged;straw;ring",209,0.5,209.5,28,0.5,28.5);
  cdc_ped_badt   = new TH1I("cdc_ped_badt","CDC pedestal, events with bad time flagged;straw;pedestal",256,0,PMAX);
  cdc_raw_t_badt = new TH1I("cdc_raw_t_badt",Form("CDC raw time (0.8ns), events with bad time flagged;straw;raw time (0.8ns)"),RTBINS,RTMIN,RTMAX);
  cdc_amp_badt   = new TH1I("cdc_amp_badt","CDC amplitude (ADC units), events with bad time flagged;ADC units",256,0,AMAX);
  //cdc_int_badt   = new TH1I("cdc_int_badt","CDC integral (ADC units), pedestal subtracted, events with bad time flagged;ADC units",100,0,IMAX);
  //cdc_intpp_badt   = new TH1I("cdc_intpp_badt","CDC integral (ADC units), including pedestal, events with bad time flagged;ADC units",128,0,IPPMAX);

  cdc_qf = new TH1D("cdc_qf","CDC time quality factor;time quality factor (0:good, 1:zero, 2:hi ped, 3: below TH, 4:late TCL, 5: neg ups, 9:hi ups)",10,0,10);
  cdc_qf_vs_n = new TH2D("cdc_qf_vs_n","CDC time quality factor vs straw number;straw;time quality factor",NSTRAWS,HALF,NSTRAWSPH,10,0,10);
  cdc_qf_vs_a = new TH2D("cdc_qf_vs_a","CDC time quality factor vs amplitude;amplitude;time quality factor",128,0,AMAX,10,0,10);
  cdc_qf_vs_rt = new TH2D("cdc_qf_vs_raw_t","CDC time quality factor vs raw time;time;time quality factor",RTBINS,RTMIN,RTMAX,10,0,10);

  xd->cd();


  gDirectory->mkdir("overflows","CDC overflow flagged")->cd();

  cdc_o_overflow = new TH2I("cdc_o_overflow","CDC overflow occupancy by straw,ring;straw;ring",209,0.5,209.5,28,0.5,28.5);
  cdc_ped_overflow  = new TH1I("cdc_ped_overflow","CDC pedestal, events with ADC overflow;pedestal",256,0,PMAX);
  cdc_raw_t_overflow = new TH1I("cdc_raw_t_overflow",Form("CDC raw time (0.8ns), events with ADC overflow;raw time (0.8ns)"),RTBINS,RTMIN,RTMAX);

  xd->cd();


  Int_t i;

  gDirectory->mkdir("rings_e_vs_t","CDC rings: charge vs time")->cd();
  for (i=1; i<29; i++) {
    cdc_e_vs_t_ring[i] = new TH2D(Form("cdc_e_vs_t_ring[%i]",i),"CDC charge (fC) vs time (ns);time (ns);charge (fC)",TBINS,TMIN,TMAX,100,0,EMAX);
  }
  xd->cd();


  gDirectory->mkdir("rings_int_vs_raw_t","CDC rings: integral vs raw time (pedestal subtracted)")->cd();
  for (i=1; i<29; i++) {
    cdc_int_vs_raw_t_ring[i] = new TH2I(Form("cdc_int_vs_raw_t_ring[%i]",i),Form("CDC integral (ADC units), pedestal subtracted, vs raw time (0.8ns);raw time (0.8ns);integral, pedestal subtracted (ADC units)"),RTBINS,RTMIN,RTMAX,100,0,IMAX);
  }
  xd->cd();


  // gDirectory->mkdir("rings_e","CDC rings: charge vs straw")->cd();
  // for (i=1; i<29; i++) {
  //   cdc_e_ring[i] = new TH2D(Form("cdc_e_ring[%i]",i),Form("CDC charge (fC), ring %i;straw;charge (fC)",i),straws[i],HALF,straws[i]+HALF,100,0,EMAX);
  // }
  // xd->cd();


  gDirectory->mkdir("rings_t","CDC rings: time vs straw")->cd();
  for (i=1; i<29; i++) {
    cdc_t_ring[i] = new TH2D(Form("cdc_t_ring[%i]",i),Form("CDC time (ns), ring %i;straw;time (ns)",i),straws[i],HALF,straws[i]+HALF,TBINS,TMIN,TMAX);
  }
  cdc_t_ring[11]->GetXaxis()->SetTitle(Form("pedestal %s",deadrow11));
  cdc_t_ring[23]->GetXaxis()->SetTitle(Form("pedestal %s",deadrow23));
  xd->cd();


  gDirectory->mkdir("rings_pedestal","CDC rings: pedestal vs straw")->cd();
  for (i=1; i<29; i++) {
    cdc_ped_ring[i]   = new TH2I(Form("cdc_ped_ring[%i]",i),Form("CDC pedestal (ADC units), ring %i;straw;pedestal",i),straws[i],HALF,straws[i]+HALF,(Int_t)PMAX/2,0,PMAX);
  }
  cdc_ped_ring[11]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow11));
  cdc_ped_ring[23]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow23));
  xd->cd();


  gDirectory->mkdir("rings_initpedestal","CDC rings: initial pedestal from raw window data vs straw")->cd();
  for (i=1; i<29; i++) {
    cdc_initped_ring[i]   = new TH2I(Form("cdc_initped_ring[%i]",i),Form("CDC initial pedestal (ADC units) from raw window data, ring %i;straw;pedestal",i),straws[i],HALF,straws[i]+HALF,(Int_t)PMAX/2,0,PMAX);
  }
  cdc_initped_ring[11]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow11));
  cdc_initped_ring[23]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow23));
  xd->cd();



  gDirectory->mkdir("rings_raw_t","CDC rings: raw time vs straw")->cd();
  for (i=1; i<29; i++) {
    cdc_raw_t_ring[i] = new TH2I(Form("cdc_raw_t_ring[%i]",i),Form("CDC raw time (0.8ns), ring %i;straw;raw time (0.8ns)",i),straws[i],HALF,straws[i]+HALF,256,0,RTVSNMAX);
  }
  cdc_raw_t_ring[11]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow11));
  cdc_raw_t_ring[23]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow23));
  xd->cd();
    

  gDirectory->mkdir("rings_amp","CDC rings: amplitude")->cd();
  for (i=1; i<29; i++) {
    cdc_amp_ring[i]   = new TH2I(Form("cdc_amp_ring[%i]",i),Form("CDC amplitude (ADC units), ring %i",i),straws[i],HALF,straws[i]+HALF,256,0,AMAX);
  }
  cdc_amp_ring[11]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow11));
  cdc_amp_ring[23]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow23));
  xd->cd();


  // gDirectory->mkdir("rings_integral","CDC rings: integral vs straw (pedestal subtracted)")->cd();
  // for (i=1; i<29; i++) {
  //   cdc_int_ring[i]   = new TH2I(Form("cdc_int_ring[%i]",i),Form("CDC integral (ADC units), pedestal subtracted, ring %i",i),straws[i],HALF,straws[i]+HALF,128,0,IMAX);
  // }
  // xd->cd();


  gDirectory->mkdir("rings_integral_incl_ped","CDC rings: integral vs straw (including pedestal)")->cd();
  for (i=1; i<29; i++) {
    cdc_intpp_ring[i]   = new TH2I(Form("cdc_intpp_ring[%i]",i),Form("CDC integral (ADC units), including pedestal, ring %i",i),straws[i],HALF,straws[i]+HALF,100,0,IPPMAX);
  }
  cdc_intpp_ring[11]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow11));
  cdc_intpp_ring[23]->GetXaxis()->SetTitle(Form("Straw number, %s",deadrow23));


  main->cd();    // back to main dir

  return NOERROR;


}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert_2::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes

  return NOERROR;

}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert_2::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  float q,t;         // dcdchits quantities charge, time

  uint32_t rocid;
  uint32_t slot;
  uint32_t channel;

  uint16_t ring,straw; // ring and straw numbers from either dcdchits or dcdcdigihits
  uint16_t n;         // straw number, 1 to 3522

  uint32_t qf,ocount; // time quality factor and overflow count from new firmware
  uint32_t rt,p,a; // dcdcdigihits raw quantities: time, pedestal, amplitude, quality factor, overflow count
  uint32_t integral; // dcdcdigihits integral, includes pedestal
  uint32_t integ;    // dcdcdigihits integral minus pedestal

  uint16_t originalq; //last digit of le_time if qf=1

  uint32_t total_ped;  //total pedestal during integration period
  uint32_t initped; //pedestal calculated from WRD at start of window

  // default scaling factors will be overridden by Df125Config if present
  uint16_t ISCALE = 16;  //scaling factor for integral
  uint16_t ASCALE = 8;   //amplitude
  uint16_t PSCALE = 1;   //ped
  uint16_t NW = 200;
  uint16_t IE = 200;


  const uint16_t NPEDSAMPLES=16; //number of samples to use for initial pedestal initped calculated from window raw data if present

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


	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK


  for (uint32_t i=0; i<hits.size(); i++) {

    const DCDCHit *hit = hits[i]; 
    
    if (hit->q>0.0) {

      q      = hit->q;     // in fC
      t      = hit->t;      // in nanoseconds
      ring   = hit->ring;
      straw  = hit->straw;
      
      n = straw_offset[ring] + straw;

      if (q > 0.0) {
        cdc_e->Fill(q);
        cdc_e_vs_n->Fill(n,q);    
      }
  
      cdc_t->Fill(t);
      cdc_t_vs_n->Fill(n,t);    

      cdc_e_vs_t->Fill(t,q);
      cdc_e_vs_t_ring[ring]->Fill(t,q);

      //cdc_e_ring[ring]->Fill(straw,q);
      cdc_t_ring[ring]->Fill(straw,t);
    }
  }



  const DCDCDigiHit *digihit = NULL;
  const Df125CDCPulse *cp = NULL;
  const Df125WindowRawData *wrd = NULL;
  const Df125Config *cf = NULL;

  for (uint32_t i=0; i<digihits.size(); i++) {

    digihit = digihits[i];

    ring     = digihit->ring;
    straw    = digihit->straw;

    n = straw_offset[ring] + straw;
  

    total_ped = 0;
    originalq = 0;


    //new firmware uses Df125CDCPulseData

    cp = NULL; 
    digihit->GetSingle(cp);

    if (!cp) continue; //no CDCPulseData (happens occasionally)

    cp->GetSingle(cf);
    if (cf) {
      ISCALE = 1<<cf->IBIT;
      ASCALE = 1<<cf->ABIT;
      PSCALE = 1<<cf->PBIT;

      NW = cf->NW;
      IE = cf->IE;

    }

    rocid = cp->rocid;
    slot = cp->slot;
    channel = cp->channel;

    rt = cp->le_time;
    qf = cp->time_quality_bit;
    ocount = cp->overflow_count;

    a = ASCALE*cp->first_max_amp;
    p = PSCALE*cp->pedestal;
    integral = ISCALE*cp->integral;


    int lastsample = NW-20-1;  //eg sample 179 is the last sample integrated for NW=200
    int timesample = int(0.1*rt);
    if (timesample+IE < lastsample) lastsample = timesample+IE;

    int pulselength = 1 + lastsample - timesample;

    integ = integral - p*pulselength;

    if (qf==0) {

      originalq = 0;

      cdc_rt_qf0->Fill(rt);
    
    } else {

      originalq = rt - 10*int(0.1*rt);
   
      cdc_qf_vs_n->Fill(n,originalq);
      cdc_qf_vs_a->Fill(a,originalq);
      cdc_qf_vs_rt->Fill(rt,originalq);

      cdc_o_badt->Fill(straw,ring);
      cdc_ped_badt->Fill(p);
      cdc_raw_t_badt->Fill(rt);
      cdc_amp_badt->Fill(a); 

      //cdc_int_badt->Fill(integ); 
      //cdc_intpp_badt->Fill(integral); 

    }

    cdc_qf->Fill(originalq);

    cdc_rt->Fill(rt);
    cdc_rt_vs_n->Fill(n,rt);

    cdc_amp->Fill(a);
    cdc_amp_vs_n->Fill(n,a);
                 
    cdc_int_vs_raw_t->Fill(rt,integ);
    cdc_int_vs_raw_t_ring[ring]->Fill(rt,integ);
      
    cdc_ped_ring[ring]->Fill(straw,p);
    cdc_raw_t_ring[ring]->Fill(straw,rt);
    cdc_amp_ring[ring]->Fill(straw,a);   //no ped subtraction 
    //cdc_int_ring[ring]->Fill(straw,integ);
    cdc_intpp_ring[ring]->Fill(straw,integral);

 
      
    if (ocount>0) {  // overflow samples present
      cdc_o_overflow->Fill(straw,ring);
      cdc_ped_overflow->Fill(p);
      cdc_raw_t_overflow->Fill(rt); 
    } 


    if (rocid == 25) cdc_ped_roc25->Fill(100*slot + channel,p);
    if (rocid == 26) cdc_ped_roc26->Fill(100*slot + channel,p);
    if (rocid == 27) cdc_ped_roc27->Fill(100*slot + channel,p);
    if (rocid == 28) cdc_ped_roc28->Fill(100*slot + channel,p);


    // initial pedestals from window raw data samples if available

    wrd = NULL;
    cp->GetSingle(wrd);
    if (!wrd) continue;


    initped=0;

    if (wrd->samples.size()>=NPEDSAMPLES) {

      for (uint16_t j=0; j<NPEDSAMPLES; j++) initped += (uint32_t)wrd->samples[j]; 
          
      initped = (uint32_t)initped/16.0;

      if (initped > 0) {

        if (rocid == 25) cdc_initped_roc25->Fill(100*slot + channel,initped);
        if (rocid == 26) cdc_initped_roc26->Fill(100*slot + channel,initped);
        if (rocid == 27) cdc_initped_roc27->Fill(100*slot + channel,initped);
        if (rocid == 28) cdc_initped_roc28->Fill(100*slot + channel,initped);

      }

    }

    if (initped) cdc_initped_ring[ring]->Fill(straw,initped);

  } //end for each digihit


	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK


  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert_2::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_expert_2::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
