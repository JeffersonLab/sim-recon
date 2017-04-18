// $Id$
//
//    File: JEventProcessor_CDC_drift.cc
// Created: Wed Oct 22 2014
// Creator: Naomi Jarvis


#include <stdint.h>
#include <vector>

#include <TMath.h>


#include "JEventProcessor_CDC_drift.h"
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>


using namespace std;
using namespace jana;


#include "CDC/DCDCHit.h"
#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125CDCPulse.h"

#include "TRIGGER/DTrigger.h"

#include <TDirectory.h>
#include <TH2.h>
#include <TH1.h>
#include <TF1.h>
#include <TTree.h>
#include <TBranch.h>

// root hist pointers

static TH1F *cdc_time = NULL;
static TH1I *cdc_rawtime = NULL;

static TTree *tfit = NULL;
static TTree *rtfit = NULL;

static bool FIT_TIME = true;


//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_CDC_drift());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_drift::JEventProcessor_CDC_drift() {
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_drift::~JEventProcessor_CDC_drift() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_CDC_drift::init(void) {

  /*
  // max values for histogram scales, modified fa250-format readout
  const Int_t RTMAX = 12000; //max for raw time, less than full field width
  const Char_t rtunits[8] = "0.125ns";  //raw time is in units of sample/64 = ns/8
  */


  // raw quantities for read out (125 format) are
  //   time                    field max 2047   scaled x 1, units 0.8ns
  //   time qf                 field max 1 
  //   overflow count          field max 7
  //   pedestal                field max 255    scaled x 1/4 initially
  //   max amplitude 9 bits,   field max 511    scaled x 1/8
  //   integral                field max 16383  scaled x 1/14


  // max values for histogram scales, fa125-format readout


  const Char_t rtunits[6] = "0.8ns";  //raw time is in units of sample/10 = 0.8n

  const Int_t RTMAX = 2048;  //max for raw time histo, fa125-format, 11 bits
  const Int_t TMAX = 2000;  //max for time


  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  // create root folder for cdc and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("CDC_drift")->cd();


  // book histograms

  cdc_time = new TH1F("cdc_time","CDC time (units of ns); time (ns)",TMAX,-500,TMAX-500);
  cdc_rawtime = new TH1I("cdc_rawtime",Form("CDC raw time (units of %s); raw time (%s)",rtunits,rtunits),RTMAX,0,RTMAX);


  if (FIT_TIME) {

    rtfit = new TTree("rawtimefit","raw drift time fit params");

    Long64_t tentries;
    rtfit->Branch("entries",&tentries,"entries/L");

    Double_t t0;
    rtfit->Branch("t0",&t0,"t0/D");

    Double_t tmax;
    rtfit->Branch("tmax",&tmax,"tmax/D");

    Double_t tmax_slope;
    rtfit->Branch("tmax_slope",&tmax_slope,"tmax_slope/D");

    Double_t tdiff;
    rtfit->Branch("tdiff_ns",&tdiff,"tdiff/D");


    tfit = new TTree("timefit","drift time fit params");

    tfit->Branch("entries",&tentries,"entries/L");

    tfit->Branch("t0",&t0,"t0/D");

    tfit->Branch("tmax",&tmax,"tmax/D");

    tfit->Branch("tmax_slope",&tmax_slope,"tmax_slope/D");

    tfit->Branch("tdiff_ns",&tdiff,"tdiff/D");

  } 

  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  main->cd();

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_drift::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_drift::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  // cosmics, estimate 15 mins ~ 4.4e5 events ~ 4.4e5*82/372 ~ 1e5 useful hits



  const uint32_t MIN_EVENTS = 50000;            //min events to collect before fitting drift time
  const uint32_t UPDATE_INTERVAL = 25000;   //incremental events required to update fit

  const Bool_t RESET = kFALSE;             // if true, zero histos after fitting
  const Bool_t VERBOSE = kFALSE;           // if true, print fits to stdout

  uint16_t ring,straw; // ring and straw numbers from either dcdchits or dcdcdigihits
  uint16_t n;         // straw number, 1 to 3522
  uint16_t j;

  Long64_t nentries;  // current number of entries

  int64_t previous,tprevious;  // number of entries when histo was last fitted


  //array to make straw number n; add extra 0 at front to use offset[1] for ring 1
  int straw_offset[29] = {0,0,42,84,138,192,258,324,404,484,577,670,776,882,1005,1128,1263,1398,1544,1690,1848,2006,2176,2346,2528,2710,2907,3104,3313};

  const uint16_t nstraws = 77;  //size of strawlist - list of n of straws to include in fit

  const uint16_t strawlist[] = {176, 237, 496, 497, 775, 776, 777, 782, 879, 881, 882, 895, 900, 1021, 1026, 1047, 1052, 1056, 1057, 1130, 1241, 1252, 1266, 1318, 1340, 1376, 1567, 1568, 1679, 1682, 1701, 1849, 1853, 1864, 1918, 1998, 2088, 2242, 2244, 2248, 2255, 2256, 2430, 2445, 2556, 2585, 2748, 2767, 2770, 2772, 2774, 2782, 2788, 2789, 2793, 2796, 2943, 2951, 2952, 2962, 2963, 2965, 2969, 2973, 2985, 3159, 3160, 3176, 3177, 3184, 3214, 3361, 3363, 3365, 3369, 3428, 3429};


  Bool_t fillhisto;    // fill histo if true
  Bool_t fithisto;     // fit histo if true

  const DTrigger* locTrigger = NULL; 
  eventLoop->GetSingle(locTrigger); 
  if(locTrigger->Get_L1FrontPanelTriggerBits() != 0)
    return NOERROR;
  if (!locTrigger->Get_IsPhysicsEvent()){ // do not look at PS triggers
    return NOERROR;
  }
	

  // get raw data for cdc
  vector<const DCDCDigiHit*> digihits;
  eventLoop->Get(digihits);


  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!

  fithisto = kFALSE;

  previous =  (uint32_t)(cdc_rawtime->GetEntries()/UPDATE_INTERVAL);


  for (uint32_t i=0; i<digihits.size(); i++) {

    const DCDCDigiHit *digihit = digihits[i];  
    const Df125CDCPulse *cp = NULL;

    digihit->GetSingle(cp);

    if (!cp) continue;

    ring   = digihit->ring;
    straw  = digihit->straw;
    n      = straw_offset[ring] + straw;

    fillhisto = kFALSE;

    if ((digihit->pulse_time) && (!cp->time_quality_bit)) {

      j=0;

      while ((!fillhisto) && (j<nstraws)) {

        if (n == strawlist[j]) fillhisto = kTRUE;
        j++;

      }


      if (fillhisto) { 

        cdc_rawtime->Fill(digihit->pulse_time);
        nentries = cdc_rawtime->GetEntries();
        if ((nentries > MIN_EVENTS) && (uint32_t(nentries/UPDATE_INTERVAL) > previous)) fithisto = kTRUE;
        if (!FIT_TIME) fithisto = kFALSE;

      }

    } 

  }

  Bool_t fitthisto;     // fit histo if true

  // get raw data for cdc
  vector<const DCDCHit*> hits;
  eventLoop->Get(hits);

  fitthisto = kFALSE;

  tprevious =  (uint32_t)(cdc_time->GetEntries()/UPDATE_INTERVAL);


  for (uint32_t i=0; i<hits.size(); i++) {

    const DCDCHit *hit = hits[i];  

    ring   = hit->ring;
    straw  = hit->straw;
    n      = straw_offset[ring] + straw;

    fillhisto = kFALSE;

    j=0;

    while ((!fillhisto) && (j<nstraws)) {

      if (n == strawlist[j]) fillhisto = kTRUE;
      j++;

    }


    if (fillhisto) { 

      cdc_time->Fill(hit->t);
      nentries = cdc_time->GetEntries();
      if ((nentries > MIN_EVENTS) && (uint32_t(nentries/UPDATE_INTERVAL) > tprevious)) fitthisto = kTRUE;
      if (!FIT_TIME) fitthisto = kFALSE;

    }

  }



  if (fithisto && FIT_TIME) {

    //***  the quantities to save are fitstatus, fitparams 0 to 9 and tdiff  ***

    if (VERBOSE) printf("\n\nFitting cdc_rawtime\n");

    const float TUNITS = 0.8;  // fa125 - time is in units of 0.8ns
    //  const float TUNITS = 0.125;  // fa250 - time is in units of 0.125ns

    Double_t fitparams[10]; 
    Float_t startpar[10];
    Int_t fitstatus;
    Double_t tdiff=0;  // max drift time in ns

    Int_t imin = cdc_rawtime->FindFirstBinAbove(); // first histogram bin with counts
    Int_t imax;      // last bin with counts 
    Int_t ipeak=0;  // bin with t0 peak in, used to find startparam

    Double_t xmin;   // x value of imin - used for fit range
    Double_t xmax;  // x value of imax
    
    Int_t bgpeakwidth = 10; // width in bins of initial background peak, if there is one
    Int_t bgrange = 40; //scan 4 full samples - include this many bins in background estimate

    Int_t i;

    Int_t chunk1 = 0; 
    for (i=0; i<bgpeakwidth; i++) chunk1 += cdc_rawtime->GetBinContent(imin+i);

    Int_t chunk2 = 0; 
    for (i=0; i<bgpeakwidth; i++) chunk2 += cdc_rawtime->GetBinContent(imin+bgpeakwidth+i);

    // skip over noise peak
    if (chunk1>1.5*chunk2) imin += bgpeakwidth;
    xmin = cdc_rawtime->GetXaxis()->GetBinLowEdge(imin);

    // find max content bin
    // search on from peak to find counts=0, end of fit range

    Int_t maxcontent=0;
    for (i=imin; i<cdc_rawtime->GetNbinsX(); i++) {
      if (cdc_rawtime->GetBinContent(i) > maxcontent) {
        maxcontent = cdc_rawtime->GetBinContent(i);
        ipeak = i;
      }
    }

    Double_t xpeak = cdc_rawtime->GetXaxis()->GetBinLowEdge(ipeak);

    i=ipeak; 
    while (cdc_rawtime->GetBinContent(i) > 0 && i<cdc_rawtime->GetNbinsX() ) i++;

    imax = i;
    xmax = imax*cdc_rawtime->GetBinWidth(imax);

    //starting point for background height
    Double_t bg = 0;      
    for (i=imin; i<imin+bgrange; i++) bg += cdc_rawtime->GetBinContent(i);
    bg = bg/(Double_t)bgrange;


    TF1 *f = new TF1("f","[9] + [0] * (1 + [1]*exp(([3]-x)/[2]) + [7]*exp(([3]-x)/[8]) ) / ( (1+exp(([3]-x)/[5])) * (1+exp((x-[4])/[6])) )",xmin,xmax);

    f->SetLineWidth(1);
    f->SetLineColor(6);

    // set start values and limits here for all fit params except 0,3,4

    startpar[1] = 15; //amplitude of first exp contrib to peak
    startpar[7] = 3; //amplitude of second exp contrib to peak

    f->SetParLimits(1,0,startpar[1]*2);  //prev *10
    f->SetParLimits(7,0,startpar[7]*2);  //prev *10

    startpar[5] = 5*0.8/TUNITS; //slope up of t0 edge
    startpar[6] = 25*0.8/TUNITS; //slope down of tmax edge

    f->SetParLimits(5,0,startpar[5]*2.5);   //prev *2
    f->SetParLimits(6,0,startpar[6]*2.5);   //prev *2

    startpar[2] = 20*0.8/TUNITS; //first exp fall-off
    startpar[8] = 200*0.8/TUNITS; //second exp fall-off

    f->SetParLimits(2,0,startpar[2]*3);
    f->SetParLimits(8,startpar[2]*3,startpar[8]*3);

    for (j=1;j<3;j++) f->SetParameter(j,startpar[j]);
    for (j=5;j<9;j++) f->SetParameter(j,startpar[j]);

    previous = (uint32_t)(nentries/UPDATE_INTERVAL);

    // start values & limits for fit params 0,3,4 depend on nentries

    startpar[0] = 0.0005*nentries;  //overall scaling factor
    startpar[9] = bg; //noise background

    f->SetParLimits(0,0,startpar[0]*100);
    f->SetParameter(0,startpar[0]);

    f->SetParLimits(9,0,bg*2);
    f->SetParameter(9,startpar[9]);

    startpar[3] = xpeak;  //t0
    startpar[4] = xmax; //xpeak+500*0.8/TUNITS; //tmax  //prev 550

    f->SetParLimits(3,startpar[3]-(50*0.8/TUNITS),startpar[3]);
    f->SetParLimits(4,startpar[3]+(500*0.8/TUNITS),xmax);  //min 0.5us

    f->SetParameter(3,startpar[3]);
    f->SetParameter(4,startpar[3] + (700*0.8/TUNITS));

    if (!VERBOSE) fitstatus = cdc_rawtime->Fit("f","QRLL");
    if (VERBOSE) fitstatus = cdc_rawtime->Fit("f","RLL");

    if (fitstatus == 0 || fitstatus == 2) {   //fitstatus 0=good, 2=error matrix not posdef, fit params are correlated

      f->GetParameters(fitparams);

      tdiff = (fitparams[4] - fitparams[3])*TUNITS;

      //cdc_rawtime->SetTitle(Form("Estimated max drift time is %3.2f ns",tdiff));

    } else { 

      tdiff = 0;

    }

    if (VERBOSE) printf("fitstatus:%1i nentries:%5lli [0] %2.0f  [1] %2.0f  [2] %3.0f  [3] %4.0f  [4] %4.0f  [5] %4.1f  [6] %3.0f  [7] %3.1f  [8] %4.0f [9] %4.0f [tmax] %3.0f\n",fitstatus,nentries,fitparams[0],fitparams[1],fitparams[2],fitparams[3],fitparams[4],fitparams[5],fitparams[6],fitparams[7],fitparams[8],fitparams[9],tdiff);


    rtfit->SetBranchAddress("entries",&nentries);
    rtfit->SetBranchAddress("t0",&fitparams[3]);
    rtfit->SetBranchAddress("tmax",&fitparams[4]);
    rtfit->SetBranchAddress("tmax_slope",&fitparams[6]);
    rtfit->SetBranchAddress("tdiff_ns",&tdiff);

    rtfit->Fill();

    // **** reset histogram ****
    if (RESET) cdc_rawtime->Reset();

  }


  if (fitthisto && FIT_TIME) {

    //***  the quantities to save are fitstatus, fitparams 0 to 9 and tdiff  ***

    if (VERBOSE) printf("\n\nFitting cdc_time\n");

    const float TUNITS = 1.0;  // time in ns

    Double_t fitparams[10]; 
    Float_t startpar[10];
    Int_t fitstatus;
    Double_t tdiff=0;  // max drift time in ns

    Int_t imin = cdc_time->FindFirstBinAbove(); // first histogram bin with counts
    Int_t imax;      // last bin with counts 
    Int_t ipeak=0;  // bin with t0 peak in, used to find startparam

    Double_t xmin;   // x value of imin - used for fit range
    Double_t xmax;  // x value of imax
    
    Int_t bgpeakwidth = 8; // width in bins of initial background peak, if there is one
    Int_t bgrange = 32; //scan 4 x 8ns samples - include this many bins in background estimate

    Int_t i;

    Double_t chunk1 = 0; 
    for (i=0; i<bgpeakwidth; i++) chunk1 += cdc_time->GetBinContent(imin+i);

    Double_t chunk2 = 0; 
    for (i=0; i<bgpeakwidth; i++) chunk2 += cdc_time->GetBinContent(imin+bgpeakwidth+i);
    if (chunk1>1.5*chunk2) imin += bgpeakwidth;
    xmin = cdc_time->GetXaxis()->GetBinLowEdge(imin);

    // find max content bin
    // search on from peak to find counts=0, end of fit range

    Int_t maxcontent=0;
    for (i=imin; i<cdc_time->GetNbinsX(); i++) {
      if (cdc_time->GetBinContent(i) > maxcontent) {
        maxcontent = cdc_time->GetBinContent(i);
        ipeak = i;
      }
    }

    Double_t xpeak = cdc_time->GetXaxis()->GetBinLowEdge(ipeak);

    i=ipeak; 
    while (cdc_time->GetBinContent(i) > 0 && i<cdc_time->GetNbinsX() ) i++;

    imax = i;
    xmax = cdc_time->GetBinLowEdge(imax);

    //starting point for background height
    Double_t bg = 0; 
    for (i=imin; i<imin+bgrange; i++) bg += cdc_time->GetBinContent(i);
    bg = bg/(Double_t)bgrange;


    TF1 *f = new TF1("f","[9] + [0] * (1 + [1]*exp(([3]-x)/[2]) + [7]*exp(([3]-x)/[8]) ) / ( (1+exp(([3]-x)/[5])) * (1+exp((x-[4])/[6])) )",xmin,xmax);

    f->SetLineWidth(1);
    f->SetLineColor(6);

    // set start values and limits here for all fit params except 0,3,4

    startpar[1] = 15; //amplitude of first exp contrib to peak
    startpar[7] = 3; //amplitude of second exp contrib to peak

    f->SetParLimits(1,0,startpar[1]*2);  //prev *10
    f->SetParLimits(7,0,startpar[7]*2);  //prev *10

    startpar[5] = 5*0.8/TUNITS; //slope up of t0 edge
    startpar[6] = 25*0.8/TUNITS; //slope down of tmax edge

    f->SetParLimits(5,0,startpar[5]*2.5);   //prev *2
    f->SetParLimits(6,0,startpar[6]*2.5);   //prev *2

    startpar[2] = 20*0.8/TUNITS; //first exp fall-off
    startpar[8] = 200*0.8/TUNITS; //second exp fall-off

    f->SetParLimits(2,0,startpar[2]*3);
    f->SetParLimits(8,startpar[2]*3,startpar[8]*3);

    for (j=1;j<3;j++) f->SetParameter(j,startpar[j]);
    for (j=5;j<9;j++) f->SetParameter(j,startpar[j]);

    previous = (uint32_t)(nentries/UPDATE_INTERVAL);

    // start values & limits for fit params 0,3,4 depend on nentries

    startpar[0] = 0.0005*nentries;  //overall scaling factor
    startpar[9] = bg; //noise background

    f->SetParLimits(0,0,startpar[0]*100);
    f->SetParameter(0,startpar[0]);

    f->SetParLimits(9,0,bg*2);
    f->SetParameter(9,startpar[9]);

    startpar[3] = xpeak;  //t0
    startpar[4] = xmax; //xpeak + 500*0.8/TUNITS; //tmax  //prev 550

    f->SetParLimits(3,startpar[3]-(50*0.8/TUNITS),startpar[3]);
    f->SetParLimits(4,startpar[3]+(500*0.8/TUNITS),xmax);  // min 0.5us

    f->SetParameter(3,startpar[3]);
    f->SetParameter(4,startpar[3] + (700*0.8/TUNITS));

    if (!VERBOSE) fitstatus = cdc_time->Fit("f","QRLL");
    if (VERBOSE) fitstatus = cdc_time->Fit("f","RLL");

    if (fitstatus == 0 || fitstatus == 2) {   //fitstatus 0=good, 2=error matrix not posdef, fit params are correlated

      f->GetParameters(fitparams);

      tdiff = (fitparams[4] - fitparams[3])*TUNITS;

      //cdc_time->SetTitle(Form("Estimated max drift time is %3.2f ns",tdiff));

    } else { 

      tdiff = 0;

    }

    if (VERBOSE) printf("fitstatus:%1i nentries:%5lli [0] %2.0f  [1] %2.0f  [2] %3.0f  [3] %4.0f  [4] %4.0f  [5] %4.1f  [6] %3.0f  [7] %3.1f  [8] %4.0f [9] %4.0f [tmax] %3.0f\n",fitstatus,nentries,fitparams[0],fitparams[1],fitparams[2],fitparams[3],fitparams[4],fitparams[5],fitparams[6],fitparams[7],fitparams[8],fitparams[9],tdiff);


    tfit->SetBranchAddress("entries",&nentries);
    tfit->SetBranchAddress("t0",&fitparams[3]);
    tfit->SetBranchAddress("tmax",&fitparams[4]);
    tfit->SetBranchAddress("tmax_slope",&fitparams[6]);
    tfit->SetBranchAddress("tdiff_ns",&tdiff);

    tfit->Fill();

    // **** reset histogram ****
    if (RESET) cdc_time->Reset();

  }





  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_drift::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_drift::fini(void) {
  // Called before program exit after event processing is finished.


  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
