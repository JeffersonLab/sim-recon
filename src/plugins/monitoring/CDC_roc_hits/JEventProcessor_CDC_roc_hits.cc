// $Id$
//
//    File: JEventProcessor_CDC_roc_hits.cc
// Created: 18 May 2015
// Creator: Naomi Jarvis
//
// Plot channel number vs slot for any hit in each roc
// Useful to see which channels are very noisy, or dead
// Same purpose as occupancy, but organized by roc & slot for easier diagnostics

#include <stdint.h>
#include <vector>

#include "JEventProcessor_CDC_roc_hits.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "CDC/DCDCDigiHit.h"
#include "DAQ/Df125PulseIntegral.h"
#include "DAQ/Df125PulsePedestal.h"

#include <TDirectory.h>
#include <TH2.h>

// root hist pointers

static TH1I *cdc_nevents;

static TH2D *cdc_hits_roc25;  
static TH2D *cdc_hits_roc26;  
static TH2D *cdc_hits_roc27;  
static TH2D *cdc_hits_roc28;  

static TH2D *cdc_amp_roc25;  
static TH2D *cdc_amp_roc26;  
static TH2D *cdc_amp_roc27;  
static TH2D *cdc_amp_roc28;  


//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_CDC_roc_hits());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_roc_hits::JEventProcessor_CDC_roc_hits() {
}


//----------------------------------------------------------------------------------


JEventProcessor_CDC_roc_hits::~JEventProcessor_CDC_roc_hits() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_CDC_roc_hits::init(void) {

  // create root folder for cdc and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("CDC_roc_hits")->cd();

  // book histograms
  cdc_nevents = new TH1I("cdc_nevents","CDC number of events",1, 0.5, 1.5);
  cdc_hits_roc25   = new TH2D("cdc_hits_roc25","CDC hits in ROC 25, channel vs slot;slot;channel",15,2.5,17.5,72,0.5,72.5);
  cdc_hits_roc26   = new TH2D("cdc_hits_roc26","CDC hits in ROC 26, channel vs slot;slot;channel",14,2.5,16.5,72,0.5,72.5);
  cdc_hits_roc27   = new TH2D("cdc_hits_roc27","CDC hits in ROC 27, channel vs slot;slot;channel",14,2.5,16.5,72,0.5,72.5);
  cdc_hits_roc28   = new TH2D("cdc_hits_roc28","CDC hits in ROC 28, channel vs slot;slot;channel",15,2.5,17.5,72,0.5,72.5);

  cdc_amp_roc25   = new TH2D("cdc_amp_roc25","CDC pulse peak amplitude in ROC 25;slot*100+channel;pulse peak (ADC units)",2000,0,2000,4096,0,4096);
  cdc_amp_roc26   = new TH2D("cdc_amp_roc26","CDC pulse peak amplitude in ROC 26;slot*100+channel;pulse peak (ADC units)",2000,0,2000,4096,0,4096);
  cdc_amp_roc27   = new TH2D("cdc_amp_roc27","CDC pulse peak amplitude in ROC 27;slot*100+channel;pulse peak (ADC units)",2000,0,2000,4096,0,4096);
  cdc_amp_roc28   = new TH2D("cdc_amp_roc28","CDC pulse peak amplitude in ROC 28;slot*100+channel;pulse peak (ADC units)",2000,0,2000,4096,0,4096);
  
  main->cd();

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.


 
  // get raw data for cdc
  vector<const DCDCDigiHit*> digihits;
  eventLoop->Get(digihits);

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

  if(digihits.size() > 0) cdc_nevents->Fill(1);

  for(uint32_t i=0; i<digihits.size(); i++) {

    const DCDCDigiHit *digihit = digihits[i];  // avoids having to use the uglier “cdcdigihits[0]->” syntax
    const Df125PulseIntegral *pi = NULL;       // Get pointer to the underlying object of interest

    digihit->GetSingle(pi);

    if (pi->rocid == 25) cdc_hits_roc25->Fill(pi->slot-0.5,pi->channel);
    if (pi->rocid == 26) cdc_hits_roc26->Fill(pi->slot,pi->channel);
    if (pi->rocid == 27) cdc_hits_roc27->Fill(pi->slot,pi->channel);
    if (pi->rocid == 28) cdc_hits_roc28->Fill(pi->slot,pi->channel);

    const Df125PulsePedestal *pp = NULL;       // Get pointer to the underlying object of interest

    digihit->GetSingle(pp);

    if (!pp) continue;

    if (pp->pulse_peak) {
      if (pp->rocid == 25) cdc_amp_roc25->Fill(pp->slot*100+pp->channel,pp->pulse_peak);
      if (pp->rocid == 26) cdc_amp_roc26->Fill(pp->slot*100+pp->channel,pp->pulse_peak);
      if (pp->rocid == 27) cdc_amp_roc27->Fill(pp->slot*100+pp->channel,pp->pulse_peak);
      if (pp->rocid == 28) cdc_amp_roc28->Fill(pp->slot*100+pp->channel,pp->pulse_peak);
    }

  }

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_CDC_roc_hits::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
