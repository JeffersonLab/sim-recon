// $Id$
//
//    File: JEventProcessor_PSC_TW.cc
// Created: Fri Aug 21 10:42:28 EDT 2015
// Creator: aebarnes (on Linux ifarm1102 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_PSC_TW.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
// ROOT header fiels
#include <TH2.h>
#include <TDirectory.h>
// RF header files
#include <RF/DRFTime_factory.h>
// PSC header files
#include <PAIR_SPECTROMETER/DPSCPair.h>

// Define constants
const uint32_t NMODULES = 8;

const double TMIN = -5;
const double TMAX = 5;
const double TBINS = (TMAX - TMIN)/0.1;

const double PMIN = 0;
const double PMAX = 2000;
const double PBINS = (PMAX - PMIN)/2;

// Define histograms
static TH2F* h_dt_vs_pp_tdc_l[NMODULES];
static TH2F* h_dt_vs_pp_tdc_r[NMODULES];
static TH2F* h_dt_vs_pp_t_l[NMODULES];
static TH2F* h_dt_vs_pp_t_r[NMODULES];

// Define variables
Int_t psc_mod_l;
Int_t psc_mod_r;
Double_t pp_l;
Double_t pp_r;
Double_t adc_l;
Double_t adc_r;
Double_t tdc_l;
Double_t tdc_r;
Double_t t_l;
Double_t t_r;
Double_t rf_l;
Double_t rf_r;

// Define RFTime_factory
DRFTime_factory* dRFTimeFactory;

extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_PSC_TW());
}
} // "C"


//------------------
// JEventProcessor_PSC_TW (Constructor)
//------------------
JEventProcessor_PSC_TW::JEventProcessor_PSC_TW()
{

}

//------------------
// ~JEventProcessor_PSC_TW (Destructor)
//------------------
JEventProcessor_PSC_TW::~JEventProcessor_PSC_TW()
{

}

//------------------
// init //------------------
jerror_t JEventProcessor_PSC_TW::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

   TDirectory *pscDir = gDirectory->mkdir("PSC_TW");
   pscDir->cd();

   gDirectory->mkdir("tdc-rf")->cd();
   for (uint32_t i = 0; i < NMODULES; ++i)
   {
      h_dt_vs_pp_tdc_l[i] = new TH2F(Form("h_dt_vs_pp_tdc_l_%i",i+1),
                                     Form("#Delta t vs. pulse peak, raw TDC, left PSC %i;\
                                           Pulse peak; #Delta t (raw TDC - RF)",i+1),
                                     PBINS, PMIN, PMAX, TBINS, TMIN, TMAX);
      h_dt_vs_pp_tdc_r[i] = new TH2F(Form("h_dt_vs_pp_tdc_r_%i",i+1),
                                     Form("#Delta t vs. pulse peak, raw TDC, left PSC %i;\
                                           Pulse peak; #Delta t (raw TDC - RF)",i+1),
                                     PBINS, PMIN, PMAX, TBINS, TMIN, TMAX);
   }
   pscDir->cd();

   gDirectory->mkdir("t-rf")->cd();
   for (uint32_t i = 0; i < NMODULES; ++i)
   {
      h_dt_vs_pp_t_l[i] = new TH2F(Form("h_dt_vs_pp_t_l_%i",i+1),
                                   Form("#Delta t vs. pulse peak, corrected TDC, left PSC %i;\
                                         Pulse Peak; #Delta t (raw TDC - RF)",i+1),
                                   PBINS, PMIN, PMAX, TBINS, TMIN, TMAX);
      h_dt_vs_pp_t_r[i] = new TH2F(Form("h_dt_vs_pp_t_r_%i",i+1),
                                   Form("#Delta t vs. pulse peak, corrected TDC, left PSC %i;\
                                         Pulse peak; #Delta t (raw TDC - RF)",i+1),
                                   PBINS, PMIN, PMAX, TBINS, TMIN, TMAX);
   }
   pscDir->cd();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_PSC_TW::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes

   //////////////
   // RF
   //////////////

   // Initialize RF time factory
   dRFTimeFactory = static_cast<DRFTime_factory*>(eventLoop->GetFactory("DRFTime"));

   // be sure that DRFTime_factory::init() and brun() are called
   vector<const DRFTime*> locRFTimes;

   eventLoop->Get(locRFTimes);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_PSC_TW::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();

   vector<const DRFTime*>	locRFTimes;
   vector<const DPSCPair*>	pairs;

   loop->Get(pairs);
   loop->Get(locRFTimes,"PSC");
   const DRFTime* locRFTime = NULL;

   if (locRFTimes.size() > 0)
      locRFTime = locRFTimes[0];
   else
      return NOERROR;

   // Since we are filling histograms local to this plugin, 
   // it will not interfere with other ROOT operations:
   // can use plugin-wide ROOT fill lock
   japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

   for (uint32_t i = 0; i < pairs.size(); ++i)
   {
      // Left modules
      psc_mod_l = pairs[i]->ee.first->module;
      pp_l = pairs[i]->ee.first->pulse_peak;
      adc_l = pairs[i]->ee.first->time_fadc;
      tdc_l = pairs[i]->ee.first->time_tdc;
      t_l = pairs[i]->ee.first->t;
      // Right modules
      psc_mod_r = pairs[i]->ee.second->module;
      pp_r = pairs[i]->ee.second->pulse_peak;
      adc_r = pairs[i]->ee.second->time_fadc;
      tdc_r = pairs[i]->ee.second->time_tdc;
      t_r = pairs[i]->ee.second->t;

      // Use the ADC time to find the closest RF time.
      rf_l = dRFTimeFactory->Step_TimeToNearInputTime(locRFTime->dTime,adc_l);
      rf_r = dRFTimeFactory->Step_TimeToNearInputTime(locRFTime->dTime,adc_r);

      h_dt_vs_pp_tdc_l[psc_mod_l - 1]->Fill(pp_l, tdc_l - rf_l);
      h_dt_vs_pp_t_l[psc_mod_l - 1]->Fill(pp_l, t_l - rf_l);
      h_dt_vs_pp_tdc_r[psc_mod_r - 1]->Fill(pp_r, tdc_r - rf_r);
      h_dt_vs_pp_t_r[psc_mod_r - 1]->Fill(pp_r, t_r - rf_r);
   }

   japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_PSC_TW::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_PSC_TW::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}
