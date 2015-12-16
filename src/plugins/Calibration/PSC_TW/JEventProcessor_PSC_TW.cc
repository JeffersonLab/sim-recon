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
// C++ header files
#include <sstream>
#include <stdio.h>
#include <iostream>
#include <cmath>
#include <limits>
#include <vector>
// RF header files
#include <RF/DRFTime_factory.h>
// PSC header files
#include <PAIR_SPECTROMETER/DPSCPair.h>

// Define constants
const uint32_t NMODULES = 8;

// Define histograms
static TH2F* h_dt_vs_pp_l[NMODULES];
static TH2F* h_dt_vs_pp_r[NMODULES];

// Define variables
Int_t psc_mod_l;
Int_t psc_mod_r;
Double_t pp_l;
Double_t pp_r;
Double_t tdc_l;
Double_t tdc_r;
//Double_t rf_l;
double rf_l;
//Double_t rf_r;
double rf_r;

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

   japp->RootWriteLock();

   // Name histograms
   for (uint32_t i = 0; i < NMODULES; ++i) {
      h_dt_vs_pp_l[i] = new TH2F(Form("h_dt_vs_pp_l_%i",i+1),Form("Time difference vs. pulse peak for left PSC module %i",i+1),1000,0,1000,100,-5,5);
      h_dt_vs_pp_r[i] = new TH2F(Form("h_dt_vs_pp_r_%i",i+1),Form("Time difference vs. pulse peak for right PSC module %i",i+1),1000,0,1000,100,-5,5);
   }

   japp->RootUnLock();

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

   japp->RootWriteLock();

   if (pairs.size() > 0) {
      psc_mod_l = pairs[0]->ee.first->module;
      pp_l = pairs[0]->ee.first->pulse_peak;
      tdc_l = pairs[0]->ee.first->t;
      rf_l = dRFTimeFactory->Step_TimeToNearInputTime(locRFTime->dTime, tdc_l);
      h_dt_vs_pp_l[psc_mod_l - 1]->Fill(pp_l,tdc_l - rf_l);

      psc_mod_r = pairs[0]->ee.second->module;
      pp_r = pairs[0]->ee.second->pulse_peak;
      tdc_r = pairs[0]->ee.second->t;
      rf_r = dRFTimeFactory->Step_TimeToNearInputTime(locRFTime->dTime, tdc_r);
      h_dt_vs_pp_r[psc_mod_r - 1]->Fill(pp_r,tdc_r - rf_r);
   }

   japp->RootUnLock();

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
