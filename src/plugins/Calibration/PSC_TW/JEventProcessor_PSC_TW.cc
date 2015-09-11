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
#include <TTree.h>
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

// Define TTree
static TTree *PSC_tree;

// Define branches
Double_t tdc_l;		// PSC TDC time for the left modules
Double_t tdc_r;		// PSC TDC time for the right modules
Double_t rf_l;		// RF time associated with the left PSC TDC time
Double_t rf_r;		// RF time associated with the right PSC TDC time
Double_t pp_l;		// PSC pulse peak for the left modules
Double_t pp_r;		// PSC pulse peak for the right modules
Int_t psc_mod_l;	// PSC left module number
Int_t psc_mod_r;	// PSC right modules number

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
// init
//------------------
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

   // Create the tree
   PSC_tree = new TTree("PSC_tree","PSC_tree");
   // Create branches
   PSC_tree->Branch("tdc_l", &tdc_l, "tdc_l/D");
   PSC_tree->Branch("tdc_r", &tdc_r, "tdc_r/D");
   PSC_tree->Branch("rf_l", &rf_l, "rf_l/D");
   PSC_tree->Branch("rf_r", &rf_r, "rf_r/D");
   PSC_tree->Branch("pp_l", &pp_l, "pp_l/D");
   PSC_tree->Branch("pp_r", &pp_r, "pp_r/D");
   PSC_tree->Branch("psc_mod_l", &psc_mod_l, "psc_mod_l/I");
   PSC_tree->Branch("psc_mod_r", &psc_mod_r, "psc_mod_r/I");

   japp->RootUnLock();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_PSC_TW::brun(JEventLoop *eventLoop, int runnumber)
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
jerror_t JEventProcessor_PSC_TW::evnt(JEventLoop *loop, int eventnumber)
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
   loop->Get(locRFTimes);

   japp->RootWriteLock();

   if (pairs.size() > 0) {
      psc_mod_l = pairs[0]->ee.first->module;
      pp_l = pairs[0]->ee.first->pulse_peak;
      tdc_l = pairs[0]->ee.first->t;
      rf_l = dRFTimeFactory->Step_TimeToNearInputTime(locRFTimes[0]->dTime, tdc_l);

      psc_mod_r = pairs[0]->ee.second->module;
      pp_r = pairs[0]->ee.second->pulse_peak;
      tdc_r = pairs[0]->ee.second->t;
      rf_r = dRFTimeFactory->Step_TimeToNearInputTime(locRFTimes[0]->dTime, tdc_r);

      PSC_tree->Fill();
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
