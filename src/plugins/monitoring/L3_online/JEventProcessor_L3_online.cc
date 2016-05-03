// $Id$
//
//    File: JEventProcessor_L3_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>
#include <random>


#include "JEventProcessor_L3_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;


// for root
#include <TDirectory.h>
#include <TH1.h>


// root hist pointers
static TH1I * l3;



//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_L3_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_L3_online::JEventProcessor_L3_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_L3_online::~JEventProcessor_L3_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_L3_online::init(void) {

  // create root folder for evnt and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("l3")->cd();


  // book hist
  l3 = new TH1I("l3","l3",100,0,100);


  // back to main dir
  main->cd();

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_L3_online::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_L3_online::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {

  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.


  japp->RootWriteLock();

  // fill hist

  japp->RootUnLock(); 

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_L3_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_L3_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
