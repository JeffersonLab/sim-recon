// $Id$
//
//    File: JEventProcessor_bigevents_skim.cc
// Created: Thu May 12 08:01:59 EDT 2016
// Creator: zihlmann (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_bigevents_skim.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_bigevents_skim());
  }
} // "C"


//------------------
// JEventProcessor_bigevents_skim (Constructor)
//------------------
JEventProcessor_bigevents_skim::JEventProcessor_bigevents_skim()
{

}

//------------------
// ~JEventProcessor_bigevents_skim (Destructor)
//------------------
JEventProcessor_bigevents_skim::~JEventProcessor_bigevents_skim()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_bigevents_skim::init(void)
{
  // This is called once at program startup. If you are creating
  // and filling historgrams in this plugin, you should lock the
  // ROOT mutex like this:
  //
  // japp->RootWriteLock();
  //  ... fill historgrams or trees ...
  // japp->RootUnLock();
  //
  
  return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_bigevents_skim::brun(JEventLoop *eventLoop, int32_t runnumber)
{
  // This is called whenever the run number changes
  return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_bigevents_skim::evnt(JEventLoop *loop, uint64_t eventnumber)
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
  
  const DEventWriterEVIO* locEventWriterEVIO = NULL;
  loop->GetSingle(locEventWriterEVIO);
  // write out BOR events
  if(loop->GetJEvent().GetStatusBit(kSTATUS_BOR_EVENT)) {
    locEventWriterEVIO->Write_EVIOEvent(loop, "bigevents");
    return NOERROR;
  }
  // write out EPICS events
  if(loop->GetJEvent().GetStatusBit(kSTATUS_EPICS_EVENT)) {
    locEventWriterEVIO->Write_EVIOEvent(loop, "bigevents");
    return NOERROR;
  }
  
  // get trigger types
  const DL1Trigger *trig_words = NULL;
  uint32_t trig_mask, fp_trig_mask;
  try {
    loop->GetSingle(trig_words);
  } catch(...) {};
  if (trig_words) {
    trig_mask = trig_words->trig_mask;
    fp_trig_mask = trig_words->fp_trig_mask;
  }
  else {
    trig_mask = 0;
    fp_trig_mask = 0;
  }
  
  if (!(trig_mask & 0x4)){
    return NOERROR;
  }

  vector <const DCDCDigiHit*> CDCHits;
  loop->Get(CDCHits);
  
  if (CDCHits.size() > 2000) { // huge event lets keep it!

    locEventWriterEVIO->Write_EVIOEvent(loop, "bigevents");

  }



  return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_bigevents_skim::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_bigevents_skim::fini(void)
{
  // Called before program exit after event processing is finished.
  return NOERROR;
}

