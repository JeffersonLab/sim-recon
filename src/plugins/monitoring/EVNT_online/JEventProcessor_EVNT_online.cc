// $Id$
//
//    File: JEventProcessor_EVNT_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>
#include <random>


#include <DAQ/JEventSource_EVIO.h>

#include "JEventProcessor_EVNT_online.h"
#include <JANA/JApplication.h>


using namespace std;
using namespace jana;


// for random numbers
// static mt19937 generator;
// static normal_distribution<float> normal_dist(150.,150.);

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>


// root hist pointers
static TH1I * evntdata;
static TH2I * rocdata;



//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_EVNT_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_EVNT_online::JEventProcessor_EVNT_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_EVNT_online::~JEventProcessor_EVNT_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_EVNT_online::init(void) {

  // create root folder for evnt and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("evnt")->cd();


  // book hist
  evntdata = new TH1I("evntdata","Total data words in event",100,0,30000);
  rocdata  = new TH2I("rocdata","ROC vs data words",80,0.5,80.5,2000,0.,2000.);


  // back to main dir
  main->cd();

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_EVNT_online::brun(JEventLoop *eventLoop, int32_t runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_EVNT_online::evnt(JEventLoop *eventLoop, uint64_t eventnumber) {

  int nword,ntot;
  // uint32_t *buff;
  // uint32_t buff_size;


  // Get pointer to JEventSource base class
  JEvent &jevent = eventLoop->GetJEvent();
  JEventSource *source = jevent.GetJEventSource();
  
  // Cast source to type JEventSource_EVIO. You could use a
  // dynamic_cast here but that sometimes has problems when
  // used through a plugin. It is also notoriously inefficient so
  // checking a string is probably no worse.
  if( string("JEventSource_EVIO") == source->className() ) {
    JEventSource_EVIO *eviosource = (JEventSource_EVIO*)source;
    
    // Get EVIO buffer pointer and size (yes, I know, size is redundant)
    //    eviosource->GetEVIOBuffer(jevent,buff,buff_size);
    
    // Get evioDOMTree pointer and list of data banks
    evioDOMTree *dom = eviosource->GetEVIODOMTree(jevent);
    evioDOMNodeListP bankList = dom->getNodeList([](evioDOMNodeP n) {
        return (n->tag==1)&&(n->getContentType()==0x1)&&(n->getParent()->tag>0)&&(n->getSize()>0);
      });

    // loop over data banks and hist bank sizes, etc.
	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
    ntot=0;
    for(auto bank : *bankList.get()) {
      nword=bank->getSize();
      rocdata->Fill(bank->getParent()->tag,nword);
      ntot+=nword;
    }
    evntdata->Fill(ntot);
	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK
  }
  
  
  // fake data for the moment...
  // japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK
  // ntot=0;
  // for(auto roc=1; roc<=60; roc++) {
  //   nword = 25+normal_dist(generator);
  //   rochist->Fill(roc,nword);
  //   ntot+=nword;
  // }
  // evntsize->Fill(ntot);
  // japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_EVNT_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_EVNT_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
