// $Id$
//
//    File: JEventProcessor_rawevent.cc
// Created: Fri Jun 24 12:05:19 EDT 2011
// Creator: wolin (on Linux stan.jlab.org 2.6.18-194.11.1.el5 x86_64)
//

// JANA event processor plugin "rawevent" is run as follows:
//
//     $ hd_ana --plugin=rawevent input_file.hddm
//
//
//  ejw, 27-jun-2011



#include "JEventProcessor_rawevent.h"


static pthread_mutex_t rawMutex = PTHREAD_MUTEX_INITIALIZER;


// for evio output file
static evioFileChannel *chan = NULL;
static int evioBufSize       = 750000;


// current run number
static int runNumber;



//----------------------------------------------------------------------------


// required for JANA framework to find and identify this plugin
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_rawevent());
}
} // "C"


//----------------------------------------------------------------------------


// JEventProcessor_rawevent (Constructor) invoked once only
JEventProcessor_rawevent::JEventProcessor_rawevent() {
}


//----------------------------------------------------------------------------


// ~JEventProcessor_rawevent (Destructor) once only
JEventProcessor_rawevent::~JEventProcessor_rawevent() {
}


//----------------------------------------------------------------------------


// init called once-only at beginning, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::init(void) {
	return NOERROR;
}

//----------------------------------------------------------------------------


// brun called once-only at beginning of run, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::brun(JEventLoop *eventLoop, int runnumber) {

  runNumber=runnumber;
  jout << endl << "   brun called for run " << runNumber << endl << endl;


  // close old and open new output file channel
  if(chan!=NULL)chan->close();
  chan = new evioFileChannel("fileName.evio","w",evioBufSize);
  chan->open();

  return NOERROR;
}


//----------------------------------------------------------------------------


// evnt called once per event, in every processing thread
// MUST be thread-safe!
jerror_t JEventProcessor_rawevent::evnt(JEventLoop *eventLoop, int eventnumber) {

  unsigned short tag;
  unsigned char num;


  // create evio output event tree
  evioDOMTree eventTree(tag=1,num=0);


  // create and fill header bank with vector of ints and add to tree
  evioDOMNodeP head =  evioDOMNode::createEvioDOMNode<int>(tag=1,num=1);
  vector<int> v(10,1);
  *head << v;
  eventTree << head;


  // get DTOFRAWHIT banks
  vector<const DTOFRawHit*> dtofrawhits; 
  eventLoop->Get(dtofrawhits);
  for(unsigned int i=0; i<dtofrawhits.size(); i++) {
    int bar   = dtofrawhits[i]->bar;
    int plane = dtofrawhits[i]->plane;
    int lr    = dtofrawhits[i]->lr;
    float dE  = dtofrawhits[i]->dE;
    float t   = dtofrawhits[i]->t;

    // now do something with hit information
  }


  // also need to process DBCALHIT,DFCALHIT,DFDCHIT,DCDCHIT,DSCHIT


  // lock and write to file
  pthread_mutex_lock(&rawMutex);
  chan->write(eventTree);
  pthread_mutex_unlock(&rawMutex);


  return NOERROR;
}


//----------------------------------------------------------------------------


// erun called once-only at end of run, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.

  jout << endl << "   erun called for run " << runNumber << endl << endl;


  // close evio output file
  if(chan!=NULL)chan->close();
  chan=NULL;

  return NOERROR;
}


//----------------------------------------------------------------------------


// fini called once-only when done, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::fini(void) {
	// Called before program exit after event processing is finished.
	return NOERROR;
}


//----------------------------------------------------------------------------
