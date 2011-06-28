// $Id$
//
//    File: JEventProcessor_rawevent.cc
// Created: Fri Jun 24 12:05:19 EDT 2011
// Creator: wolin (on Linux stan.jlab.org 2.6.18-194.11.1.el5 x86_64)
//
//
//
// JANA event processor plugin "rawevent"
//
//  Gets raw hit data from hddm file, converts to (crate,slot,channel), then creates
//   simulated raw output evio event and writes to file.
//
//  Opens new file for every run.
//  Default output file name is rawevent_xxxxxx.evio where xxxxxx is the run number.
//
//
// To run:
//
//     $ hd_ana --plugin=rawevent inputFile.hddm
//
//
//
// still to do:
//    translation table
//
//
//
//  ejw, 27-jun-2011



#include "JEventProcessor_rawevent.h"

#include<sstream>
#include<iomanip>



// to protect writing to output file
static pthread_mutex_t rawMutex = PTHREAD_MUTEX_INITIALIZER;


// for evio output
static evioFileChannel *chan   = NULL;
static int evioBufSize         = 750000;
static string fileBase         = "rawevent";
static string outputFileName;


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

  // get fileBase from command line params
  gPARMS->SetDefaultParameter("RAW:FILEBASE",fileBase);

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


  // close old output file
  if(chan!=NULL) {
    chan->close();
    delete(chan);
    chan=NULL;
  }


  // get new file name
  stringstream ss;
  ss << fileBase << "_" << setw(6) << setfill('0') << runNumber << ".evio" << ends;
  outputFileName=ss.str();


  // open new output file
  chan = new evioFileChannel(outputFileName,"w",evioBufSize);
  chan->open();
  jout << endl << "   opening output file:   " << outputFileName << endl << endl << endl;


  // add header event if required
  // ...


  return NOERROR;
}


//----------------------------------------------------------------------------


// called once per event in many different processing threads, so:
//
//    *** MUST be thread-safe ***
//
jerror_t JEventProcessor_rawevent::evnt(JEventLoop *eventLoop, int eventnumber) {

  unsigned short tag;
  unsigned char num;
  unsigned int i;


  // create evio output event tree in CODA raw data format
  evioDOMTree eventTree(tag=1,num=0);



  // example how to create and fill bank and add to event tree
  evioDOMNodeP testBank =  evioDOMNode::createEvioDOMNode<int>(tag=1,num=1);
  *testBank << vector<int>(10,1);
  eventTree << testBank;



  // DTOFRawHit
  vector<const DTOFRawHit*> dtofrawhits; 
  eventLoop->Get(dtofrawhits);
  for(i=0; i<dtofrawhits.size(); i++) {
    int bar   = dtofrawhits[i]->bar;
    int plane = dtofrawhits[i]->plane;
    int lr    = dtofrawhits[i]->lr;
    float dE  = dtofrawhits[i]->dE;
    float t   = dtofrawhits[i]->t;

    // translate to crate/slot/channel
    cscVal csc  = DTOFRawHitTranslation(bar,plane,lr);
    int crate   = csc.get<0>();
    int slot    = csc.get<1>();
    int channel = csc.get<2>();

    // do something...
  }


  // DBCALHit
  vector<const DBCALHit*> dbcalhits;
  eventLoop->Get(dbcalhits);
  for(i=0; i<dbcalhits.size(); i++) {
    int module   = dbcalhits[i]->module;
    int layer    = dbcalhits[i]->layer;
    int sector   = dbcalhits[i]->sector;
    int end      = dbcalhits[i]->end;       // 0 for upstream
    float E      = dbcalhits[i]->E;
    float t      = dbcalhits[i]->t;

    // translate to crate/slot/channel
    cscVal csc  = DBCALHitTranslation(module,layer,sector,end);
    int crate   = csc.get<0>();
    int slot    = csc.get<1>();
    int channel = csc.get<2>();

    // do something...
  }      



  // DFCALHit
  vector<const DFCALHit*> dfcalhits;
  eventLoop->Get(dfcalhits);
  for(i=0; i<dfcalhits.size(); i++) {
    int row      = dfcalhits[i]->row;
    int column   = dfcalhits[i]->column;
    float x      = dfcalhits[i]->x;
    float y      = dfcalhits[i]->y;
    float E      = dfcalhits[i]->E;
    float t      = dfcalhits[i]->t;
    
    // translate to crate/slot/channel
    cscVal csc  = DFCALHitTranslation(row,column,x,y);
    int crate   = csc.get<0>();
    int slot    = csc.get<1>();
    int channel = csc.get<2>();

    // do something...
  }      


  // DFDCHit
  vector<const DFDCHit*> dfdchits; 
  eventLoop->Get(dfdchits);
  for(unsigned int i=0; i<dfdchits.size(); i++) {
    int layer    = dfdchits[i]->layer;
    int module   = dfdchits[i]->module;
    int element  = dfdchits[i]->element;
    int plane    = dfdchits[i]->plane;
    int gPlane   = dfdchits[i]->gPlane;
    int gLayer   = dfdchits[i]->gLayer;
    float q      = dfdchits[i]->q;
    float t      = dfdchits[i]->t;
    float r      = dfdchits[i]->r;
    int type     = dfdchits[i]->type;

    // translate to crate/slot/channel
    cscVal csc  = DFDCHitTranslation(layer,module,element,plane,gPlane,gLayer);
    int crate   = csc.get<0>();
    int slot    = csc.get<1>();
    int channel = csc.get<2>();

    // do something...
  }


  // DFDCHit
  vector<const DCDCHit*> dcdchits; 
  eventLoop->Get(dcdchits);
  for(i=0; i<dcdchits.size(); i++) {
    int ring     = dcdchits[i]->ring;
    int straw    = dcdchits[i]->straw;
    float dE     = dcdchits[i]->dE;
    float t      = dcdchits[i]->t;

    // translate to crate/slot/channel
    cscVal csc  = DCDCHitTranslation(ring,straw);
    int crate   = csc.get<0>();
    int slot    = csc.get<1>();
    int channel = csc.get<2>();

    // do something...
  }      


  // DSCHit
  vector<const DSCHit*> dschits;
  eventLoop->Get(dschits);
  for(unsigned int i=0; i<dschits.size(); i++) {
    float dE     = dschits[i]->dE;
    float t      = dschits[i]->t;
    int sector   = dschits[i]->sector;

    // translate to crate/slot/channel
    cscVal csc  = DSCHitTranslation(sector);
    int crate   = csc.get<0>();
    int slot    = csc.get<1>();
    int channel = csc.get<2>();

    // do something...
  }      


  // DTagger
  vector<const DTagger*> dtaggerhits;
  eventLoop->Get(dtaggerhits);
  for(i=0; i<dtaggerhits.size(); i++) {
    int row      = dtaggerhits[i]->row;
    int column   = dtaggerhits[i]->column;
    float E      = dtaggerhits[i]->E;
    float t      = dtaggerhits[i]->t;

    // translate to crate/slot/channel
    cscVal csc  = DTaggerTranslation(row,column);
    int crate   = csc.get<0>();
    int slot    = csc.get<1>();
    int channel = csc.get<2>();

    // do something...
  }      




  // construct evio banks from hit data collected earlier and add to event tree
  // ...




  // write out event tree:  get lock, write to file, unlock
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


  // add end event if required
  // ...


  // close evio output file and delete channel
  if(chan!=NULL) {
    chan->close();
    delete(chan);
    chan=NULL;
    jout << endl << "  output file " << outputFileName << " closed" << endl << endl;
  }

  return NOERROR;
}


//----------------------------------------------------------------------------


// fini called once-only when done, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::fini(void) {
  return NOERROR;
}


//----------------------------------------------------------------------------






//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
// the following routines access the translation table
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DTOFRawHitTranslation(int bar,int plane, int lr) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DBCALHitTranslation(int module,int layer,int sector,int end) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DFCALHitTranslation(int row,int column,float x,float y) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DFDCHitTranslation(int layer,int module,int element,int plane,int gPlane,int gLayer) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DCDCHitTranslation(int ring,int straw) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DSCHitTranslation(int sector) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DTaggerTranslation(int row,int column) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


