// JEventProcessor_danaevio.cc
//
//
// JANA event processor plugin writes out evio DOM tree to file
//
//
// Implements JANA command-line parameters:
//
//    EVIO::FILENAME     output file name, default "dana_events.evio"
//    EVIO::BUFSIZE      serialized event internal buffer size, default 200000 words
//
//
//  dana_evio_dict.xml is corresponding evio2xml dictionary
//
//  E.g. to run:
//     $HALLD_HOME/bin/Linux_CentOS5-x86_64-gcc4.1.2/hd_ana --plugin=danaevio -PEVIO:DANAEVIO="all" ../Event.hddm
//
//
//  Elliott Wolin, 22-Feb-2010


#include <DDANAEVIODOMTree.h>
#include <DDANAEVIO_factory.h>
#include <JEventProcessor_danaevio.h>
#include <dana_evio_dict.h>


// evio output file name, use EVIO:FILENAME command-line parameter to override
static string evioFileName = "dana_events.evio";


// internal evio buffer size, use EVIO:BUFSIZE command-line parameter to override
static int evioBufSize=200000;


// mutex for serializing writing to file
static pthread_mutex_t evioMutex = PTHREAD_MUTEX_INITIALIZER;



//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


JEventProcessor_danaevio::JEventProcessor_danaevio() {
    

  jout << endl << "  Default JEventProcessor_danaevio invoked" << endl << endl;


  // check for EVIO:FILENAME output file name parameter
  gPARMS->SetDefaultParameter("EVIO:FILENAME",evioFileName);
  jout << endl << "  EVIO output file name is " << evioFileName << endl << endl;
  
  
  // check for EVIO:BUFSIZE internal buffer size parameter
  gPARMS->SetDefaultParameter("EVIO:BUFSIZE",evioBufSize);
  jout << endl << "  EVIO internal buf size is " << evioBufSize << endl << endl;
  
  
  // create and open file channel
  try {
    chan = new evioFileChannel(evioFileName,"w",evioBufSize);
    chan->open();

  } catch (evioException e) {
    jerr << endl << "  ?evioException in JEventProcessor_danaevio" << endl << endl 
         << e.toString() << endl;
  } catch (...) {
    jerr << endl << "  ?unknown exception in JEventProcessor_danaevio, unable to open output file" << endl << endl;
  }
}  


//----------------------------------------------------------------------------


JEventProcessor_danaevio::~JEventProcessor_danaevio() {
    
  // close file and delete file channel object
  try {
    
    chan->close();
    
  } catch (evioException e) {
    jerr << endl << "  ?evioException in ~JEventProcessor_danaevio" << endl << endl 
         << e.toString() << endl;
  } catch (...) {
    jerr << endl << "  ?unknown exception in ~JEventProcessor_danaevio, unable to close output file" << endl << endl;
  }
  
  delete(chan);
}


//----------------------------------------------------------------------------


jerror_t JEventProcessor_danaevio::brun(JEventLoop *eventLoop, int runnumber) {

  static bool first_time = true;


  // get write lock
  pthread_mutex_lock(&evioMutex);
  
  // create dictionary banks from DDANAEVIO factory tagMap<string, pair<uint16_t,uint8_t> >
  //  and write out as first event in file
  if(first_time) {
    first_time=false;
    
    evioDOMTree tree(1,0);
    evioDOMNodeP name =  evioDOMNode::createEvioDOMNode<string>   (1,1);
    evioDOMNodeP tag  =  evioDOMNode::createEvioDOMNode<uint16_t> (1,2);
    evioDOMNodeP num  =  evioDOMNode::createEvioDOMNode<uint8_t>  (1,3);
    tree << name << tag << num;

    map< string, pair<uint16_t,uint8_t> > *theMap = DDANAEVIO_factory::getTagMapPointer();
    map< string, pair<uint16_t,uint8_t> >::iterator iter;
    for(iter=theMap->begin(); iter!=theMap->end(); iter++) {
      *name << iter->first;
      *tag  << iter->second.first;
      *num  << iter->second.second;
    }    

    chan->write(tree);
  }
  
  // unlock
  pthread_mutex_unlock(&evioMutex);
  
  return(NOERROR);
}


//----------------------------------------------------------------------------


jerror_t JEventProcessor_danaevio::evnt(JEventLoop *eventLoop, int eventnumber) {
    

  // get evio trees
  vector<const DDANAEVIODOMTree*> evioTrees; 
  eventLoop->Get(evioTrees);
  if(evioTrees.size()<=0)return(NOERROR);


  // get write lock
  pthread_mutex_lock(&evioMutex);


  // write out all evio trees
  for(unsigned int i=0; i<evioTrees.size(); i++) {
    try {
      chan->write(evioTrees[i]->tree);

    } catch (evioException e) {
      jerr << endl << "  ?evioException in JEventProcessor_danaevio::evnt" << endl << endl 
           << e.toString() << endl;
    } catch (...) {
      jerr << endl << "  ?unknown exception in JEventProcessor_danaevio::evnt, unable to write to file" << endl << endl;
    }
  }


  // unlock
  pthread_mutex_unlock(&evioMutex);

  
  // done
  return NOERROR;
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
