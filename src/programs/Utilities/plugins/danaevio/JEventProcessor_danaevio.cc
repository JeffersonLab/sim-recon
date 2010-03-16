// JEventProcessor_danaevio.cc
//
//
// JANA plugin writes out DANA objects to EVIO files
//
//
// Implements JANA command-line parameters:
//
//    EVIOFILENAME     output file name, default "dana_events.evio"
//    EVIOBUFSIZE      serialized event internal buffer size, default 100000 words
//    DANAEVIO         specify which objects to serialize, see below for defaults
//    WRITEOUT         specify which objects to serialize, see below for defaults
//
//
//  Specifying which objects to serialize:
//    
//    Comma-separated, case-insensitive list (best not to included embedded whitespace)
//    Accepts "all", "none", "truth" and "hits", as well as individual DANA object names
//    Use prefix "-" to invert selection, "+" also accepted
//
//
//  dana_evio_dict.xml is corresponding evio2xml dictionary
//
//
//
//  Elliott Wolin, 9-Feb-2010
//
//
// still to do:
//    bank tags definitions
//    optimize multi-threading
//    add evio to external packages
//    associated objects, tag and num, etc.
//    tagged factories



#include <map>

#include <JANA/JEventProcessor.h>
#include <JANA/JFactory_base.h>
#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include "DANA/DApplication.h"

#include "TRACKING/DMCThrown.h"
#include "TRACKING/DMCTrackHit.h"
#include "TRACKING/DTrackWireBased.h"
#include "TRACKING/DTrackTimeBased.h"
#include "TRACKING/DMCTrajectoryPoint.h"

#include "FCAL/DFCALTruthShower.h"
#include "FCAL/DFCALHit.h"

#include "BCAL/DBCALTruthShower.h"
#include "BCAL/DHDDMBCALHit.h"

#include "TOF/DTOFTruth.h"
#include "TOF/DHDDMTOFHit.h"

#include "CDC/DCDCHit.h"

#include "FDC/DFDCHit.h"

#include "PID/DBeamPhoton.h"
#include "PID/DPhoton.h"
#include "PID/DChargedTrack.h"

#include "START_COUNTER/DSCHit.h"
#include "START_COUNTER/DSCTruthHit.h"

#include "SplitString.h"

#include "evioFileChannel.hxx"
#include "evioUtil.hxx"



using namespace std;
using namespace jana;
using namespace evio;



// evio output file name, use EVIOFILENAME command-line parameter to override
static string evioFileName = "dana_events.evio";


// internal evio buffer size, use EVIOBUFSIZE command-line parameter to override
static int evioBufSize=200000;


// list of dana objects that can be written out in evio format including default output flag
// use -PDANAEVIO or -PWRITEOUT override
// *** NOTE:  if you add to this list be sure to modify decode_object_parameters() appropriately ***
static pair<string,bool> danaObs[] =  {
  pair<string,bool> ("dmctrackhit",          false),
  pair<string,bool> ("dbeamphoton",          true),
  pair<string,bool> ("dmcthrown",            true),
  pair<string,bool> ("dfcaltruthshower",     true),
  pair<string,bool> ("dbcaltruthshower",     true),
  pair<string,bool> ("dtoftruth",            true),
  pair<string,bool> ("dsctruthhit",          true),
  pair<string,bool> ("dmctrajectorypoint",   false),
  pair<string,bool> ("dcdchit",              true),
  pair<string,bool> ("dfdchit",              true),
  pair<string,bool> ("dfcalhit",             true),
  pair<string,bool> ("dhddmbcalhit",         true),
  pair<string,bool> ("dhddmtofhit",          true),
  pair<string,bool> ("dschit",               true),
  pair<string,bool> ("dtrackwirebased",      false),
  pair<string,bool> ("dtracktimebased",      false),
  pair<string,bool> ("dchargedtrack",        false),
  pair<string,bool> ("dphoton",              false),
  pair<string,bool> ("dcdctrackhit",         false),
  pair<string,bool> ("dfdcpseudo",           false),
};
static map<string,bool> evioMap(danaObs,danaObs+sizeof(danaObs)/sizeof(danaObs[0]));
static map<string,bool> processedMap(danaObs,danaObs+sizeof(danaObs)/sizeof(danaObs[0]));


// indices to associated objects
static map<int,int> emptyMap;
static pair<string, map<int,int> > idPairs[] = {
  pair<string, map<int,int> > ("dmctrackhit",        emptyMap),
  pair<string, map<int,int> > ("dbeamphoton",        emptyMap),
  pair<string, map<int,int> > ("dmcthrown",          emptyMap),
  pair<string, map<int,int> > ("dfcalthruthshower",  emptyMap),
  pair<string, map<int,int> > ("dbcalthruthshower",  emptyMap),
  pair<string, map<int,int> > ("dtoftruth",          emptyMap),
  pair<string, map<int,int> > ("dsctruthhit",        emptyMap),
  pair<string, map<int,int> > ("dmctrajectorypoint", emptyMap),
  pair<string, map<int,int> > ("dcdchit",            emptyMap),
  pair<string, map<int,int> > ("dfdchit",            emptyMap),
  pair<string, map<int,int> > ("dfcalhit",           emptyMap),
  pair<string, map<int,int> > ("dhddmbcalhit",       emptyMap),
  pair<string, map<int,int> > ("dhddmtofhit",        emptyMap),
  pair<string, map<int,int> > ("dschit",             emptyMap),
  pair<string, map<int,int> > ("dtrackwirebased",    emptyMap),
  pair<string, map<int,int> > ("dtracktimebased",    emptyMap),
  pair<string, map<int,int> > ("dchargedtrack",      emptyMap),
  pair<string, map<int,int> > ("dphoton",            emptyMap),
  pair<string, map<int,int> > ("dcdctrackhit",       emptyMap),
  pair<string, map<int,int> > ("dfdcpseudo",         emptyMap),
};
static map<string, map<int,int> > idMap(idPairs,idPairs+sizeof(idPairs)/sizeof(idPairs[0]));


// bank tag definitions (totally arbitrary at the moment)
static pair<string,int> tagPairs[] = {
  pair<string,int> ("danaevent",            10000),
  pair<string,int> ("dmctrackhit",          10100),
  pair<string,int> ("dbeamphoton",          10200),
  pair<string,int> ("dmcthrown",            10300),
  pair<string,int> ("dfcaltruthshower",     10400),
  pair<string,int> ("dbcaltruthshower",     10500),
  pair<string,int> ("dtoftruth",            10600),
  pair<string,int> ("dsctruthhit",          10700),
  pair<string,int> ("dmctrajectorypoint",   10800),
  pair<string,int> ("dcdchit",              10900),
  pair<string,int> ("dfdchit",              11000),
  pair<string,int> ("dfcalhit",             11100),
  pair<string,int> ("dhddmbcalhit",         11200),
  pair<string,int> ("dhddmtofhit",          11300),
  pair<string,int> ("dschit",               11400),
  pair<string,int> ("dtracktimebased",      11500),
  pair<string,int> ("dchargedtrack",        11600),
  pair<string,int> ("dphoton",              11700),
  pair<string,int> ("dtrackwirebased",      11800),
  pair<string,int> ("dcdctrackhit",         11900),
  pair<string,int> ("dfdcpseudo",           12000),
};
static map<string,int> tagMap(tagPairs,tagPairs+sizeof(tagPairs)/sizeof(tagPairs[0]));


// mutex for serializing writing to file
static pthread_mutex_t evioMutex = PTHREAD_MUTEX_INITIALIZER;



//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


class JEventProcessor_danaevio : public JEventProcessor {
  

private:

  // file channel handle
  evioFileChannel *chan;
  

//----------------------------------------------------------------------------


public:  
  const char* className(void){
    return "JEventProcessor_danaevio";
  }
  
  
//----------------------------------------------------------------------------
  

  JEventProcessor_danaevio() {
    

    // check for EVIOFILENAME output file name parameter
    gPARMS->SetDefaultParameter("EVIOFILENAME",evioFileName);
    jout << endl << "  EVIO output file name is " << evioFileName << endl << endl;


    // check for DANAEVIO and WRITEOUT parameters which control which objects to output
    decode_object_parameters();
    map<string,bool>::iterator iter;
    jout << endl << "  DANA object output flags:" << endl << endl;
    for(iter=evioMap.begin(); iter!=evioMap.end(); iter++) {
      jout << "     " << setiosflags(ios::left) << setw(25) << iter->first
           << setw(2) << iter->second << endl;
    }
    jout << endl << endl;


    // check for EVIOBUFSIZE internal buffer size parameter
    gPARMS->SetDefaultParameter("EVIOBUFSIZE",evioBufSize);
    jout << endl << "  EVIO internal buf size is " << evioBufSize << endl << endl;


    // create file channel and open file
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


  ~JEventProcessor_danaevio() {
    
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


private:
  
  // needed until bug fixed in jana framework ???
  jerror_t fini() {
    chan->close();
    return(NOERROR);
  }


//----------------------------------------------------------------------------



  jerror_t evnt(JEventLoop *eventLoop, int eventnumber) {
    

    // clear all object processed flags
    map<string,bool>::iterator pIter;
    for(pIter=processedMap.begin(); pIter!=processedMap.end(); pIter++) pIter->second=false;


    // clear all id maps
    map< string, map<int,int> >::iterator idIter;
    for(idIter=idMap.begin(); idIter!=idMap.end(); idIter++) idIter->second.clear();


    // create evio DOM tree
    evioDOMTree tree(tagMap["danaevent"],0);
  

    // add various banks to tree
    if(evioMap["dmctrackhit"           ])  addDMCTrackHit(eventLoop,tree);
    if(evioMap["dbeamphoton"           ])  addDBeamPhoton(eventLoop,tree);
    if(evioMap["dmcthrown"             ])  addDMCThrown(eventLoop,tree);
    if(evioMap["dfcaltruthshower"      ])  addDFCALTruthShower(eventLoop,tree);
    if(evioMap["dbcaltruthshower"      ])  addDBCALTruthShower(eventLoop,tree);
    if(evioMap["dtoftruth"             ])  addDTOFTruth(eventLoop,tree);
    if(evioMap["dsctruthhit"           ])  addDSCTruthHit(eventLoop,tree);
    if(evioMap["dmctrajectorypoint"    ])  addDMCTrajectoryPoint(eventLoop,tree);

    if(evioMap["dcdchit"               ])  addDCDCHit(eventLoop,tree);
    if(evioMap["dfdchit"               ])  addDFDCHit(eventLoop,tree);
    if(evioMap["dfcalhit"              ])  addDFCALHit(eventLoop,tree);
    if(evioMap["dhddmbcalhit"          ])  addDHDDMBCALHit(eventLoop,tree);
    if(evioMap["dhddmtofhit"           ])  addDHDDMTOFHit(eventLoop,tree);
    if(evioMap["dschit"                ])  addDSCHit(eventLoop,tree);
  
    if(evioMap["dcdctrackhit"          ])  addDCDCTrackHit(eventLoop,tree);
    if(evioMap["dfdcpseudo"            ])  addDFDCPseudo(eventLoop,tree);
    if(evioMap["dtrackwirebased"       ])  addDTrackWireBased(eventLoop,tree);
    if(evioMap["dtracktimebased"       ])  addDTrackTimeBased(eventLoop,tree);
    if(evioMap["dchargedtrack"         ])  addDChargedTrack(eventLoop,tree);
    if(evioMap["dphoton"               ])  addDPhoton(eventLoop,tree);

  
    // get lock, write out evio tree, unlock
    pthread_mutex_lock(&evioMutex);
    try {
      chan->write(tree);
    } catch (evioException e) {
      jerr << endl << "  ?evioException in JEventProcessor_danaevio::evnt" << endl << endl 
           << e.toString() << endl;
    } catch (...) {
      jerr << endl << "  ?unknown exception in JEventProcessor_danaevio::evnt, unable to write to file" << endl << endl;
    }
    pthread_mutex_unlock(&evioMutex);

  
    // done
    return NOERROR;
  }


//----------------------------------------------------------------------------


  void addDMCThrown(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DMCThrown*> mcthrowns; 
    eventLoop->Get(mcthrowns); 
    if(mcthrowns.size()<=0)return;


    // create mcthrown bank and add to event tree
    evioDOMNodeP mcthrown = evioDOMNode::createEvioDOMNode(tagMap["dmcthrown"],0);
    tree << mcthrown;


    // create data banks and add to mcthrown
    evioDOMNodeP typeBank      =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dmcthrown"],1);   
    evioDOMNodeP pdgtypeBank   =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dmcthrown"],2);
    evioDOMNodeP myidBank      =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dmcthrown"],3);   
    evioDOMNodeP parentidBank  =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dmcthrown"],4);
    evioDOMNodeP mechBank      =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dmcthrown"],5);   
    evioDOMNodeP xBank         =  evioDOMNode::createEvioDOMNode<float>(tagMap["dmcthrown"],6);      
    evioDOMNodeP yBank         =  evioDOMNode::createEvioDOMNode<float>(tagMap["dmcthrown"],7);      
    evioDOMNodeP zBank         =  evioDOMNode::createEvioDOMNode<float>(tagMap["dmcthrown"],8);      
    evioDOMNodeP pxBank        =  evioDOMNode::createEvioDOMNode<float>(tagMap["dmcthrown"],9);     
    evioDOMNodeP pyBank        =  evioDOMNode::createEvioDOMNode<float>(tagMap["dmcthrown"],10);    
    evioDOMNodeP pzBank        =  evioDOMNode::createEvioDOMNode<float>(tagMap["dmcthrown"],11);    
    evioDOMNodeP energyBank    =  evioDOMNode::createEvioDOMNode<float>(tagMap["dmcthrown"],12);
    *mcthrown << typeBank << pdgtypeBank << myidBank << parentidBank << mechBank 
              << xBank << yBank << zBank << pxBank << pyBank << pzBank << energyBank;


    // add track data to banks
    for(unsigned int i=0; i<mcthrowns.size(); i++) {
      *typeBank     << mcthrowns[i]->type;
      *pdgtypeBank  << mcthrowns[i]->pdgtype;
      *myidBank     << mcthrowns[i]->myid;
      *parentidBank << mcthrowns[i]->parentid;
      *mechBank     << mcthrowns[i]->mech;

      DVector3 pos = mcthrowns[i]->position();
      *xBank        << pos.X();
      *yBank        << pos.Y();
      *zBank        << pos.Z();

      DVector3 mom = mcthrowns[i]->momentum();
      *pxBank       << mom.X();
      *pyBank       << mom.Y();
      *pzBank       << mom.Z();
      
      *energyBank   << mcthrowns[i]->energy();

      idMap["dmcthrown"][mcthrowns[i]->id]=i;
    }

    // done
    processedMap["dmcthrown"]=true;
  }



//------------------------------------------------------------------------------


  void addDMCTrackHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DMCTrackHit*> mctrackhits; 
    eventLoop->Get(mctrackhits); 
    if(mctrackhits.size()<=0)return;


    // create mctrackhit bank and add to event tree
    evioDOMNodeP mctrackhit = evioDOMNode::createEvioDOMNode(tagMap["dmctrackhit"],0);
    tree << mctrackhit;


    // create data banks and add to mctrackhit bank
    evioDOMNodeP rBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dmctrackhit"],1);
    evioDOMNodeP phiBank     = evioDOMNode::createEvioDOMNode<float>(tagMap["dmctrackhit"],2);
    evioDOMNodeP zBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dmctrackhit"],3);
    evioDOMNodeP trackBank   = evioDOMNode::createEvioDOMNode  <int>(tagMap["dmctrackhit"],4);
    evioDOMNodeP primaryBank = evioDOMNode::createEvioDOMNode  <int>(tagMap["dmctrackhit"],5);
    evioDOMNodeP ptypeBank   = evioDOMNode::createEvioDOMNode  <int>(tagMap["dmctrackhit"],6);
    evioDOMNodeP systemBank  = evioDOMNode::createEvioDOMNode  <int>(tagMap["dmctrackhit"],7);
    *mctrackhit << rBank << phiBank << zBank << trackBank << primaryBank << ptypeBank << systemBank;


    // add track data to banks
    for(unsigned int i=0; i<mctrackhits.size(); i++) {
      *rBank        << mctrackhits[i]->r;
      *phiBank      << mctrackhits[i]->phi;
      *zBank        << mctrackhits[i]->z;
      *trackBank    << mctrackhits[i]->track;
      *primaryBank  << mctrackhits[i]->primary;
      *ptypeBank    << mctrackhits[i]->ptype;
      *systemBank   << mctrackhits[i]->system;

      idMap["dmctrackhit"][mctrackhits[i]->id]=i;
    }


    // done
    processedMap["dmctrackhit"]=true;
  }


//------------------------------------------------------------------------------


  void addDTOFTruth(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DTOFTruth*> toftruths; 
    eventLoop->Get(toftruths); 
    if(toftruths.size()<=0)return;


    // create toftruth bank and add to event tree
    evioDOMNodeP toftruth = evioDOMNode::createEvioDOMNode(tagMap["dtoftruth"],0);
    tree << toftruth;


    // create data banks and add to toftruth bank
    evioDOMNodeP trackBank   =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dtoftruth"],1);
    evioDOMNodeP primaryBank =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dtoftruth"],2);
    evioDOMNodeP xBank       =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],3);
    evioDOMNodeP yBank       =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],4);
    evioDOMNodeP zBank       =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],5);
    evioDOMNodeP pxBank      =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],6);
    evioDOMNodeP pyBank      =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],7);
    evioDOMNodeP pzBank      =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],8);
    evioDOMNodeP tBank       =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],9);
    evioDOMNodeP EBank       =  evioDOMNode::createEvioDOMNode<float>(tagMap["dtoftruth"],10);
    evioDOMNodeP ptypeBank   =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dtoftruth"],11);
    *toftruth << trackBank << primaryBank << xBank << yBank << zBank << pxBank << pyBank << pzBank
              << tBank << EBank << ptypeBank;

    
    // add track data to banks
    for(unsigned int i=0; i<toftruths.size(); i++) {

      *trackBank   << toftruths[i]->track;
      *primaryBank << toftruths[i]->primary;
      *xBank       << toftruths[i]->x;
      *yBank       << toftruths[i]->y;
      *zBank       << toftruths[i]->z;
      *pxBank      << toftruths[i]->px;
      *pyBank      << toftruths[i]->py;
      *pzBank      << toftruths[i]->pz;
      *tBank       << toftruths[i]->t;
      *EBank       << toftruths[i]->E;
      *ptypeBank   << toftruths[i]->ptype;

      idMap["dtoftruth"][toftruths[i]->id]=i;
    }


    // done
    processedMap["dtoftruth"]=true;
  }


//------------------------------------------------------------------------------



  void addDFCALTruthShower(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DFCALTruthShower*> fcaltruthshowers; 
    eventLoop->Get(fcaltruthshowers); 
    if(fcaltruthshowers.size()<=0)return;


    // create fcaltruthshower bank and add to event tree
    evioDOMNodeP fcaltruthshower = evioDOMNode::createEvioDOMNode(tagMap["dfcaltruthshower"],0);
    tree << fcaltruthshower;


    // create data banks and add to fcaltruthshower
    evioDOMNodeP xBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],1);
    evioDOMNodeP yBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],2);
    evioDOMNodeP zBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],3);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],4);
    evioDOMNodeP pxBank      = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],5);
    evioDOMNodeP pyBank      = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],6);
    evioDOMNodeP pzBank      = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],7);
    evioDOMNodeP EBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],8);
    evioDOMNodeP primaryBank = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],9);
    evioDOMNodeP trackBank   = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],10);
    evioDOMNodeP typeBank    = evioDOMNode::createEvioDOMNode<float>(tagMap["dfcaltruthshower"],11);
    *fcaltruthshower << xBank << yBank << zBank << tBank << pxBank << pyBank<< pzBank<< EBank 
                     << primaryBank<< trackBank<< typeBank;


    // add track data to banks
    for(unsigned int i=0; i<fcaltruthshowers.size(); i++) {
      *xBank       << fcaltruthshowers[i]->x();
      *yBank       << fcaltruthshowers[i]->y();
      *zBank       << fcaltruthshowers[i]->z();
      *tBank       << fcaltruthshowers[i]->t();
      *pxBank      << fcaltruthshowers[i]->px();
      *pyBank      << fcaltruthshowers[i]->py();
      *pzBank      << fcaltruthshowers[i]->pz();
      *EBank       << fcaltruthshowers[i]->E();
      *primaryBank << fcaltruthshowers[i]->primary();
      *trackBank   << fcaltruthshowers[i]->track();
      *typeBank    << fcaltruthshowers[i]->type();

      idMap["dfcaltruthshower"][fcaltruthshowers[i]->id]=i;
    }


    // done
    processedMap["dfcaltruthshower"]=true;
  }


  //------------------------------------------------------------------------------


  void addDBCALTruthShower(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DBCALTruthShower*> bcaltruthshowers; 
    eventLoop->Get(bcaltruthshowers); 
    if(bcaltruthshowers.size()<=0)return;


    // create bcaltruthshower bank and add to event tree
    evioDOMNodeP bcaltruthshower = evioDOMNode::createEvioDOMNode(tagMap["dbcaltruthshower"],0);
    tree << bcaltruthshower;


    // create data banks and add to bcaltruthshower
    evioDOMNodeP trackBank    =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dbcaltruthshower"],1);
    evioDOMNodeP primaryBank  =  evioDOMNode::createEvioDOMNode<int>  (tagMap["dbcaltruthshower"],2);
    evioDOMNodeP phiBank      =  evioDOMNode::createEvioDOMNode<float>(tagMap["dbcaltruthshower"],3);
    evioDOMNodeP rBank        =  evioDOMNode::createEvioDOMNode<float>(tagMap["dbcaltruthshower"],4);
    evioDOMNodeP zBank        =  evioDOMNode::createEvioDOMNode<float>(tagMap["dbcaltruthshower"],5);
    evioDOMNodeP tBank        =  evioDOMNode::createEvioDOMNode<float>(tagMap["dbcaltruthshower"],6);
    evioDOMNodeP EBank        =  evioDOMNode::createEvioDOMNode<float>(tagMap["dbcaltruthshower"],7);
    *bcaltruthshower << trackBank << primaryBank << phiBank << rBank << zBank << tBank << EBank;


    // add track data to banks
    for(unsigned int i=0; i<bcaltruthshowers.size(); i++) {
      *trackBank   << bcaltruthshowers[i]->track;
      *primaryBank << bcaltruthshowers[i]->primary;
      *phiBank     << bcaltruthshowers[i]->phi;
      *rBank       << bcaltruthshowers[i]->r;
      *zBank       << bcaltruthshowers[i]->z;
      *tBank       << bcaltruthshowers[i]->t;
      *EBank       << bcaltruthshowers[i]->E;
      
      idMap["dbcaltruthshower"][bcaltruthshowers[i]->id]=i;
    }


    // done
    processedMap["dbcaltruthshower"]=true;
  }


//------------------------------------------------------------------------------


  void addDCDCHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DCDCHit*> cdchits; 
    eventLoop->Get(cdchits); 
    if(cdchits.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP cdchit = evioDOMNode::createEvioDOMNode(tagMap["dcdchit"],0);
    tree << cdchit;


    // create data banks and add to cdchit bank
    evioDOMNodeP ringBank  = evioDOMNode::createEvioDOMNode<int>  (tagMap["dcdchit"],1);
    evioDOMNodeP strawBank = evioDOMNode::createEvioDOMNode<int>  (tagMap["dcdchit"],2);
    evioDOMNodeP dEBank    = evioDOMNode::createEvioDOMNode<float>(tagMap["dcdchit"],3);
    evioDOMNodeP tBank     = evioDOMNode::createEvioDOMNode<float>(tagMap["dcdchit"],4);
    *cdchit << ringBank << strawBank << dEBank << tBank;


    // add track data to banks
    for(unsigned int i=0; i<cdchits.size(); i++) {
      *ringBank  << cdchits[i]->ring;
      *strawBank << cdchits[i]->straw;
      *dEBank    << cdchits[i]->dE;
      *tBank     << cdchits[i]->t;

      idMap["dcdchit"][cdchits[i]->id]=i;
    }


    // done
    processedMap["dcdchit"]=true;
  }


//------------------------------------------------------------------------------


  void addDMCTrajectoryPoint(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DMCTrajectoryPoint*> mctrajectorypoints; 
    eventLoop->Get(mctrajectorypoints); 
    if(mctrajectorypoints.size()<=0)return;


    // create bank and add to event tree
    evioDOMNodeP mctrajectorypoint = evioDOMNode::createEvioDOMNode(tagMap["dmctrajectorypoint"],0);
    tree << mctrajectorypoint;


    // create data banks and add to bank
    evioDOMNodeP xBank              = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],1);
    evioDOMNodeP yBank              = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],2);
    evioDOMNodeP zBank              = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],3);
    evioDOMNodeP pxBank             = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],4);
    evioDOMNodeP pyBank             = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],5);
    evioDOMNodeP pzBank             = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],6);
    evioDOMNodeP EBank              = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],7);
    evioDOMNodeP dEBank             = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],8);
    evioDOMNodeP primary_trackBank  = evioDOMNode::createEvioDOMNode<int>   (tagMap["dmctrajectorypoint"],9);
    evioDOMNodeP trackBank          = evioDOMNode::createEvioDOMNode<int>   (tagMap["dmctrajectorypoint"],10);
    evioDOMNodeP partBank           = evioDOMNode::createEvioDOMNode<int>   (tagMap["dmctrajectorypoint"],11);
    evioDOMNodeP radlenBank         = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],12);
    evioDOMNodeP stepBank           = evioDOMNode::createEvioDOMNode<float> (tagMap["dmctrajectorypoint"],13);
    evioDOMNodeP mechBank           = evioDOMNode::createEvioDOMNode<int>   (tagMap["dmctrajectorypoint"],14);
    *mctrajectorypoint << xBank << yBank << zBank << pxBank << pyBank << pzBank << EBank << dEBank
                       << primary_trackBank << trackBank << partBank << radlenBank << stepBank << mechBank;


    // add track data to banks
    for(unsigned int i=0; i<mctrajectorypoints.size(); i++) {
      *xBank              << mctrajectorypoints[i]->x;
      *yBank              << mctrajectorypoints[i]->y;
      *zBank              << mctrajectorypoints[i]->z;
      *pxBank             << mctrajectorypoints[i]->px;
      *pyBank             << mctrajectorypoints[i]->py;
      *pzBank             << mctrajectorypoints[i]->pz;
      *EBank              << mctrajectorypoints[i]->E;
      *dEBank             << mctrajectorypoints[i]->dE;
      *primary_trackBank  << mctrajectorypoints[i]->primary_track;
      *trackBank          << mctrajectorypoints[i]->track;
      *partBank           << mctrajectorypoints[i]->part;
      *radlenBank         << mctrajectorypoints[i]->radlen;
      *stepBank           << mctrajectorypoints[i]->step;
      *mechBank           << mctrajectorypoints[i]->mech;

      idMap["dmctrajectorypoint"][mctrajectorypoints[i]->id]=i;
    }


    // done
    processedMap["dmctrajectorypoint"]=true;
  }


//------------------------------------------------------------------------------


  void addDFDCHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DFDCHit*> fdchits; 
    eventLoop->Get(fdchits); 
    if(fdchits.size()<=0)return;


    // create fdchit bank and add to event tree
    evioDOMNodeP fdchit = evioDOMNode::createEvioDOMNode(tagMap["dfdchit"],0);
    tree << fdchit;


    // create data banks and add to fdchit bank
    evioDOMNodeP layerBank   = evioDOMNode::createEvioDOMNode<int>  (tagMap["dfdchit"],1);
    evioDOMNodeP moduleBank  = evioDOMNode::createEvioDOMNode<int>  (tagMap["dfdchit"],2);
    evioDOMNodeP elementBank = evioDOMNode::createEvioDOMNode<int>  (tagMap["dfdchit"],3);
    evioDOMNodeP planeBank   = evioDOMNode::createEvioDOMNode<int>  (tagMap["dfdchit"],4);
    evioDOMNodeP gPlaneBank  = evioDOMNode::createEvioDOMNode<int>  (tagMap["dfdchit"],5);
    evioDOMNodeP gLayerBank  = evioDOMNode::createEvioDOMNode<int>  (tagMap["dfdchit"],6);
    evioDOMNodeP qBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfdchit"],7);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfdchit"],8);
    evioDOMNodeP rBank       = evioDOMNode::createEvioDOMNode<float>(tagMap["dfdchit"],9);
    evioDOMNodeP typeBank    = evioDOMNode::createEvioDOMNode<int>  (tagMap["dfdchit"],10);
    *fdchit << layerBank << moduleBank << elementBank << planeBank << gPlaneBank << gLayerBank
            << qBank << tBank << rBank << typeBank;


    // add track data to banks
    for(unsigned int i=0; i<fdchits.size(); i++) {
      *layerBank    << fdchits[i]->layer;
      *moduleBank   << fdchits[i]->module;
      *elementBank  << fdchits[i]->element;
      *planeBank    << fdchits[i]->plane;
      *gPlaneBank   << fdchits[i]->gPlane;
      *gLayerBank   << fdchits[i]->gLayer;
      *qBank        << fdchits[i]->q;
      *tBank        << fdchits[i]->t;
      *rBank        << fdchits[i]->r;
      *typeBank     << fdchits[i]->type;

      idMap["dfdchit"][fdchits[i]->id]=i;
    }


    // done
    processedMap["fdchit"]=true;
  }


//----------------------------------------------------------------------------


  void addDBeamPhoton(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DBeamPhoton*> beamphotons; 
    eventLoop->Get(beamphotons); 
    if(beamphotons.size()<=0)return;


    // create bank and add to event tree
    evioDOMNodeP beamphoton = evioDOMNode::createEvioDOMNode(tagMap["dbeamphoton"],0);
    tree << beamphoton;


    // create data banks and add to bank
    evioDOMNodeP xBank   = evioDOMNode::createEvioDOMNode<float>  (tagMap["dbeamphoton"],1);
    evioDOMNodeP yBank   = evioDOMNode::createEvioDOMNode<float>  (tagMap["dbeamphoton"],2);
    evioDOMNodeP zBank   = evioDOMNode::createEvioDOMNode<float>  (tagMap["dbeamphoton"],3);
    evioDOMNodeP pxBank  = evioDOMNode::createEvioDOMNode<float>  (tagMap["dbeamphoton"],4);
    evioDOMNodeP pyBank  = evioDOMNode::createEvioDOMNode<float>  (tagMap["dbeamphoton"],5);
    evioDOMNodeP pzBank  = evioDOMNode::createEvioDOMNode<float>  (tagMap["dbeamphoton"],6);
    evioDOMNodeP tBank   = evioDOMNode::createEvioDOMNode<int>    (tagMap["dbeamphoton"],7);
    *beamphoton << xBank<< yBank<< zBank<< pxBank<< pyBank<< pzBank << tBank;


    // add track data to banks
    for(unsigned int i=0; i<beamphotons.size(); i++) {
      DVector3 pos = beamphotons[i]->position();
      *xBank   <<  pos.X();
      *yBank   <<  pos.Y();
      *zBank   <<  pos.Z();
      
      DVector3 mom = beamphotons[i]->momentum();
      *pxBank  <<  mom.X();
      *pyBank  <<  mom.Y();
      *pzBank  <<  mom.Z();

      *tBank  << beamphotons[i]->t;

      idMap["dbeamphoton"][beamphotons[i]->id]=i;
    }


    // done
    processedMap["dbeamphoton"]=true;
  }


//------------------------------------------------------------------------------


  void addDSCTruthHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DSCTruthHit*> sctruthhits;
    eventLoop->Get(sctruthhits); 
    if(sctruthhits.size()<=0)return;


    // create bank and add to event tree
    evioDOMNodeP sctruthhit = evioDOMNode::createEvioDOMNode(tagMap["dsctruthhit"],0);
    tree << sctruthhit;


    // create data banks and add to bank
    evioDOMNodeP dEdxBank     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dsctruthhit"],1);
    evioDOMNodeP primaryBank  = evioDOMNode::createEvioDOMNode<int8_t> (tagMap["dsctruthhit"],2);
    evioDOMNodeP trackBank    = evioDOMNode::createEvioDOMNode<int>    (tagMap["dsctruthhit"],3);
    evioDOMNodeP ptypeBank    = evioDOMNode::createEvioDOMNode<int>    (tagMap["dsctruthhit"],4);
    evioDOMNodeP rBank        = evioDOMNode::createEvioDOMNode<float>  (tagMap["dsctruthhit"],5);
    evioDOMNodeP phiBank      = evioDOMNode::createEvioDOMNode<float>  (tagMap["dsctruthhit"],6);
    evioDOMNodeP zBank        = evioDOMNode::createEvioDOMNode<float>  (tagMap["dsctruthhit"],7);
    evioDOMNodeP tBank        = evioDOMNode::createEvioDOMNode<float>  (tagMap["dsctruthhit"],8);
    evioDOMNodeP sectorBank   = evioDOMNode::createEvioDOMNode<int>    (tagMap["dsctruthhit"],9);
    *sctruthhit << dEdxBank << primaryBank << trackBank << ptypeBank << rBank << phiBank 
                << zBank << tBank << sectorBank;


    // add track data to banks
    for(unsigned int i=0; i<sctruthhits.size(); i++) {
      *dEdxBank     << sctruthhits[i]->dEdx;
      *primaryBank  << (int8_t)sctruthhits[i]->primary;
      *trackBank    << sctruthhits[i]->track;
      *ptypeBank    << sctruthhits[i]->ptype;
      *rBank        << sctruthhits[i]->r;
      *phiBank      << sctruthhits[i]->phi;
      *zBank        << sctruthhits[i]->z;
      *tBank        << sctruthhits[i]->t;
      *sectorBank   << sctruthhits[i]->sector;

      idMap["dsctruthhit"][sctruthhits[i]->id]=i;
    }


    // done
    processedMap["dsctruthhit"]=true;
  }


//------------------------------------------------------------------------------


  void addDFCALHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DFCALHit*> fcalhits;
    eventLoop->Get(fcalhits); 
    if(fcalhits.size()<=0)return;


    // create bank and add to event tree
    evioDOMNodeP fcalhit = evioDOMNode::createEvioDOMNode(tagMap["dfcalhit"],0);
    tree << fcalhit;


    // create data banks and add to bank
    evioDOMNodeP rowBank     = evioDOMNode::createEvioDOMNode<int>   (tagMap["dfcalhit"],1);
    evioDOMNodeP columnBank  = evioDOMNode::createEvioDOMNode<int>   (tagMap["dfcalhit"],2);
    evioDOMNodeP xBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dfcalhit"],3);
    evioDOMNodeP yBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dfcalhit"],4);
    evioDOMNodeP EBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dfcalhit"],5);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dfcalhit"],6);
    *fcalhit << rowBank << columnBank << xBank << yBank << EBank << tBank;


    // add track data to banks
    for(unsigned int i=0; i<fcalhits.size(); i++) {
      *rowBank     << fcalhits[i]->row;
      *columnBank  << fcalhits[i]->column;
      *xBank       << fcalhits[i]->x;
      *yBank       << fcalhits[i]->y;
      *EBank       << fcalhits[i]->E;
      *tBank       << fcalhits[i]->t;

      idMap["dfcalhit"][fcalhits[i]->id]=i;
    }


    // done
    processedMap["dfcalhit"]=true;
  }


//------------------------------------------------------------------------------


  void addDHDDMBCALHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DHDDMBCALHit*> hddmbcalhits;
    eventLoop->Get(hddmbcalhits); 
    if(hddmbcalhits.size()<=0)return;


    // create hddm bcal bank and add to event tree
    evioDOMNodeP hddmbcalhit = evioDOMNode::createEvioDOMNode(tagMap["dhddmbcalhit"],0);
    tree << hddmbcalhit;


    // create data banks and add to bank
    evioDOMNodeP moduleBank  = evioDOMNode::createEvioDOMNode<int>   (tagMap["dhddmbcalhit"],1);
    evioDOMNodeP layerBank   = evioDOMNode::createEvioDOMNode<int>   (tagMap["dhddmbcalhit"],2);
    evioDOMNodeP sectorBank  = evioDOMNode::createEvioDOMNode<int>   (tagMap["dhddmbcalhit"],3);
    evioDOMNodeP EBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmbcalhit"],4);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmbcalhit"],5);
    evioDOMNodeP zLocalBank  = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmbcalhit"],6);
    *hddmbcalhit << moduleBank << layerBank << sectorBank << EBank << tBank << zLocalBank;


    // add track data to banks
    for(unsigned int i=0; i<hddmbcalhits.size(); i++) {
      *moduleBank  << hddmbcalhits[i]->module;
      *layerBank   << hddmbcalhits[i]->layer;
      *sectorBank  << hddmbcalhits[i]->sector;
      *EBank       << hddmbcalhits[i]->E;
      *tBank       << hddmbcalhits[i]->t;
      *zLocalBank  << hddmbcalhits[i]->zLocal;

      idMap["dhddmbcalhit"][hddmbcalhits[i]->id]=i;
    }


    // done
    processedMap["dhddmbcalhit"]=true;
  }


//------------------------------------------------------------------------------


  void addDHDDMTOFHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DHDDMTOFHit*> hddmtofhits;
    eventLoop->Get(hddmtofhits); 
    if(hddmtofhits.size()<=0)return;


    // create bank and add to event tree
    evioDOMNodeP hddmtofhit = evioDOMNode::createEvioDOMNode(tagMap["dhddmtofhit"],0);
    tree << hddmtofhit;


    // create data banks and add to bank
    evioDOMNodeP planeBank    = evioDOMNode::createEvioDOMNode<int>   (tagMap["dhddmtofhit"],1);
    evioDOMNodeP barBank      = evioDOMNode::createEvioDOMNode<int>   (tagMap["dhddmtofhit"],2);
    evioDOMNodeP ptypeBank    = evioDOMNode::createEvioDOMNode<int>   (tagMap["dhddmtofhit"],3);
    evioDOMNodeP t_northBank  = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],4);
    evioDOMNodeP dE_northBank = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],5);
    evioDOMNodeP t_southBank  = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],6);
    evioDOMNodeP dE_southBank = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],7);
    evioDOMNodeP xBank        = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],8);
    evioDOMNodeP yBank        = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],9);
    evioDOMNodeP zBank        = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],10);
    evioDOMNodeP pxBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],11);
    evioDOMNodeP pyBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],12);
    evioDOMNodeP pzBank       = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],13);
    evioDOMNodeP EBank        = evioDOMNode::createEvioDOMNode<float> (tagMap["dhddmtofhit"],14);
    *hddmtofhit << planeBank << barBank << ptypeBank << t_northBank << dE_northBank 
                << t_southBank << dE_southBank << xBank << yBank << zBank 
                << pxBank << pyBank << pzBank<< EBank;


    // add track data to banks
    for(unsigned int i=0; i<hddmtofhits.size(); i++) {
      *planeBank    << hddmtofhits[i]->plane;
      *barBank      << hddmtofhits[i]->bar;
      *ptypeBank    << hddmtofhits[i]->ptype;
      *t_northBank  << hddmtofhits[i]->t_north;
      *dE_northBank << hddmtofhits[i]->dE_north;
      *t_southBank  << hddmtofhits[i]->t_south;
      *dE_southBank << hddmtofhits[i]->dE_south;
      *xBank        << hddmtofhits[i]->x;
      *yBank        << hddmtofhits[i]->y;
      *zBank        << hddmtofhits[i]->z;
      *pxBank       << hddmtofhits[i]->px;
      *pyBank       << hddmtofhits[i]->py;
      *pzBank       << hddmtofhits[i]->pz;
      *EBank        << hddmtofhits[i]->E;

      idMap["dhddmtofhit"][hddmtofhits[i]->id]=i;
    }


    // done
    processedMap["dhddmtofhit"]=true;
  }


//------------------------------------------------------------------------------


  void addDSCHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DSCHit*> schits;
    eventLoop->Get(schits); 
    if(schits.size()<=0)return;


    // create bank and add to event tree
    evioDOMNodeP schit = evioDOMNode::createEvioDOMNode(tagMap["dschit"],0);
    tree << schit;


    // create data banks and add to bank
    evioDOMNodeP dEBank      = evioDOMNode::createEvioDOMNode<float>  (tagMap["dschit"],1);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float>  (tagMap["dschit"],2);
    evioDOMNodeP sectorBank  = evioDOMNode::createEvioDOMNode<int>    (tagMap["dschit"],3);
    *schit << dEBank << tBank << sectorBank;


    // add track data to banks
    for(unsigned int i=0; i<schits.size(); i++) {
      *dEBank      << schits[i]->dE;
      *tBank       << schits[i]->t;
      *sectorBank  << schits[i]->sector;

      idMap["dschit"][schits[i]->id]=i;
    }


    // done
    processedMap["dschit"]=true;
  }


//------------------------------------------------------------------------------


  void addDTrackWireBased(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DTrackWireBased*> wirebasedtracks;
    eventLoop->Get(wirebasedtracks); 
    if(wirebasedtracks.size()<=0)return;


    // create wirebasedtrack bank and add to event tree
    evioDOMNodeP wirebasedtrack = evioDOMNode::createEvioDOMNode(tagMap["dtrackwirebased"],0);
    tree << wirebasedtrack;


    // create data banks and add to bank (n.b. time based track has FOM, wire based doesn't)
    evioDOMNodeP chisq = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],1);
    evioDOMNodeP Ndof  = evioDOMNode::createEvioDOMNode<int>    (tagMap["dtrackwirebased"],2);
    evioDOMNodeP x     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],4);
    evioDOMNodeP y     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],5);
    evioDOMNodeP z     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],6);
    evioDOMNodeP px    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],7);
    evioDOMNodeP py    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],8);
    evioDOMNodeP pz    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],9);
    evioDOMNodeP q     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],10);
    evioDOMNodeP E     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],11);
    evioDOMNodeP mass  = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtrackwirebased"],12);
    *wirebasedtrack << chisq << Ndof << x << y << z << px << py << pz 
                    << q << E << mass;


    // add track data to banks
    for(unsigned int i=0; i<wirebasedtracks.size(); i++) {
      *chisq      << wirebasedtracks[i]->chisq;
      *Ndof       << wirebasedtracks[i]->Ndof;
      *x          << wirebasedtracks[i]->x();
      *y          << wirebasedtracks[i]->y();
      *z          << wirebasedtracks[i]->z();
      *px         << wirebasedtracks[i]->px();
      *py         << wirebasedtracks[i]->py();
      *pz         << wirebasedtracks[i]->pz();
      *q          << wirebasedtracks[i]->charge();
      *E          << wirebasedtracks[i]->energy();
      *mass       << wirebasedtracks[i]->mass();
      idMap["dtrackwirebased"][wirebasedtracks[i]->id]=i;
    }


    // done
    processedMap["dtrackwirebased"]=true;
  }


//------------------------------------------------------------------------------


  void addDTrackTimeBased(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DTrackTimeBased*> timebasedtracks;
    eventLoop->Get(timebasedtracks); 
    if(timebasedtracks.size()<=0)return;


    // create timebasedtrack bank and add to event tree
    evioDOMNodeP timebasedtrack = evioDOMNode::createEvioDOMNode(tagMap["dtracktimebased"],0);
    tree << timebasedtrack;


    // create data banks and add to bank
    evioDOMNodeP chisq = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],1);
    evioDOMNodeP Ndof  = evioDOMNode::createEvioDOMNode<int>    (tagMap["dtracktimebased"],2);
    evioDOMNodeP FOM   = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],3);
    evioDOMNodeP x     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],4);
    evioDOMNodeP y     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],5);
    evioDOMNodeP z     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],6);
    evioDOMNodeP px    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],7);
    evioDOMNodeP py    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],8);
    evioDOMNodeP pz    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],9);
    evioDOMNodeP q     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],10);
    evioDOMNodeP E     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],11);
    evioDOMNodeP mass  = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],12);
    evioDOMNodeP t0    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dtracktimebased"],13);
    *timebasedtrack << chisq << Ndof << FOM << x << y << z << px << py << pz 
                    << q << E << mass << t0;


    // add track data to banks
    for(unsigned int i=0; i<timebasedtracks.size(); i++) {
      *chisq      << timebasedtracks[i]->chisq;
      *Ndof       << timebasedtracks[i]->Ndof;
      *FOM        << timebasedtracks[i]->FOM;
      *x          << timebasedtracks[i]->x();
      *y          << timebasedtracks[i]->y();
      *z          << timebasedtracks[i]->z();
      *px         << timebasedtracks[i]->px();
      *py         << timebasedtracks[i]->py();
      *pz         << timebasedtracks[i]->pz();
      *q          << timebasedtracks[i]->charge();
      *E          << timebasedtracks[i]->energy();
      *mass       << timebasedtracks[i]->mass();
      *t0         << timebasedtracks[i]->t0();
      idMap["dtracktimebased"][timebasedtracks[i]->id]=i;
    }


    //  ??? must check if index is available ???
    // add associated object banks
    

    // DTrackTimeBased
    evioDOMNodeP wireBasedIndex = evioDOMNode::createEvioDOMNode<int>(tagMap["dtracktimebased"],100);
    *timebasedtrack << wireBasedIndex;
    for(unsigned int i=0; i<timebasedtracks.size(); i++) {
      vector<const DTrackWireBased*> wireTracks;
      timebasedtracks[i]->GetT(wireTracks); 
      for(unsigned int j=0; j<wireTracks.size(); j++) *wireBasedIndex << idMap["dtrackwirebased"][wireTracks[j]->id];
    }


    // DCDCTrackHit
    evioDOMNodeP cdctrackhitIndexBank  = evioDOMNode::createEvioDOMNode(tagMap["dtracktimebased"],101);
    *timebasedtrack << cdctrackhitIndexBank;
    for(unsigned int i=0; i<timebasedtracks.size(); i++) {

      evioDOMNodeP cdctrackhitIndex = evioDOMNode::createEvioDOMNode<int> (101,i);  // ??? tag,num ???
      *cdctrackhitIndexBank << cdctrackhitIndex;

      vector<const DCDCTrackHit*> cdctrackhits;
      timebasedtracks[i]->GetT(cdctrackhits); 
      for(unsigned int j=0; j<cdctrackhits.size(); j++) {
        map<int,int>::iterator iter = idMap["dcdctrackhit"].find(cdctrackhits[j]->id);
        if(iter==idMap["dcdctrackhit"].end())cerr << "unable to find cdctrackhit id " << cdctrackhits[j]->id << endl;
        *cdctrackhitIndex << idMap["dcdctrackhit"][cdctrackhits[j]->id];
      }
    }


    // DFDCPseudo
    evioDOMNodeP fdcpseudoIndexBank  = evioDOMNode::createEvioDOMNode(tagMap["dtracktimebased"],102);
    *timebasedtrack << fdcpseudoIndexBank;
    for(unsigned int i=0; i<timebasedtracks.size(); i++) {

      evioDOMNodeP fdcpseudoIndex = evioDOMNode::createEvioDOMNode<int> (101,i);
      *fdcpseudoIndexBank << fdcpseudoIndex;

      vector<const DFDCPseudo*> fdcpseudos;
      timebasedtracks[i]->GetT(fdcpseudos); 
      for(unsigned int j=0; j<fdcpseudos.size(); j++) {
        map<int,int>::iterator iter = idMap["dfdcpseudo"].find(fdcpseudos[j]->id);
        if(iter==idMap["dfdcpseudo"].end())cerr << "unable to find fdcpseudo id " << fdcpseudos[j]->id << endl;
        *fdcpseudoIndex << idMap["dfdcpseudo"][fdcpseudos[j]->id];
      }
    }


    // done
    processedMap["dtracktimebased"]=true;
  }


//------------------------------------------------------------------------------


  void addDChargedTrack(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DChargedTrack*> chargedtracks;
    eventLoop->Get(chargedtracks); 
    if(chargedtracks.size()<=0)return;


    // create chargedtrack bank and add to event tree
    evioDOMNodeP chargedtrack = evioDOMNode::createEvioDOMNode(tagMap["dchargedtrack"],0);
    tree << chargedtrack;

    // create index bank for each charged track and add to chargedtrack bank
    for(unsigned int i=0; i<chargedtracks.size(); i++) {
      evioDOMNodeP hypotheses = evioDOMNode::createEvioDOMNode<int> (tagMap["dchargedtrack"],1);
      *chargedtrack << hypotheses;
      for(unsigned int j=0; j<chargedtracks[i]->hypotheses.size(); j++) {
        *hypotheses<< idMap["dtracktimebased"][chargedtracks[i]->hypotheses[j]->id];
      }
      idMap["dchargedtrack"][chargedtracks[i]->id]=i;
    }


    // done
    processedMap["dchargedtrack"]=true;
  }


//------------------------------------------------------------------------------


  void addDPhoton(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DPhoton*> photons;
    eventLoop->Get(photons); 
    if(photons.size()<=0)return;


    // create photon bank and add to event tree
    evioDOMNodeP photon = evioDOMNode::createEvioDOMNode(tagMap["dphoton"],0);
    tree << photon;


    // create data banks and add to photon bank
    evioDOMNodeP E     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],1);
    evioDOMNodeP px    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],2);
    evioDOMNodeP py    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],3);
    evioDOMNodeP pz    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],4);
    evioDOMNodeP x     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],5);
    evioDOMNodeP y     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],6);
    evioDOMNodeP z     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],7);
    evioDOMNodeP t     = evioDOMNode::createEvioDOMNode<float>  (tagMap["dphoton"],8);
    evioDOMNodeP Tag   = evioDOMNode::createEvioDOMNode<int>    (tagMap["dphoton"],9);
    *photon << E << px << py << pz << x << y << z << t << Tag;


    // add track data to banks
    for(unsigned int i=0; i<photons.size(); i++) {
      *E          << photons[i]->energy();
      *px         << photons[i]->px();
      *py         << photons[i]->py();
      *pz         << photons[i]->pz();
      *x          << photons[i]->x();
      *y          << photons[i]->y();
      *z          << photons[i]->z();
      *t          << photons[i]->getTime();
      *Tag        << photons[i]->getTag();

      idMap["dphoton"][photons[i]->id]=i;
    }


    // done
    processedMap["dphoton"]=true;
  }


//------------------------------------------------------------------------------


  void addDCDCTrackHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DCDCTrackHit*> cdctrackhits; 
    eventLoop->Get(cdctrackhits); 
    if(cdctrackhits.size()<=0)return;


    // create cdctrackhit bank and add to event tree
    evioDOMNodeP cdctrackhit = evioDOMNode::createEvioDOMNode(tagMap["dcdctrackhit"],0);
    tree << cdctrackhit;


    // create data banks and add to cdctrackhit bank
    evioDOMNodeP ringBank    = evioDOMNode::createEvioDOMNode<int>    (tagMap["dcdctrackhit"],1);
    evioDOMNodeP strawBank   = evioDOMNode::createEvioDOMNode<int>    (tagMap["dcdctrackhit"],2);
    evioDOMNodeP xBank       = evioDOMNode::createEvioDOMNode<float>  (tagMap["dcdctrackhit"],3);
    evioDOMNodeP yBank       = evioDOMNode::createEvioDOMNode<float>  (tagMap["dcdctrackhit"],4);
    evioDOMNodeP stereoBank  = evioDOMNode::createEvioDOMNode<float>  (tagMap["dcdctrackhit"],5);
    evioDOMNodeP tdriftBank  = evioDOMNode::createEvioDOMNode<float>  (tagMap["dcdctrackhit"],6);
    evioDOMNodeP distBank    = evioDOMNode::createEvioDOMNode<float>  (tagMap["dcdctrackhit"],7);
    evioDOMNodeP dEBank      = evioDOMNode::createEvioDOMNode<float>  (tagMap["dcdctrackhit"],8);
    *cdctrackhit << ringBank << strawBank << xBank << yBank << stereoBank << tdriftBank << distBank << dEBank;


    // add track data to banks
    for(unsigned int i=0; i<cdctrackhits.size(); i++) {
      *ringBank    << cdctrackhits[i]->wire->ring;
      *strawBank   << cdctrackhits[i]->wire->straw;
      *xBank       << cdctrackhits[i]->wire->origin.x();
      *yBank       << cdctrackhits[i]->wire->origin.y();
      *stereoBank  << cdctrackhits[i]->wire->stereo;
      *tdriftBank  << cdctrackhits[i]->tdrift;
      *distBank    << cdctrackhits[i]->dist;
      *dEBank      << cdctrackhits[i]->dE;

      idMap["dcdctrackhit"][cdctrackhits[i]->id]=i;
    }


    // done
    processedMap["dcdctrackhit"]=true;
  }


//------------------------------------------------------------------------------


  void addDFDCPseudo(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DFDCPseudo*> fdcpseudos; 
    eventLoop->Get(fdcpseudos); 
    if(fdcpseudos.size()<=0)return;


    // create fdcpseudo bank and add to event tree
    evioDOMNodeP fdcpseudo = evioDOMNode::createEvioDOMNode(tagMap["dfdcpseudo"],0);
    tree << fdcpseudo;


    // create data banks and add to fdcpseudo bank
    evioDOMNodeP wBank         = evioDOMNode::createEvioDOMNode<float> (tagMap["dfdcpseudo"],1);
    evioDOMNodeP sBank         = evioDOMNode::createEvioDOMNode<float> (tagMap["dfdcpseudo"],2);
    evioDOMNodeP layerBank     = evioDOMNode::createEvioDOMNode<int>   (tagMap["dfdcpseudo"],3);
    evioDOMNodeP wireBank      = evioDOMNode::createEvioDOMNode<int>   (tagMap["dfdcpseudo"],4);
    evioDOMNodeP timeBank      = evioDOMNode::createEvioDOMNode<float> (tagMap["dfdcpseudo"],5);
    evioDOMNodeP distBank      = evioDOMNode::createEvioDOMNode<float> (tagMap["dfdcpseudo"],6);
    evioDOMNodeP statusBank    = evioDOMNode::createEvioDOMNode<int>   (tagMap["dfdcpseudo"],7);
    evioDOMNodeP xBank         = evioDOMNode::createEvioDOMNode<float> (tagMap["dfdcpseudo"],8);
    evioDOMNodeP yBank         = evioDOMNode::createEvioDOMNode<float> (tagMap["dfdcpseudo"],9);
    evioDOMNodeP dEBank        = evioDOMNode::createEvioDOMNode<float> (tagMap["dfdcpseudo"],10);
    *fdcpseudo << wBank << sBank << layerBank << wireBank << timeBank << distBank
               << statusBank << xBank << yBank << dEBank;


    // add track data to banks
    for(unsigned int i=0; i<fdcpseudos.size(); i++) {
      *wBank        << fdcpseudos[i]->w;
      *sBank        << fdcpseudos[i]->s;
      *layerBank    << fdcpseudos[i]->wire->layer;
      *wireBank     << fdcpseudos[i]->wire->wire;
      *timeBank     << fdcpseudos[i]->time;
      *distBank     << fdcpseudos[i]->dist;
      *statusBank   << fdcpseudos[i]->status;
      *xBank        << fdcpseudos[i]->x;
      *yBank        << fdcpseudos[i]->y;
      *dEBank       << fdcpseudos[i]->dE;

      idMap["dfdcpseudo"][fdcpseudos[i]->id]=i;
    }


    // done
    processedMap["fdcpseudo"]=true;
  }


//----------------------------------------------------------------------------


  void decode_object_parameters(void) {
    
    // params are comma-separated and case-insensitive
    // "-" means invert
    // "+" is ignored
    // also supported: "all", "none", "truth" "hits", "tracks"
    // otherwise parameter must be the name of a DANA object that is processed by this program
    

    // first check for DANAEVIO tag, then WRITEOUT tag
    map<string,bool>::iterator iter;   
    string danaevio= "";
    gPARMS->SetDefaultParameter("DANAEVIO",danaevio);
    if(danaevio=="") gPARMS->SetDefaultParameter("WRITEOUT",danaevio);

    if(danaevio!="") {

      vector<string> params;
      SplitString<string>(danaevio,params,",");
      for(unsigned int i=0; i<params.size(); i++) {
        std::transform(params[i].begin(), params[i].end(), params[i].begin(), (int(*)(int)) tolower);
        bool plus=(params[i][0]=='+');
        bool minus=(params[i][0]=='-');
        string value=params[i].substr((plus||minus)?1:0);

        if(value=="all") {
          for(iter=evioMap.begin(); iter!=evioMap.end(); iter++) iter->second=!minus;

        } else if(value=="none") {
          for(iter=evioMap.begin(); iter!=evioMap.end(); iter++) iter->second=minus;

        } else if(value=="truth") {
          evioMap["dbeamphoton"]=!minus;
          evioMap["dmcthrown"]=!minus;
          evioMap["dmctrackhit"]=!minus;
          evioMap["dfcaltruthshower"]=!minus;
          evioMap["dbcaltruthshower"]=!minus;
          evioMap["dtoftruth"]=!minus;
          evioMap["dsctruth"]=!minus;

        } else if(value=="hits") {
          evioMap["dcdchit"]=!minus;
          evioMap["dfdchit"]=!minus;
          evioMap["dfcalchit"]=!minus;
          evioMap["dhddmbcalchit"]=!minus;
          evioMap["dhddmtofchit"]=!minus;
          evioMap["dschit"]=!minus;

        } else if(value=="tracks") {
          evioMap["dtrackwirebased"]=!minus;
          evioMap["dtracktimebased"]=!minus;
          evioMap["dchargedtrack"]=!minus;
          evioMap["dphoton"]=!minus;
          evioMap["dcdctrackhit"]=!minus;
          evioMap["dfdcpseudo"]=!minus;

        } else {
          map<string,bool>::iterator found = evioMap.find(value);
          if(found!=evioMap.end()) {
            found->second=!minus;
          } else {
            jerr << endl << "  ?unknown DANAEVIO or WRITEOUT parameter: " << params[i] << endl;
          }

        }
      }
    }


    // if DChargedTrack requested then DTrackTimeBased must be also be present
    if(evioMap["dchargedtrack"])evioMap["dtracktimebased"]=true;
    
  }


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------
};


// for initializing plugin
extern "C" {
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_danaevio());
  }
} // "extern C"


//----------------------------------------------------------------------------
