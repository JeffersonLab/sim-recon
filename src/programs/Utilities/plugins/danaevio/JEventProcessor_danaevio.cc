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
//    bank tags
//    optimize multi-threading
//    add evio to external packages
//    create makefile.evio
//    bool vs int8_t in evio?



#include <map>

#include <JANA/JEventProcessor.h>
#include <JANA/JFactory_base.h>
#include <JANA/JApplication.h>
#include <JANA/JEventLoop.h>

#include "DANA/DApplication.h"

#include "TRACKING/DMCThrown.h"
#include "TRACKING/DMCTrackHit.h"
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


// dana objects that can be written out in evio format and default output flag
// use DANAEVIO or WRITEOUT command-line parameters to override
// *** if you add to this list be sure to modify decode_object_parameters() appropriately ***
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
  pair<string,bool> ("dtracktimebased",      true),
  pair<string,bool> ("dchargedtrack",        true),
  pair<string,bool> ("dphoton",              true),
};
static map<string,bool> evioMap(danaObs,danaObs+sizeof(danaObs)/sizeof(danaObs[0]));


// to allow DChargedTrack to cross-index DTrackTimeBased
map<int,int> DTrackTimeBasedMap;


// evio bank tag definitions (totally arbitrary at the moment)
static const int danaeventTag            = 10000;

static const int dmctrackhitTag          = 10100;
static const int dbeamphotonTag          = 10200;
static const int dmcthrownTag            = 10300;
static const int dfcaltruthshowerTag     = 10400;
static const int dbcaltruthshowerTag     = 10500;
static const int dtoftruthTag            = 10600;
static const int dsctruthhitTag          = 10700;
static const int dmctrajectorypointTag   = 10800;

static const int dcdchitTag              = 10900;
static const int dfdchitTag              = 11000;
static const int dfcalhitTag             = 11100;
static const int dhddmbcalhitTag         = 11200;
static const int dhddmtofhitTag          = 11300;
static const int dschitTag               = 11400;
static const int dtracktimebasedTag      = 11500;
static const int dchargedtrackTag        = 11600;
static const int dphotonTag              = 11700;


// mutex needed to protect serialization and writing to file
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
    
    static int count = 0;
    count++;

  
    // create evio DOM tree
    evioDOMTree tree(danaeventTag,0);
  

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
    evioDOMNodeP mcthrown = evioDOMNode::createEvioDOMNode(dmcthrownTag,0);
    tree << mcthrown;


    // create data banks and add to mcthrown
    evioDOMNodeP typeBank      =  evioDOMNode::createEvioDOMNode<int>  (dmcthrownTag,1);   
    evioDOMNodeP pdgtypeBank   =  evioDOMNode::createEvioDOMNode<int>  (dmcthrownTag,2);
    evioDOMNodeP myidBank      =  evioDOMNode::createEvioDOMNode<int>  (dmcthrownTag,3);   
    evioDOMNodeP parentidBank  =  evioDOMNode::createEvioDOMNode<int>  (dmcthrownTag,4);
    evioDOMNodeP mechBank      =  evioDOMNode::createEvioDOMNode<int>  (dmcthrownTag,5);   
    evioDOMNodeP xBank         =  evioDOMNode::createEvioDOMNode<float>(dmcthrownTag,6);      
    evioDOMNodeP yBank         =  evioDOMNode::createEvioDOMNode<float>(dmcthrownTag,7);      
    evioDOMNodeP zBank         =  evioDOMNode::createEvioDOMNode<float>(dmcthrownTag,8);      
    evioDOMNodeP pxBank        =  evioDOMNode::createEvioDOMNode<float>(dmcthrownTag,9);     
    evioDOMNodeP pyBank        =  evioDOMNode::createEvioDOMNode<float>(dmcthrownTag,10);    
    evioDOMNodeP pzBank        =  evioDOMNode::createEvioDOMNode<float>(dmcthrownTag,11);    
    evioDOMNodeP energyBank    =  evioDOMNode::createEvioDOMNode<float>(dmcthrownTag,12);
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
    }
  }



//------------------------------------------------------------------------------


  void addDMCTrackHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DMCTrackHit*> mctrackhits; 
    eventLoop->Get(mctrackhits); 
    if(mctrackhits.size()<=0)return;


    // create mctrackhit bank and add to event tree
    evioDOMNodeP mctrackhit = evioDOMNode::createEvioDOMNode(dmctrackhitTag,0);
    tree << mctrackhit;


    // create data banks and add to mctrackhit bank
    evioDOMNodeP rBank       = evioDOMNode::createEvioDOMNode<float>(dmctrackhitTag,1);
    evioDOMNodeP phiBank     = evioDOMNode::createEvioDOMNode<float>(dmctrackhitTag,2);
    evioDOMNodeP zBank       = evioDOMNode::createEvioDOMNode<float>(dmctrackhitTag,3);
    evioDOMNodeP trackBank   = evioDOMNode::createEvioDOMNode  <int>(dmctrackhitTag,4);
    evioDOMNodeP primaryBank = evioDOMNode::createEvioDOMNode  <int>(dmctrackhitTag,5);
    evioDOMNodeP ptypeBank   = evioDOMNode::createEvioDOMNode  <int>(dmctrackhitTag,6);
    evioDOMNodeP systemBank  = evioDOMNode::createEvioDOMNode  <int>(dmctrackhitTag,7);
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
    }

  }


//------------------------------------------------------------------------------


  void addDTOFTruth(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DTOFTruth*> toftruths; 
    eventLoop->Get(toftruths); 
    if(toftruths.size()<=0)return;


    // create toftruth bank and add to event tree
    evioDOMNodeP toftruth = evioDOMNode::createEvioDOMNode(dtoftruthTag,0);
    tree << toftruth;


    // create data banks and add to toftruth bank
    evioDOMNodeP trackBank   =  evioDOMNode::createEvioDOMNode<int>  (dtoftruthTag,1);
    evioDOMNodeP primaryBank =  evioDOMNode::createEvioDOMNode<int>  (dtoftruthTag,2);
    evioDOMNodeP xBank       =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,3);
    evioDOMNodeP yBank       =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,4);
    evioDOMNodeP zBank       =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,5);
    evioDOMNodeP pxBank      =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,6);
    evioDOMNodeP pyBank      =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,7);
    evioDOMNodeP pzBank      =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,8);
    evioDOMNodeP tBank       =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,9);
    evioDOMNodeP EBank       =  evioDOMNode::createEvioDOMNode<float>(dtoftruthTag,10);
    evioDOMNodeP ptypeBank   =  evioDOMNode::createEvioDOMNode<int>  (dtoftruthTag,11);
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
    }


  }


//------------------------------------------------------------------------------



  void addDFCALTruthShower(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DFCALTruthShower*> fcaltruthshowers; 
    eventLoop->Get(fcaltruthshowers); 
    if(fcaltruthshowers.size()<=0)return;


    // create fcaltruthshower bank and add to event tree
    evioDOMNodeP fcaltruthshower = evioDOMNode::createEvioDOMNode(dfcaltruthshowerTag,0);
    tree << fcaltruthshower;


    // create data banks and add to fcaltruthshower
    evioDOMNodeP xBank       = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,1);
    evioDOMNodeP yBank       = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,2);
    evioDOMNodeP zBank       = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,3);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,4);
    evioDOMNodeP pxBank      = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,5);
    evioDOMNodeP pyBank      = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,6);
    evioDOMNodeP pzBank      = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,7);
    evioDOMNodeP EBank       = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,8);
    evioDOMNodeP primaryBank = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,9);
    evioDOMNodeP trackBank   = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,10);
    evioDOMNodeP typeBank    = evioDOMNode::createEvioDOMNode<float>(dfcaltruthshowerTag,11);
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
    }
  }


  //------------------------------------------------------------------------------


  void addDBCALTruthShower(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DBCALTruthShower*> bcaltruthshowers; 
    eventLoop->Get(bcaltruthshowers); 
    if(bcaltruthshowers.size()<=0)return;


    // create bcaltruthshower bank and add to event tree
    evioDOMNodeP bcaltruthshower = evioDOMNode::createEvioDOMNode(dbcaltruthshowerTag,0);
    tree << bcaltruthshower;


    // create data banks and add to bcaltruthshower
    evioDOMNodeP trackBank    =  evioDOMNode::createEvioDOMNode<int>  (dbcaltruthshowerTag,1);
    evioDOMNodeP primaryBank  =  evioDOMNode::createEvioDOMNode<int>  (dbcaltruthshowerTag,2);
    evioDOMNodeP phiBank      =  evioDOMNode::createEvioDOMNode<float>(dbcaltruthshowerTag,3);
    evioDOMNodeP rBank        =  evioDOMNode::createEvioDOMNode<float>(dbcaltruthshowerTag,4);
    evioDOMNodeP zBank        =  evioDOMNode::createEvioDOMNode<float>(dbcaltruthshowerTag,5);
    evioDOMNodeP tBank        =  evioDOMNode::createEvioDOMNode<float>(dbcaltruthshowerTag,6);
    evioDOMNodeP EBank        =  evioDOMNode::createEvioDOMNode<float>(dbcaltruthshowerTag,7);
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
    }
  }


//------------------------------------------------------------------------------


  void addDCDCHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DCDCHit*> cdchits; 
    eventLoop->Get(cdchits); 
    if(cdchits.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP cdchit = evioDOMNode::createEvioDOMNode(dcdchitTag,0);
    tree << cdchit;


    // create data banks and add to cdchit bank
    evioDOMNodeP ringBank  = evioDOMNode::createEvioDOMNode<int>  (dcdchitTag,1);
    evioDOMNodeP strawBank = evioDOMNode::createEvioDOMNode<int>  (dcdchitTag,2);
    evioDOMNodeP dEBank    = evioDOMNode::createEvioDOMNode<float>(dcdchitTag,3);
    evioDOMNodeP tBank     = evioDOMNode::createEvioDOMNode<float>(dcdchitTag,4);
    *cdchit << ringBank << strawBank << dEBank << tBank;


    // add track data to banks
    for(unsigned int i=0; i<cdchits.size(); i++) {
      *ringBank  << cdchits[i]->ring;
      *strawBank << cdchits[i]->straw;
      *dEBank    << cdchits[i]->dE;
      *tBank     << cdchits[i]->t;
    }
  }


//------------------------------------------------------------------------------


  void addDMCTrajectoryPoint(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DMCTrajectoryPoint*> mctrajectorypoints; 
    eventLoop->Get(mctrajectorypoints); 
    if(mctrajectorypoints.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP mctrajectorypoint = evioDOMNode::createEvioDOMNode(dmctrajectorypointTag,0);
    tree << mctrajectorypoint;


    // create data banks and add to cdchit bank
    evioDOMNodeP xBank              = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,1);
    evioDOMNodeP yBank              = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,2);
    evioDOMNodeP zBank              = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,3);
    evioDOMNodeP pxBank             = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,4);
    evioDOMNodeP pyBank             = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,5);
    evioDOMNodeP pzBank             = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,6);
    evioDOMNodeP EBank              = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,7);
    evioDOMNodeP dEBank             = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,8);
    evioDOMNodeP primary_trackBank  = evioDOMNode::createEvioDOMNode<int>   (dmctrajectorypointTag,9);
    evioDOMNodeP trackBank          = evioDOMNode::createEvioDOMNode<int>   (dmctrajectorypointTag,10);
    evioDOMNodeP partBank           = evioDOMNode::createEvioDOMNode<int>   (dmctrajectorypointTag,11);
    evioDOMNodeP radlenBank         = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,12);
    evioDOMNodeP stepBank           = evioDOMNode::createEvioDOMNode<float> (dmctrajectorypointTag,13);
    evioDOMNodeP mechBank           = evioDOMNode::createEvioDOMNode<int>   (dmctrajectorypointTag,14);
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
    }
  }


//------------------------------------------------------------------------------


  void addDFDCHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DFDCHit*> fdchits; 
    eventLoop->Get(fdchits); 
    if(fdchits.size()<=0)return;


    // create fdchit bank and add to event tree
    evioDOMNodeP fdchit = evioDOMNode::createEvioDOMNode(dfdchitTag,0);
    tree << fdchit;


    // create data banks and add to fdchit bank
    evioDOMNodeP layerBank   = evioDOMNode::createEvioDOMNode<int>  (dfdchitTag,1);
    evioDOMNodeP moduleBank  = evioDOMNode::createEvioDOMNode<int>  (dfdchitTag,2);
    evioDOMNodeP elementBank = evioDOMNode::createEvioDOMNode<int>  (dfdchitTag,3);
    evioDOMNodeP planeBank   = evioDOMNode::createEvioDOMNode<int>  (dfdchitTag,4);
    evioDOMNodeP gPlaneBank  = evioDOMNode::createEvioDOMNode<int>  (dfdchitTag,5);
    evioDOMNodeP gLayerBank  = evioDOMNode::createEvioDOMNode<int>  (dfdchitTag,6);
    evioDOMNodeP qBank       = evioDOMNode::createEvioDOMNode<float>(dfdchitTag,7);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float>(dfdchitTag,8);
    evioDOMNodeP rBank       = evioDOMNode::createEvioDOMNode<float>(dfdchitTag,9);
    evioDOMNodeP typeBank    = evioDOMNode::createEvioDOMNode<int>  (dfdchitTag,10);
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
    }
  }


//----------------------------------------------------------------------------


  void addDBeamPhoton(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DBeamPhoton*> beamphotons; 
    eventLoop->Get(beamphotons); 
    if(beamphotons.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP beamphoton = evioDOMNode::createEvioDOMNode(dbeamphotonTag,0);
    tree << beamphoton;


    // create data banks and add to cdchit bank
    evioDOMNodeP xBank   = evioDOMNode::createEvioDOMNode<float>  (dbeamphotonTag,1);
    evioDOMNodeP yBank   = evioDOMNode::createEvioDOMNode<float>  (dbeamphotonTag,2);
    evioDOMNodeP zBank   = evioDOMNode::createEvioDOMNode<float>  (dbeamphotonTag,3);
    evioDOMNodeP pxBank  = evioDOMNode::createEvioDOMNode<float>  (dbeamphotonTag,4);
    evioDOMNodeP pyBank  = evioDOMNode::createEvioDOMNode<float>  (dbeamphotonTag,5);
    evioDOMNodeP pzBank  = evioDOMNode::createEvioDOMNode<float>  (dbeamphotonTag,6);
    evioDOMNodeP tBank   = evioDOMNode::createEvioDOMNode<int>    (dbeamphotonTag,7);
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
    }
  }


//------------------------------------------------------------------------------


  void addDSCTruthHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DSCTruthHit*> sctruthhits;
    eventLoop->Get(sctruthhits); 
    if(sctruthhits.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP sctruthhit = evioDOMNode::createEvioDOMNode(dsctruthhitTag,0);
    tree << sctruthhit;


    // create data banks and add to cdchit bank
    evioDOMNodeP dEdxBank     = evioDOMNode::createEvioDOMNode<float>  (dsctruthhitTag,1);
    evioDOMNodeP primaryBank  = evioDOMNode::createEvioDOMNode<int8_t> (dsctruthhitTag,2);
    evioDOMNodeP trackBank    = evioDOMNode::createEvioDOMNode<int>    (dsctruthhitTag,3);
    evioDOMNodeP ptypeBank    = evioDOMNode::createEvioDOMNode<int>    (dsctruthhitTag,4);
    evioDOMNodeP rBank        = evioDOMNode::createEvioDOMNode<float>  (dsctruthhitTag,5);
    evioDOMNodeP phiBank      = evioDOMNode::createEvioDOMNode<float>  (dsctruthhitTag,6);
    evioDOMNodeP zBank        = evioDOMNode::createEvioDOMNode<float>  (dsctruthhitTag,7);
    evioDOMNodeP tBank        = evioDOMNode::createEvioDOMNode<float>  (dsctruthhitTag,8);
    evioDOMNodeP sectorBank   = evioDOMNode::createEvioDOMNode<int>    (dsctruthhitTag,9);
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
    }
  }


//------------------------------------------------------------------------------


  void addDFCALHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DFCALHit*> fcalhits;
    eventLoop->Get(fcalhits); 
    if(fcalhits.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP fcalhit = evioDOMNode::createEvioDOMNode(dfcalhitTag,0);
    tree << fcalhit;


    // create data banks and add to cdchit bank
    evioDOMNodeP rowBank     = evioDOMNode::createEvioDOMNode<int>   (dfcalhitTag,1);
    evioDOMNodeP columnBank  = evioDOMNode::createEvioDOMNode<int>   (dfcalhitTag,2);
    evioDOMNodeP xBank       = evioDOMNode::createEvioDOMNode<float> (dfcalhitTag,3);
    evioDOMNodeP yBank       = evioDOMNode::createEvioDOMNode<float> (dfcalhitTag,4);
    evioDOMNodeP EBank       = evioDOMNode::createEvioDOMNode<float> (dfcalhitTag,5);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float> (dfcalhitTag,6);
    *fcalhit << rowBank << columnBank << xBank << yBank << EBank << tBank;


    // add track data to banks
    for(unsigned int i=0; i<fcalhits.size(); i++) {
      *rowBank     << fcalhits[i]->row;
      *columnBank  << fcalhits[i]->column;
      *xBank       << fcalhits[i]->x;
      *yBank       << fcalhits[i]->y;
      *EBank       << fcalhits[i]->E;
      *tBank       << fcalhits[i]->t;
    }
  }


//------------------------------------------------------------------------------


  void addDHDDMBCALHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DHDDMBCALHit*> hddmbcalhits;
    eventLoop->Get(hddmbcalhits); 
    if(hddmbcalhits.size()<=0)return;


    // create hddm bcal bank and add to event tree
    evioDOMNodeP hddmbcalhit = evioDOMNode::createEvioDOMNode(dhddmbcalhitTag,0);
    tree << hddmbcalhit;


    // create data banks and add to cdchit bank
    evioDOMNodeP moduleBank  = evioDOMNode::createEvioDOMNode<int>   (dhddmbcalhitTag,1);
    evioDOMNodeP layerBank   = evioDOMNode::createEvioDOMNode<int>   (dhddmbcalhitTag,2);
    evioDOMNodeP sectorBank  = evioDOMNode::createEvioDOMNode<int>   (dhddmbcalhitTag,3);
    evioDOMNodeP EBank       = evioDOMNode::createEvioDOMNode<float> (dhddmbcalhitTag,4);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float> (dhddmbcalhitTag,5);
    evioDOMNodeP zLocalBank  = evioDOMNode::createEvioDOMNode<float> (dhddmbcalhitTag,6);
    *hddmbcalhit << moduleBank << layerBank << sectorBank << EBank << tBank << zLocalBank;


    // add track data to banks
    for(unsigned int i=0; i<hddmbcalhits.size(); i++) {
      *moduleBank  << hddmbcalhits[i]->module;
      *layerBank   << hddmbcalhits[i]->layer;
      *sectorBank  << hddmbcalhits[i]->sector;
      *EBank       << hddmbcalhits[i]->E;
      *tBank       << hddmbcalhits[i]->t;
      *zLocalBank  << hddmbcalhits[i]->zLocal;
    }
  }


//------------------------------------------------------------------------------


  void addDHDDMTOFHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DHDDMTOFHit*> hddmtofhits;
    eventLoop->Get(hddmtofhits); 
    if(hddmtofhits.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP hddmtofhit = evioDOMNode::createEvioDOMNode(dhddmtofhitTag,0);
    tree << hddmtofhit;


    // create data banks and add to cdchit bank
    evioDOMNodeP planeBank    = evioDOMNode::createEvioDOMNode<int>   (dhddmtofhitTag,1);
    evioDOMNodeP barBank      = evioDOMNode::createEvioDOMNode<int>   (dhddmtofhitTag,2);
    evioDOMNodeP ptypeBank    = evioDOMNode::createEvioDOMNode<int>   (dhddmtofhitTag,3);
    evioDOMNodeP t_northBank  = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,4);
    evioDOMNodeP dE_northBank = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,5);
    evioDOMNodeP t_southBank  = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,6);
    evioDOMNodeP dE_southBank = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,7);
    evioDOMNodeP xBank        = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,8);
    evioDOMNodeP yBank        = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,9);
    evioDOMNodeP zBank        = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,10);
    evioDOMNodeP pxBank       = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,11);
    evioDOMNodeP pyBank       = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,12);
    evioDOMNodeP pzBank       = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,13);
    evioDOMNodeP EBank        = evioDOMNode::createEvioDOMNode<float> (dhddmtofhitTag,14);
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
    }
  }


//------------------------------------------------------------------------------


  void addDSCHit(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DSCHit*> schits;
    eventLoop->Get(schits); 
    if(schits.size()<=0)return;


    // create cdchit bank and add to event tree
    evioDOMNodeP schit = evioDOMNode::createEvioDOMNode(dschitTag,0);
    tree << schit;


    // create data banks and add to cdchit bank
    evioDOMNodeP dEBank      = evioDOMNode::createEvioDOMNode<float>  (dschitTag,1);
    evioDOMNodeP tBank       = evioDOMNode::createEvioDOMNode<float>  (dschitTag,2);
    evioDOMNodeP sectorBank  = evioDOMNode::createEvioDOMNode<int>    (dschitTag,3);
    *schit << dEBank << tBank << sectorBank;


    // add track data to banks
    for(unsigned int i=0; i<schits.size(); i++) {
      *dEBank      << schits[i]->dE;
      *tBank       << schits[i]->t;
      *sectorBank  << schits[i]->sector;
    }
  }


//------------------------------------------------------------------------------


  void addDTrackTimeBased(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DTrackTimeBased*> timebasedtracks;
    eventLoop->Get(timebasedtracks); 
    if(timebasedtracks.size()<=0)return;


    // create timebasedtrack bank and add to event tree
    evioDOMNodeP timebasedtrack = evioDOMNode::createEvioDOMNode(dtracktimebasedTag,0);
    tree << timebasedtrack;


    // create data banks and add to cdchit bank
    evioDOMNodeP chisq = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,1);
    evioDOMNodeP Ndof  = evioDOMNode::createEvioDOMNode<int>    (dtracktimebasedTag,2);
    evioDOMNodeP FOM   = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,3);
    evioDOMNodeP x     = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,4);
    evioDOMNodeP y     = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,5);
    evioDOMNodeP z     = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,6);
    evioDOMNodeP px    = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,7);
    evioDOMNodeP py    = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,8);
    evioDOMNodeP pz    = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,9);
    evioDOMNodeP q     = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,10);
    evioDOMNodeP E     = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,11);
    evioDOMNodeP mass  = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,12);
    evioDOMNodeP t0    = evioDOMNode::createEvioDOMNode<float>  (dtracktimebasedTag,13);
    *timebasedtrack << chisq << Ndof << FOM << x << y << z << px << py << pz 
                    << q << E << mass << t0;


    // add track data to banks
    // also create map if DChargedTrack is requested
    if(evioMap["dtracktimebased"])DTrackTimeBasedMap.clear();
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
      DTrackTimeBasedMap[timebasedtracks[i]->id]=i;
    }
  }


//------------------------------------------------------------------------------


  void addDChargedTrack(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DChargedTrack*> chargedtracks;
    eventLoop->Get(chargedtracks); 
    if(chargedtracks.size()<=0)return;


    // create chargedtrack bank and add to event tree
    evioDOMNodeP chargedtrack = evioDOMNode::createEvioDOMNode(dchargedtrackTag,0);
    tree << chargedtrack;

    // create index bank for each charged track and add to chargedtrack bank
    for(unsigned int i=0; i<chargedtracks.size(); i++) {
      evioDOMNodeP hypotheses = evioDOMNode::createEvioDOMNode<int> (dchargedtrackTag,1);
      *chargedtrack << hypotheses;
      for(unsigned int j=0; j<chargedtracks[i]->hypotheses.size(); j++) {
        *hypotheses<< DTrackTimeBasedMap[chargedtracks[i]->hypotheses[j]->id];
      }
    }
  }


//------------------------------------------------------------------------------


  void addDPhoton(JEventLoop *eventLoop, evioDOMTree &tree) {


    // is there any data
    vector<const DPhoton*> photons;
    eventLoop->Get(photons); 
    if(photons.size()<=0)return;


    // create photon bank and add to event tree
    evioDOMNodeP photon = evioDOMNode::createEvioDOMNode(dphotonTag,0);
    tree << photon;


    // create data banks and add to photon bank
    evioDOMNodeP E     = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,1);
    evioDOMNodeP px    = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,2);
    evioDOMNodeP py    = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,3);
    evioDOMNodeP pz    = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,4);
    evioDOMNodeP x     = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,5);
    evioDOMNodeP y     = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,6);
    evioDOMNodeP z     = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,7);
    evioDOMNodeP t     = evioDOMNode::createEvioDOMNode<float>  (dphotonTag,8);
    evioDOMNodeP Tag   = evioDOMNode::createEvioDOMNode<int>    (dphotonTag,9);
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
    }
  }


//------------------------------------------------------------------------------


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
          evioMap["dmcthrown"]=!minus;
          evioMap["dtoftruth"]=!minus;
          evioMap["dfcaltruthshower"]=!minus;
          evioMap["dbcaltruthshower"]=!minus;

        } else if(value=="hits") {
          evioMap["dcdchit"]=!minus;
          evioMap["dfdchit"]=!minus;
          evioMap["dmctrackhit"]=!minus;

        } else if(value=="tracks") {
          evioMap["dtracktimebased"]=!minus;
          evioMap["dchargedtrack"]=!minus;
          evioMap["dphoton"]=!minus;

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
