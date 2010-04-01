// $Id$

//    File: DDANAEVIO_factory.h
// Created: Mon Mar 15 09:08:37 EDT 2010
// Creator: wolin (on Linux stan.jlab.org 2.6.18-164.el5 x86_64)



#ifndef _DDANAEVIO_factory_
#define _DDANAEVIO_factory_


#include <string>
#include <map>
#include <set>
using namespace std;


#include <JANA/JFactory.h>
using namespace jana;


#include <DDANAEVIODOMTree.h>


// list of factory tags for dana objects that to be added to tree by default
// use -PEVIO:DANAEVIO to override
// *** NOTE:  if you add to this list be sure to modify decode_object_parameters() appropriately ***
static set<string> emptySet;
static string untagged[] = {string("")};
static set<string> untaggedSet(untagged,untagged+1);
static pair< string, set<string> > danaObs[] =  {
  pair< string, set<string> > ("dmctrackhit",          emptySet),
  pair< string, set<string> > ("dbeamphoton",          untaggedSet),
  pair< string, set<string> > ("dmcthrown",            untaggedSet),
  pair< string, set<string> > ("dfcaltruthshower",     untaggedSet),
  pair< string, set<string> > ("dbcaltruthshower",     untaggedSet),
  pair< string, set<string> > ("dtoftruth",            untaggedSet),
  pair< string, set<string> > ("dsctruthhit",          untaggedSet),
  pair< string, set<string> > ("dmctrajectorypoint",   emptySet),
  pair< string, set<string> > ("dcdchit",              untaggedSet),
  pair< string, set<string> > ("dfdchit",              untaggedSet),
  pair< string, set<string> > ("dfcalhit",             untaggedSet),
  pair< string, set<string> > ("dhddmbcalhit",         untaggedSet),
  pair< string, set<string> > ("dhddmtofhit",          untaggedSet),
  pair< string, set<string> > ("dschit",               untaggedSet),
  pair< string, set<string> > ("dtrackwirebased",      emptySet),
  pair< string, set<string> > ("dtracktimebased",      emptySet),
  pair< string, set<string> > ("dchargedtrack",        emptySet),
  pair< string, set<string> > ("dphoton",              emptySet),
  pair< string, set<string> > ("dcdctrackhit",         emptySet),
  pair< string, set<string> > ("dfdcpseudo",           emptySet),
};


// global map of which factory/tags to convert
static map<string, set<string> > evioMap(danaObs,danaObs+sizeof(danaObs)/sizeof(danaObs[0]));


// holds tag/num pairs for all DANA objects
static map< string, pair<int,int> > tagMap;



//------------------------------------------------------------------------------------


class DDANAEVIO_factory : public JFactory<DDANAEVIODOMTree> {
  
 public:
  DDANAEVIO_factory();
  ~DDANAEVIO_factory() {};
  
  
 private:
  jerror_t evnt(JEventLoop *eventLoop, int eventnumber);
  void decode_DANAEVIO_parameter(void);
  void get_tagNum_dictionary(void);
  static void startElement(void *userData, const char *xmlname, const char **atts);


  void addObjIdBank(evioDOMTree &tree);

  void addDMCTrackHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDBeamPhoton(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDMCThrown(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDFCALTruthShower(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDBCALTruthShower(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDTOFTruth(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDSCTruthHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDMCTrajectoryPoint(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDCDCHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDFDCHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDFCALHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDHDDMBCALHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDHDDMTOFHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDSCHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDCDCTrackHit(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDFDCPseudo(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDTrackWireBased(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDTrackTimeBased(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDChargedTrack(JEventLoop *eventLoop, evioDOMTree &tree);
  void addDPhoton(JEventLoop *eventLoop, evioDOMTree &tree);


  // event-specific global object id map
  map<int,string> objIdMap;



  // templated methods must be in header file
  template<typename T> evioDOMNodeP createLeafNode(string nameId) {
    pair<int,int> p = tagMap[nameId];;
    return(evioDOMNode::createEvioDOMNode<T>(p.first,p.second));
  }

  // might as well put this here...
  evioDOMNodeP createContainerNode(string nameId) {
    pair<int,int> p = tagMap[nameId];;
    return(evioDOMNode::createEvioDOMNode(p.first,p.second));
  }
  
};

#endif // _DDANAEVIO_factory_


//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
