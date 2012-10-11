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
//    Can change via -PRAWEVENT:FILEBASE=newBaseName
//
//  Reads translation table in via -PRAWEVENT:TRANSLATION=fileName.xml
//    default is fakeTranslationTable.xml
//
//
// To run:
//
//     $ hd_ana -PPLUGINS=rawevent inputFile.hddm
//
//
// still to do:
//    pair spectrometer
//
//
//  ejw, 27-jun-2011



#include "rawevent/JEventProcessor_rawevent.h"
#include "mc2coda.h"


#include<sstream>
#include<iomanip>
#include<algorithm>
#include <expat.h>

#include <boost/lexical_cast.hpp>
using namespace boost;


// to protect writing to output file
static pthread_mutex_t rawMutex = PTHREAD_MUTEX_INITIALIZER;


// for evio output and translation table
static evioFileChannel *chan       = NULL;
static string fileBase             = "rawevent";
static string outputFileName;
static string translationTableName = "fakeTranslationTable.xml";


// current run number
static int runNumber;


// csc map converts from detector spec to (crate,slot,channel)
//  key is detector-dependent encoded string (e.g. "cdcadc::2:25" for CDC ring 2 straw 25)
//  value is struct containing (crate,slot,channel)
static map<string,cscVal> cscMap;


// detector map (inverse of csc map) is 3-dimensional array of strings with indices (crate,slot,channel)
//  content is detector-dependent encoded string
#define MAX_CRATE   58+1
#define MAX_SLOT    16+1
#define MAX_CHANNEL 72+1
static string detectorMap[MAX_CRATE][MAX_SLOT][MAX_CHANNEL];


// for Dave's mc2coda package
#define MAXCRATE 100
#define MAXSLOT  21
#define MAXSIZE  50000
static string expName = "HallD";
static CODA_EXP_INFO *codaExpInfo = NULL;
static int expid    = 1;

static int nCrate   = 0;
static int crateID[MAXCRATE];
static int nModules[MAXCRATE];
static int modules[MAXCRATE][MAXSLOT];
static int detID[MAXCRATE][MAXSLOT];




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
  gPARMS->SetDefaultParameter("RAWEVENT:FILEBASE",fileBase);

  // get translation table file name from command line 
  gPARMS->SetDefaultParameter("RAWEVENT:TRANSLATION",translationTableName);

  // read translation table
  readTranslationTable();

  // initialize mc2coda package
  codaExpInfo=mc2codaInitExp(nCrate,expName.c_str());

  // feed crate-specific info to mc2coda package
  int stat;
  for(int i=0; i<nCrate; i++) {
    stat = mc2codaSetCrate(codaExpInfo,crateID[i],nModules[i],modules[i],detID[i]);
    if(stat!=0) {
      cerr << "?JEventProcessor_rawevent...error return from mc2codaSetCrate()" << endl << endl;
      exit(EXIT_FAILURE);
    }
  }

}


//----------------------------------------------------------------------------


// ~JEventProcessor_rawevent (Destructor) once only
JEventProcessor_rawevent::~JEventProcessor_rawevent() {
  mc2codaFree(expID);
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
  chan = new evioFileChannel(outputFileName,"w");
  chan->open();
  jout << endl << "   opening output file:   " << outputFileName << endl << endl << endl;


  // add header event if required
  // ...


  return NOERROR;
}


//----------------------------------------------------------------------------


extern "C" {
  bool compareDTOFRawHits(const DTOFRawHit* h1, const DTOFRawHit* h2) {
    if(h1->plane!=h2->plane) {
      return(h1->plane<h2->plane);
    } else if(h1->bar!=h2->bar) {
      return(h1->bar<h2->bar);
    } else if(h1->lr!=h2->lr) {
      return(h1->lr<h2->lr);
    } else {
      return(h1->t<h2->t);
    }
  }
}


//----------------------------------------------------------------------------


// Called once per event in many different processing threads, so:
//
//
//    *** MUST be thread-safe ***
//
//
// Note:  Accesses MC hit data in DANA objects and converts to crate/slot/channel/energy/time.
//        DAQ group task is to take this data, sort and combine the hits, then create an EVIO
//          buffer in the format planned for disentangled events.
//        The buffer is written to disk using a mutex-locked EVIO write.

jerror_t JEventProcessor_rawevent::evnt(JEventLoop *eventLoop, int eventnumber) {

  unsigned int i;
  CODA_HIT_INFO hit;
  int hData[256];


  // initialize mc2coda event buffer info
  int nhit                 = 0;
  uint64_t trigTime        = 0;
  unsigned short eventType = 0;
  unsigned int *eventPtr = mc2codaOpenEvent(expID, (uint64_t)eventnumber, trigTime, eventType, MAXSIZE);


  // get all types of raw hit data, sort according to id and time order, then feed to mc2coda

  // DTOFRawHit - FADC250 and F1TDC32 (60 ps)
  vector<const DTOFRawHit*> dtofrawhits; 
  eventLoop->Get(dtofrawhits);
  sort(dtofrawhits.begin(),dtofrawhits.end(),compareDTOFRawHits);

  for(i=0; i<dtofrawhits.size(); i++) {
    nhit++;
    float dE  = dtofrawhits[i]->dE;
    float t   = dtofrawhits[i]->t;

    // feed hit and crate/slot/channel info to mc2coda package
    cscRef cscADC   = DTOFRawHitTranslationADC(dtofrawhits[i]);
    hit.hit_id      = nhit;
    hit.det_id      = ;
    hit.crate_id    = cscADC.crate;
    hit.slot_id     = cscADC.slot;
    hit.module_id   = ;
    hit.chan_id     = cscADC.channel;
    hit.nsamples    = 1;
    hit.time_offset = ;

    hData[0]        = int(dE*1000.);
    hit.hdata       = hData;
    mc2codaWrite(eventID,1,&hit);

    cscRef cscTDC  = DTOFRawHitTranslationTDC(dtofrawhits[i]);
  }


  // DBCALHit - FADC250 and F1TDC32 (60 ps)
  vector<const DBCALHit*> dbcalhits;
  eventLoop->Get(dbcalhits);
  for(i=0; i<dbcalhits.size(); i++) {
    float E      = dbcalhits[i]->E;
    float t      = dbcalhits[i]->t;

    cscRef cscADC  = DBCALHitTranslationADC(dbcalhits[i]);
    cscRef cscTDC  = DBCALHitTranslationTDC(dbcalhits[i]);

    // feed data to DAQ function...
  }      



  // DFCALHit - FADC250
  vector<const DFCALHit*> dfcalhits;
  eventLoop->Get(dfcalhits);
  for(i=0; i<dfcalhits.size(); i++) {
    float E      = dfcalhits[i]->E;
    float t      = dfcalhits[i]->t;
    
    cscRef cscADC  = DFCALHitTranslationADC(dfcalhits[i]);

    // feed data to DAQ function...
  }      


  // DFDCHit - cathode strips FADC125 or anode wires F1TDC48 (115 ps)
  vector<const DFDCHit*> dfdchits; 
  eventLoop->Get(dfdchits);
  for(unsigned int i=0; i<dfdchits.size(); i++) {
    float q      = dfdchits[i]->q;
    float t      = dfdchits[i]->t;

    int type = dfdchits[i]->type;
    if(type==0) {           // F1TDC48
      cscRef cscTDC  = DFDCAnodeHitTranslation(dfdchits[i]);

      // feed data to DAQ function...

    } else if(type==1) {    // FADC125
      cscRef cscADC  = DFDCCathodeHitTranslation(dfdchits[i]);

      // feed data to DAQ function...
    }
  }


  // DCDCHit - FADC125
  vector<const DCDCHit*> dcdchits; 
  eventLoop->Get(dcdchits);
  for(i=0; i<dcdchits.size(); i++) {
    float dE     = dcdchits[i]->dE;
    float t      = dcdchits[i]->t;

    cscRef cscADC  = DCDCHitTranslationADC(dcdchits[i]);

    // feed data to DAQ function...
  }      


  // DSCHit - FADC250 and F1TDC32 (60 ps)
  vector<const DSCHit*> dschits;
  eventLoop->Get(dschits);
  for(unsigned int i=0; i<dschits.size(); i++) {
    float dE     = dschits[i]->dE;
    float t      = dschits[i]->t;

    cscRef cscADC  = DSCHitTranslationADC(dschits[i]);
    cscRef cscTDC  = DSCHitTranslationTDC(dschits[i]);

    // feed data to DAQ function...
  }      


  // DTagger - FADC250 and F1TDC32 (60 ps)
  vector<const DTagger*> dtaggerhits;
  eventLoop->Get(dtaggerhits);
  for(i=0; i<dtaggerhits.size(); i++) {
    float E      = dtaggerhits[i]->E;
    float t      = dtaggerhits[i]->t;

    cscRef cscADC  = DTaggerTranslationADC(dtaggerhits[i]);
    cscRef cscTDC  = DTaggerTranslationTDC(dtaggerhits[i]);

    // feed data to DAQ function...
  }      



  // close evio buffer, write it out, then free buffer
  uint32_t *buffer = mc2codaCloseEvent(eventPtr, (void*)eventPtr);
  pthread_mutex_lock(&rawMutex);
  chan->write(buffer);
  pthread_mutex_unlock(&rawMutex);
  mc2codaFreeEvent(eventPtr);

  return NOERROR;
}


//----------------------------------------------------------------------------


// erun called once-only at end of run, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::erun(void) {

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
//  The following routines access the translation table
//----------------------------------------------------------------------------
//----------------------------------------------------------------------------


void JEventProcessor_rawevent::readTranslationTable(void) {

  // create parser and specify element handlers
  XML_Parser xmlParser = XML_ParserCreate(NULL);
  if(xmlParser==NULL) {
    jerr << endl << endl << "readTranslationTable...unable to create parser" << endl << endl;
    exit(EXIT_FAILURE);
  }
  XML_SetElementHandler(xmlParser,startElement,endElement);


  // clear crate summary info needed by Dave's library


  // open and parse the file
  FILE *f = fopen(translationTableName.c_str(),"r");
  if(f!=NULL) {
    int status,len;
    bool done;
    const int bufSize = 50000;
    char *buf = new char[bufSize];
    do {
      len  = fread(buf,1,bufSize,f);
      done = len!=bufSize;
      status=XML_Parse(xmlParser,buf,len,done);
      if((!done)&&(status==0)) {
        jerr << endl << endl << endl << "  ?readTranslationTable...parseXMLFile parse error for " << translationTableName
             << endl << endl << XML_ErrorString(XML_GetErrorCode(xmlParser))
             << endl << endl << endl;
        fclose(f);
        delete [] buf;
        XML_ParserFree(xmlParser);
        return;
      }
    } while (!done);
    fclose(f);
    delete [] buf;
    jout << endl << endl << " Successfully read translation table:  " << translationTableName << endl << endl << endl;

  } else {
    jerr << endl << endl << endl << "  ?readTranslationTable...unable to open " << translationTableName
         << endl << endl << endl;
  }
  XML_ParserFree(xmlParser);
}


//----------------------------------------------------------------------------


void JEventProcessor_rawevent::startElement(void *userData, const char *xmlname, const char **atts) {
  
  static int crate=0, slot=0;
  static string Type="Unknown",type="unknown";

  int channel = 0;
  string Detector;
  string end;
  string row,column,module,sector,layer,chan;
  string ring,straw,plane,bar,gPlane,element;


  // store crate summary info, fill both maps
  if(strcasecmp(xmlname,"translation_table")==0) {
    // do nothing


  } else if(strcasecmp(xmlname,"crate")==0) {
    for (int i=0; atts[i]; i+=2) {
      if(strcasecmp(atts[i],"number")==0) {
        crate = atoi(atts[i+1]);
        break;
      }
    }
    nCrate++;
    crateID[nCrate-1]=crate;
    nModules[nCrate-1]=0;


  } else if(strcasecmp(xmlname,"slot")==0) {
    for (int i=0; atts[i]; i+=2) {
      if(strcasecmp(atts[i],"number")==0) {
        slot = atoi(atts[i+1]);
      } else if(strcasecmp(atts[i],"type")==0) {
        Type = string(atts[i+1]);
        std::transform(Type.begin(), Type.end(), type.begin(), (int(*)(int)) tolower);
      }
    }
    nModules[nCrate-1]++;
    modules[nCrate-1][nModules[nCrate-1]-1] = slot;
    detID[nCrate-1][nModules[nCrate-1]-1]   = type2detID(type);


  } else if(strcasecmp(xmlname,"channel")==0) {

    for (int i=0; atts[i]; i+=2) {
      if(strcasecmp(atts[i],"number")==0) {
        channel = atoi(atts[i+1]);
      } else if(strcasecmp(atts[i],"detector")==0) {
        Detector = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"row")==0) {
        row = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"column")==0) {
        column = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"module")==0) {
        module = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"sector")==0) {
        sector = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"layer")==0) {
        layer = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"end")==0) {
        end = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"chan")==0) {
        chan = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"ring")==0) {
        ring = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"straw")==0) {
        straw = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"gPlane")==0) {
        gPlane = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"element")==0) {
        element = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"plane")==0) {
        plane = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"bar")==0) {
        bar = string(atts[i+1]);
      }
    }


    // fill maps
    cscVal csc = {crate,slot,channel};
    string detector = Detector;
    std::transform(detector.begin(), detector.end(), detector.begin(), (int(*)(int)) tolower);

    string s="unknown::";
    if(detector=="fcal") {
      if(type=="fadc250") {
        s = "fcaladc::";
      } else {
        s = "unknownFCAL::";
        jerr << endl << endl << "?startElement...illegal type for FCAL: " << Type << endl << endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;


    } else if(detector=="bcal") {
      if(type=="f1tdc32") {
        s = "bcaltdc::";
      } else if (type=="fadc250") {
        s = "bcaladc::";
      } else {
        s = "unknownBCAL::";
        jerr << endl << endl << "?startElement...illegal type for BCAL: " << Type << endl << endl;
      }
      s += module + ":" + sector + ":" + layer + ":" + end;
      cscMap[s] = csc;
      

    } else if(detector=="cdc") {
      if(type=="fadc125") {
        s = "cdcadc::";
      } else {
        s = "unknownCDC::";
        jerr << endl << endl << "?startElement...illegal type for CDC: " << Type << endl << endl;
      }
      s += ring + ":" + straw;
      cscMap[s] = csc;
        
      
    } else if(detector=="sc") {
      if(type=="f1tdc32") {
        s = "sctdc::";
      } else if (type=="fadc250") {
        s = "scadc::";
      } else {
        s = "unknownSC::";
        jerr << endl << endl << "?startElement...illegal type for SC: " << Type << endl << endl;
      }
      s += sector;
      cscMap[s] = csc;
    

    } else if(detector=="fdccathode") {
      if(type=="fadc125") {
        s = "fdccathode::";
      } else {
        s = "unknownFDCCathode::";
        jerr << endl << endl << "?startElement...illegal type for FDC Cathode: " << Type << endl << endl;
      }
      s += gPlane + ":" + element;
      cscMap[s] = csc;

      
    } else if(detector=="fdcanode") {
      if(type=="f1tdc48") {
        s = "fdcanode::";
      } else {
        s = "unknownFDCAnode::";
        jerr << endl << endl << "?startElement...illegal type for FDC Anode: " << Type << endl << endl;
      }
      s += gPlane + ":" + element;
      cscMap[s] = csc;

      
    } else if(detector=="tof") {
      if(type=="f1tdc32") {
        s = "toftdc::";
      } else if (type=="fadc250") {
        s = "tofadc::";
      } else {
        s = "unknownTOF::";
        jerr << endl << endl << "?startElement...illegal type for TOF: " << Type << endl << endl;
      }
      s += plane + ":" + bar + ":" + end;
      cscMap[s] = csc;

     
    } else if(detector=="tagger") {
      if(type=="f1tdc32") {
        s = "taggertdc::";
      } else if (type=="fadc250") {
        s = "taggeradc::";
      } else {
        s = "unknownTagger::";
        jerr << endl << endl << "?startElement...illegal type for TAGGER: " << Type << endl << endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

      
    } else {
      jerr << endl << endl << "?startElement...unknown detector " << Detector << endl << endl;
    }


    // fill detector map, index is crate,slot,channel
    detectorMap[crate][slot][channel] = s;
    


  } else {
    jerr << endl << endl << "?startElement...unknown xml tag " << xmlname << endl << endl;
  }

}
    

//--------------------------------------------------------------------------


void JEventProcessor_rawevent::endElement(void *userData, const char *xmlname) {

  // finished a crate, feed crate summary info to Dave's library
  if(strcasecmp(xmlname,"crate")==0) {
  }

  // clear crate summary info
}


//--------------------------------------------------------------------------



//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//  aux routines encode hit info into string for inverse lookup table
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTOFRawHitTranslationADC(const DTOFRawHit* hit) const {
  string s = "tofadc::" + lexical_cast<string>(hit->plane) + ":" + lexical_cast<string>(hit->bar)
    + ":" + lexical_cast<string>(hit->lr);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTOFRawHitTranslationTDC(const DTOFRawHit* hit) const {
  string s = "toftdc::" + lexical_cast<string>(hit->plane) + ":" + lexical_cast<string>(hit->bar)
    + ":" + lexical_cast<string>(hit->lr);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DBCALHitTranslationADC(const DBCALHit *hit) const {
  string s = "bcaladc::" + lexical_cast<string>(hit->module) + ":" + lexical_cast<string>(hit->sector)
    + ":" + lexical_cast<string>(hit->layer) + ":" + lexical_cast<string>(hit->end);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DBCALHitTranslationTDC(const DBCALHit *hit) const {
  string s = "bcaltdc::" + lexical_cast<string>(hit->module) + ":" + lexical_cast<string>(hit->sector)
    + ":" + lexical_cast<string>(hit->layer) + ":" + lexical_cast<string>(hit->end);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFCALHitTranslationADC(const DFCALHit* hit) const {
  string s = "fcaladc::" + lexical_cast<string>(hit->row) + ":" + lexical_cast<string>(hit->column);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFDCAnodeHitTranslation(const DFDCHit* hit) const {
  string s = "fdcanode::"  + lexical_cast<string>(hit->gPlane) + ":" + lexical_cast<string>(hit->element);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFDCCathodeHitTranslation(const DFDCHit* hit) const {
  string s = "fdccathode::"  + lexical_cast<string>(hit->gPlane) + ":" + lexical_cast<string>(hit->element);

  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DCDCHitTranslationADC(const DCDCHit* hit) const {
  string s = "cdcadc::" + lexical_cast<string>(hit->ring) + ":" + lexical_cast<string>(hit->straw);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DSCHitTranslationADC(const DSCHit* hit) const {
  string s = "scadc::" + lexical_cast<string>(hit->sector);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DSCHitTranslationTDC(const DSCHit* hit) const {
  string s = "sctdc::" + lexical_cast<string>(hit->sector);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTaggerTranslationTDC(const DTagger* hit) const {
  string s = "taggertdc::" + lexical_cast<string>(hit->row) +":" + lexical_cast<string>(hit->column);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTaggerTranslationADC(const DTagger* hit) const {
  string s = "taggeradc::" + lexical_cast<string>(hit->row) +":" + lexical_cast<string>(hit->column);
  return(cscMap[s]);
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//--------------------------------------------------------------------------
