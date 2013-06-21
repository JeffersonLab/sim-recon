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
//  mc2coda expects time in natural units of the readout module: 
//       25 ps/count for CAENTDC
//       60 ps/count for F1TDC32
//      120 ps/count for F1TDC48
//        4 ns/count for FADC250
//        8 ns/count for FADC125
//
//  trig time minimum is 8 us or 2000 counts using 4 ns/count
//
//
//
//  to compile/link:
//    source some-halld-build
//    setenv EVIOROOT /group/da/ejw/coda/Linux-xxx
//    setenv HALLD_MY /home/wolin/halld_my
//
//
// To run:
//
//     $ hd_ana -PPLUGINS=rawevent -PEVENTS_TO_KEEP=1 -EVENTS_TO_SKIP=0 /local/scratch/filtered.hddm
//
//
// still to do:
//    remove ROOT
//    negative times ok
//    energy, charge units?  bombproof values for mc2coda
//    bombproofing should be done in mcsmear
//    too many hits in one channel?  hit merging?
//    pair spectrometer, others
//
//
//  ejw, 5-Nov-2012



#include "rawevent/JEventProcessor_rawevent.h"


// root histograms...probably not thread-safe
static TH1F *tofEnergies;
static TH1F *fcalEnergies;
static TH1F *bcalEnergies;
static TH1F *stEnergies;
static TH1F *tagEnergies;

static TH1F *fdcCharges;
static TH1F *cdcCharges;

static TH1F *tofTimes;
static TH1F *fdcTimes;
static TH1F *cdcTimes;
static TH1F *fcalTimes;
static TH1F *bcalTimes;
static TH1F *stTimes;
static TH1F *tagTimes;


//  is this thread-safe?
extern "C" {
#include "mc2coda.h"
}


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
#define MAXDCRATE   72+1
#define MAXDSLOT    21+1
#define MAXDCHANNEL 72+1
static string detectorMap[MAXDCRATE][MAXDSLOT][MAXDCHANNEL];


// for Dave's mc2coda package
#define MAXCRATE      128
#define MAXSLOT       21
#define MAXEVENTSIZE  30000*4    // in bytes


static string expName = "HallD";
static CODA_EXP_INFO *expID     = NULL;
static CODA_EVENT_INFO *eventID = NULL;

static int nCrate   = 0;
static int crateID[MAXCRATE];
static int nModules[MAXCRATE];
static int modules[MAXCRATE][MAXSLOT];
static int detID[MAXCRATE][MAXSLOT];


// misc
static uint64_t trigTime     = 32000000;  // in picoseconds
static float tMin            = -100000.;  // minimum hit time in picoseconds

static int trigtick          = 4000;    // in picoseconds
static int CAENTDCtick       = 25;      // in picoseconds
static int F1TDC32tick       = 60;      // in picoseconds
static int F1TDC48tick       = 120;     // in picoseconds
static int FADC250tick       = 4000;    // in picoseconds
static int FADC125tick       = 8000;    // in picoseconds


// debug
static int dumphits  = 0;
static int nomc2coda = 0;
static int noroot    = 0;



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

  // trigger time picoseconds
  gPARMS->SetDefaultParameter("RAWEVENT:TRIGTIME",trigTime);

  // minimum time in picoseconds
  gPARMS->SetDefaultParameter("RAWEVENT:TMIN",tMin);

  // option to turn off mc2coda output
  gPARMS->SetDefaultParameter("RAWEVENT:NOMC2CODA",nomc2coda);

  // option to turn off root
  gPARMS->SetDefaultParameter("RAWEVENT:NOROOT",noroot);

  // option to dump hits
  gPARMS->SetDefaultParameter("RAWEVENT:DUMPHITS",dumphits);


  // read translation table, fill crate id arrays
  readTranslationTable();


  // initialize mc2coda package
  if(nomc2coda==0) {
    expID=mc2codaInitExp(nCrate,expName.c_str());
    if(expID==NULL) {
      jerr << "?NULL return from mc2codaInitExp()" << endl;
      exit(EXIT_FAILURE);
    }
  }


  // feed crate-specific info to mc2coda package, note that VMECPU and TID not included in module count
  int stat;
  if(nomc2coda==0) {
    for(int i=0; i<nCrate; i++) {
	  int nMod = nModules[i]-2;
	  if(nMod<1) continue; // ignore crates with no digitization modules
      stat = mc2codaSetCrate(expID,crateID[i],nModules[i]-2,modules[i],detID[i]);
      if(stat==-1) {
        jerr << "?JEventProcessor_rawevent...error return from mc2codaSetCrate()" << endl << endl;
        exit(EXIT_FAILURE);
      }
    }
  }

}


//----------------------------------------------------------------------------


// ~JEventProcessor_rawevent (Destructor) once only
JEventProcessor_rawevent::~JEventProcessor_rawevent() {
  if(nomc2coda==0) {
    mc2codaFree(expID);
  }
}


//----------------------------------------------------------------------------


// init called once-only at beginning, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::init(void) {

  // root histograms
  if(noroot==0) {
    tofEnergies  = new TH1F("tofe",  "TOF energies in keV",1000,0.,5000.);
    stEnergies   = new TH1F("ste",   "ST energies in keV",1000,0.,5000.);
    fcalEnergies = new TH1F("fcale", "FCAL energies in keV",1000,0.,500000.);
    bcalEnergies = new TH1F("bcale", "BCAL energies in keV",1000,0.,1000000.);
    tagEnergies  = new TH1F("tage",  "Tagger energies in keV",1000,0.,40000.);

    fdcCharges   = new TH1F("fdcq",  "FDC charges in fC",1000,0.,1500.);
    cdcCharges   = new TH1F("cdcq",  "CDC charges in fC",1000,0.,25000.);

    tofTimes     = new TH1F("toft","  TOF times in nsec",1000,0.,500.-tMin/1000);
    stTimes      = new TH1F("stt",   "ST times in nsec",1000,0.,250.-tMin/1000);
    fdcTimes     = new TH1F("fdct",  "FDC times in nsec",1000,0.,2000.-tMin/1000);
    cdcTimes     = new TH1F("cdct",  "CDC times in nsec",1000,0.,1000.-tMin/1000);
    bcalTimes    = new TH1F("bcalt", "BCAL times in nsec",1000,0.,200.-tMin/1000);
    fcalTimes    = new TH1F("fcalt", "FCAL times in nsec",1000,0.,200.-tMin/1000);
    tagTimes     = new TH1F("tagt",  "Tagger times in nsec",1000,0.,250.-tMin/1000);
  }

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
//  map<string,cscVal>::iterator iter = cscMap.begin();
//  for(; iter!=cscMap.end(); iter++){
//  	string key = iter->first;
//	if(key.find("toftdc") == 0) _DBG_<<key<<endl;
//  }
  
  
  return NOERROR;
}


//----------------------------------------------------------------------------


  static bool compareDTOFRawHits(const DTOFRawHit* h1, const DTOFRawHit* h2) {
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

  static bool compareDFCALHits(const DFCALHit* h1, const DFCALHit* h2) {
    if(h1->row!=h2->row) {
      return(h1->row<h2->row);
    } else if(h1->column!=h2->column) {
      return(h1->column<h2->column);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDBCALHits(const DBCALHit* h1, const DBCALHit* h2) {
    if(h1->module!=h2->module) {
      return(h1->module<h2->module);
    } else if(h1->sector!=h2->sector) {
      return(h1->sector<h2->sector);
    } else if(h1->layer!=h2->layer) {
      return(h1->layer<h2->layer);
    } else if(h1->end!=h2->end) {
      return(h1->end<h2->end);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDFDCHits(const DFDCHit* h1, const DFDCHit* h2) {
    if(h1->gPlane!=h2->gPlane) {
      return(h1->gPlane<h2->gPlane);
    } else if(h1->element!=h2->element) {
      return(h1->element<h2->element);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDCDCHits(const DCDCHit* h1, const DCDCHit* h2) {
    if(h1->ring!=h2->ring) {
      return(h1->ring<h2->ring);
    } else if(h1->straw!=h2->straw) {
      return(h1->straw<h2->straw);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDSTHits(const DSCHit* h1, const DSCHit* h2) {
    if(h1->sector!=h2->sector) {
      return(h1->sector<h2->sector);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDTaggerHits(const DTagger* h1, const DTagger* h2) {
    if(h1->row!=h2->row) {
      return(h1->row<h2->row);
    } else if(h1->column!=h2->column) {
      return(h1->column<h2->column);
    } else {
      return(h1->t<h2->t);
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
  CODA_HIT_INFO hit[10];
  uint32_t mcData[10];
  int stat,nhits,hc;
  static bool first_time = true;


  // initialize event buffer info
  int hitCount             = 0;
  int detID                = 1;
  unsigned short eventType = 0;


  // open event, default max event size is 1 MB
  if(nomc2coda==0) {
    if(first_time) {
      first_time=false;
      eventID = mc2codaOpenEvent(expID, (uint64_t)eventnumber, trigTime/trigtick, eventType, MAXEVENTSIZE);
      if(eventID==NULL) {
        jerr << "?NULL return from mc2codaOpenEvent()" << endl << endl;
        exit(EXIT_FAILURE);
      }
    } else {
      int stat = mc2codaResetEvent(eventID, (uint64_t)eventnumber, trigTime/trigtick, eventType);
      if(stat!=0) {
        jerr << "?error return from mc2codaResetEvent()" << endl << endl;
        exit(EXIT_FAILURE);
      }
    }
  }


  // get all raw hit data, sort according to id and time order before feeding to mc2coda package
  // scale data appropriately to minimize loss of precision by using full available bit range


  // DCDCHit - FADC125
  vector<const DCDCHit*> dcdchits; 
  eventLoop->Get(dcdchits);
  sort(dcdchits.begin(),dcdchits.end(),compareDCDCHits);

  hc=0;
  for(i=0; i<dcdchits.size(); i++) {
    if((dcdchits[i]->q>0)&&((dcdchits[i]->t*1000.)>tMin)&&(dcdchits[i]->t*1000.<trigTime)) {

      uint32_t q     = dcdchits[i]->q;            // in femtoCoulombs
      uint32_t t     = dcdchits[i]->t*1000.-tMin;  // in picoseconds
      
      if(noroot==0)cdcCharges->Fill(dcdchits[i]->q);
      if(noroot==0)cdcTimes->Fill(dcdchits[i]->t-tMin/1000);

      hc++;
      hitCount++;
      nhits=1;
      cscRef cscADC      = DCDCHitTranslationADC(dcdchits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC125;
      hit[0].module_mode = FADC125_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = q;  // in fC
      hit[0].hdata[1]    = t/FADC125tick;
      if(q>0x7ffff)cerr << "q too large for CDC: " << q << endl;
      
      if(dumphits>1) {
        jout << endl;
        jout << " CDC ADC ring,straw are " << dcdchits[i]->ring << ", " << dcdchits[i]->straw << endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << endl;
        jout << " q,t are " << q << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for CDC ADC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }      
  if((dumphits>=1)&&(hc>0)) {
    jout << endl << "CDC hits: " << hc << endl << endl;
  }


  // DTOFRawHit - FADC250 and F1TDC32 (60 ps)
  vector<const DTOFRawHit*> dtofrawhits; 
  eventLoop->Get(dtofrawhits);
  sort(dtofrawhits.begin(),dtofrawhits.end(),compareDTOFRawHits);

  hc=0;
  for(i=0; i<dtofrawhits.size(); i++) {
    if((dtofrawhits[i]->dE>0)&&((dtofrawhits[i]->t*1000.)>tMin)&&(dtofrawhits[i]->t*1000.<trigTime)) {

      uint32_t E  = dtofrawhits[i]->dE*1000000.;   // in keV (from GeV)
      uint32_t t  = dtofrawhits[i]->t*1000.-tMin;  // in picoseconds
    
      if(noroot==0)tofEnergies->Fill(dtofrawhits[i]->dE*1000000.);
      if(noroot==0)tofTimes->Fill(dtofrawhits[i]->t-tMin/1000);

      hc++;
      hitCount++;
      nhits=1;
      cscRef cscADC      = DTOFRawHitTranslationADC(dtofrawhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC250;
      hit[0].module_mode = FADC250_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = E;  // in keV
      hit[0].hdata[1]    = t/FADC250tick;
      if(E>0x7ffff)cerr << "E too large for TOF: " << E << endl;
      
      if(dumphits>1) {
        jout << endl;
        jout << " TOF ADC plane,bar,lr are " << dtofrawhits[i]->plane << ", " << dtofrawhits[i]->bar 
             << ", " << dtofrawhits[i]->lr << endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for TOF ADC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }

      hitCount++;
      nhits=1;
      cscRef cscTDC      = DTOFRawHitTranslationTDC(dtofrawhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscTDC.crate;
      hit[0].slot_id     = cscTDC.slot;
      hit[0].chan_id     = cscTDC.channel;
      hit[0].module_id   = F1TDC32;
      hit[0].module_mode = 0;
      hit[0].nwords      = 1;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = t/F1TDC32tick;
      
      if(dumphits>1) {
        jout << endl;
        jout << " TOF TDC plane,bar,lr are " << dtofrawhits[i]->plane << ", " << dtofrawhits[i]->bar 
             << ", " << dtofrawhits[i]->lr << endl;
        jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot << ", " << cscTDC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }

      if(nomc2coda==0) {
       stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for TOF TDC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if((dumphits>=1)&&(hc>0)) {
    jout << endl << "TOF hits: " << hc << endl << endl;
  }


  // DBCALHit - FADC250 and F1TDC32 (60 ps)
  vector<const DBCALHit*> dbcalhits;
  eventLoop->Get(dbcalhits);
  sort(dbcalhits.begin(),dbcalhits.end(),compareDBCALHits);

  hc=0;
  for(i=0; i<dbcalhits.size(); i++) {
    if((dbcalhits[i]->E>0)&&((dbcalhits[i]->t*1000.)>tMin)&&(dbcalhits[i]->t*1000.<trigTime)) {

      uint32_t E     = dbcalhits[i]->E*1000.;       // in keV (original apparently in MeV ???)
      uint32_t t     = dbcalhits[i]->t*1000.-tMin;  // in picoseconds

      if(noroot==0)bcalEnergies->Fill(dbcalhits[i]->E*1000.);
      if(noroot==0)bcalTimes->Fill(dbcalhits[i]->t-tMin/1000);

      hc++;
      hitCount++;
      nhits=1;
      cscRef cscADC      = DBCALHitTranslationADC(dbcalhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC250;
      hit[0].module_mode = FADC250_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = E/10;  // in 10 keV units
      hit[0].hdata[1]    = t/FADC250tick;
      if(E/10>0x7ffff)cerr << "E too large for BCAL: " << E << endl;
      
      if(dumphits>1) {
        jout << endl;
        jout << " BCAL ADC module,sector,layer,end are " << dbcalhits[i]->module<< ", " << dbcalhits[i]->sector
             << ", " << dbcalhits[i]->layer << ", " << dbcalhits[i]->end << endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }

      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for BCAL ADC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }

      hitCount++;
      nhits=1;
      cscRef cscTDC      = DBCALHitTranslationTDC(dbcalhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscTDC.crate;
      hit[0].slot_id     = cscTDC.slot;
      hit[0].chan_id     = cscTDC.channel;
      hit[0].module_id   = F1TDC32;
      hit[0].module_mode = 0;
      hit[0].nwords      = 1;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = t/F1TDC32tick;
      
      if(dumphits>1) {
        jout << endl;
        jout << " BCAL TDC module,sector,layer,end are " << dbcalhits[i]->module<< ", " << dbcalhits[i]->sector
             << ", " << dbcalhits[i]->layer << ", " << dbcalhits[i]->end << endl;
        jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot << ", " << cscTDC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for BCAL TDC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if((dumphits>=1)&&(hc>0)) {
    jout << endl << "BCAL hits: " << hc << endl << endl;
  }




  // DFCALHit - FADC250
  vector<const DFCALHit*> dfcalhits;
  eventLoop->Get(dfcalhits);
  sort(dfcalhits.begin(),dfcalhits.end(),compareDFCALHits);

  hc=0;
  for(i=0; i<dfcalhits.size(); i++) {
    if((dfcalhits[i]->E>0)&&((dfcalhits[i]->t*1000.)>tMin)&&(dfcalhits[i]->t*1000.<trigTime)) {

      uint32_t E     = dfcalhits[i]->E*1000.;       // in MeV (from GeV)
      uint32_t t     = dfcalhits[i]->t*1000.-tMin;  // in picoseconds
    
      if(noroot==0)fcalEnergies->Fill(dfcalhits[i]->E*1000000.);
      if(noroot==0)fcalTimes->Fill(dfcalhits[i]->t-tMin/1000);

      hc++;
      hitCount++;
      nhits=1;
      cscRef cscADC      = DFCALHitTranslationADC(dfcalhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC250;
      hit[0].module_mode = FADC250_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = E/10;  // in 10 MeV units
      hit[0].hdata[1]    = t/FADC250tick;
      if(E/10>0x7ffff)cerr << "E too large for FCAL: " << E << endl;
      
      if(dumphits>1) {
        jout << endl;
        jout << " FCAL ADC row,column are " << dfcalhits[i]->row << ", " << dfcalhits[i]->column << endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for FCAL ADC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if((dumphits>=1)&&(hc>0)) {
    jout << endl << "FCAL hits: " << hc << endl << endl;
  }


  // DFDCHit - cathode strips FADC125 or anode wires F1TDC48 (115 ps)
  vector<const DFDCHit*> dfdchits; 
  eventLoop->Get(dfdchits);
  sort(dfdchits.begin(),dfdchits.end(),compareDFDCHits);

  hc=0;
  for(i=0; i<dfdchits.size(); i++) {
    if((dfdchits[i]->q>0)&&((dfdchits[i]->t*1000.)>tMin)&&(dfdchits[i]->t*1000.<trigTime)) {

      uint32_t q      = dfdchits[i]->q;            // in femtoCoulombs
      uint32_t t      = dfdchits[i]->t*1000.-tMin; // in picoseconds
      
      if(noroot==0)fdcCharges->Fill(dfdchits[i]->q);
      if(noroot==0)fdcTimes->Fill(dfdchits[i]->t-tMin/1000);

      int type = dfdchits[i]->type;

      // FADC125
      if(type==1) {
        hc++;
        hitCount++;
        nhits=1;
        cscRef cscADC      = DFDCCathodeHitTranslation(dfdchits[i]);
        hit[0].hit_id      = hitCount;
        hit[0].det_id      = detID;
        hit[0].crate_id    = cscADC.crate;
        hit[0].slot_id     = cscADC.slot;
        hit[0].chan_id     = cscADC.channel;
        hit[0].module_id   = FADC125;
        hit[0].module_mode = FADC125_MODE_IP;
        hit[0].nwords      = 2;
        hit[0].hdata       = mcData;
        hit[0].hdata[0]    = q;  // in fC
        hit[0].hdata[1]    = t/FADC125tick;
        if(q>0x7ffff)cerr << "q too large for FDC: " << q << endl;
        
        if(dumphits>2) {
          jout << endl;
          jout << " FDC ADC gPlane,element are " << dfdchits[i]->gPlane << ", " << dfdchits[i]->element << endl;
          jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << endl;
          jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << endl;
          jout << " q,t are " << q << ", " << t << endl;
          jout << endl;
        }
        
        if(nomc2coda==0) {
          stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
          if(stat!=nhits) {
            jerr << "?error return from mc2codaWrite() for FDC ADC: " << stat << endl << endl;
            exit(EXIT_FAILURE);
          }
        }


      // F1TDC48
      } else if(type==0) {
        hitCount++;
        nhits=1;
        cscRef cscTDC      = DFDCAnodeHitTranslation(dfdchits[i]);
        hit[0].hit_id      = hitCount;
        hit[0].det_id      = detID;
        hit[0].crate_id    = cscTDC.crate;
        hit[0].slot_id     = cscTDC.slot;
        hit[0].chan_id     = cscTDC.channel;
        hit[0].module_id   = F1TDC48;
        hit[0].module_mode = 0;
        hit[0].nwords      = 1;
        hit[0].hdata       = mcData;
        hit[0].hdata[0]    = t/F1TDC48tick;
        
        if(dumphits>2) {
          jout << endl;
          jout << " FDC TDC gPlane,element are " << dfdchits[i]->gPlane << ", " << dfdchits[i]->element << endl;
          jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot << ", " << cscTDC.channel << endl;
          jout << " hdata is: " << hit[0].hdata[0] << endl;
          jout << " q,t are " << q << ", " << t << endl;
          jout << endl;
        }
      
        if(nomc2coda==0) {
          stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
          if(stat!=nhits) {
            jerr << "?error return from mc2codaWrite() for Fdc TDC: " << stat << endl << endl;
            exit(EXIT_FAILURE);
          }
        }
      }
    }
  }
  if((dumphits>=1)&&(hc>0)) {
    jout << endl << "FDC hits: " << hc << endl << endl;
  }



  // DSCHit - FADC250 and F1TDC32 (60 ps)
  vector<const DSCHit*> dsthits;
  eventLoop->Get(dsthits);
  sort(dsthits.begin(),dsthits.end(),compareDSTHits);

  hc=0;
  for(i=0; i<dsthits.size(); i++) {
    if((dsthits[i]->dE>0)&&((dsthits[i]->t*1000.)>tMin)&&(dsthits[i]->t*1000.<trigTime)) {

      uint32_t E     = dsthits[i]->dE*1000000.;   // in keV (from GeV)
      uint32_t t     = dsthits[i]->t*1000.-tMin;  // in picoseconds

      if(noroot==0)stEnergies->Fill(dsthits[i]->dE*1000000.);
      if(noroot==0)stTimes->Fill(dsthits[i]->t-tMin/1000);

      hc++;
      hitCount++;
      nhits=1;
      cscRef cscADC      = DSTHitTranslationADC(dsthits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC250;
      hit[0].module_mode = FADC250_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = E;  // in keV
      hit[0].hdata[1]    = t/FADC250tick;
      if(E>0x7ffff)cerr << "E too large for ST: " << E << endl;
      
      if(dumphits>1) {
        jout << endl;
        jout << " ST ADC sector is " << dsthits[i]->sector << endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for ST ADC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }


      hitCount++;
      nhits=1;
      cscRef cscTDC      = DSTHitTranslationTDC(dsthits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscTDC.crate;
      hit[0].slot_id     = cscTDC.slot;
      hit[0].chan_id     = cscTDC.channel;
      hit[0].module_id   = F1TDC32;
      hit[0].module_mode = 0;
      hit[0].nwords      = 1;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = t/F1TDC32tick;
      
      if(dumphits>1) {
        jout << endl;
        jout << " ST TDC sector is " << dsthits[i]->sector << endl;
        jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot << ", " << cscTDC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for ST TDC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if((dumphits>=1)&&(hc>0)) {
    jout << endl << "ST hits: " << hc << endl << endl;
  }


  // DTagger - FADC250 and F1TDC32 (60 ps)
  vector<const DTagger*> dtaggerhits;
  eventLoop->Get(dtaggerhits);
  sort(dtaggerhits.begin(),dtaggerhits.end(),compareDTaggerHits);

  hc=0;
  for(i=0; i<dtaggerhits.size(); i++) {
    if((dtaggerhits[i]->E>0)&&((dtaggerhits[i]->t*1000.)>tMin)&&(dtaggerhits[i]->t*1000.<trigTime)) {

      uint32_t E     = dtaggerhits[i]->E*1000000.;    // in keV
      uint32_t t     = dtaggerhits[i]->t*1000.-tMin;  // in picoseconds

      if(noroot==0)tagEnergies->Fill(dtaggerhits[i]->E*1000000.);
      if(noroot==0)tagTimes->Fill(dtaggerhits[i]->t-tMin/1000);

      hc++;
      hitCount++;
      nhits=1;
      cscRef cscADC      = DTaggerTranslationADC(dtaggerhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC250;
      hit[0].module_mode = FADC250_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = E/1000;  // in MeV
      hit[0].hdata[1]    = t/FADC250tick;
      if(E/1000>0x7ffff)cerr << "E too large for Tagger: " << E << endl;
      
      if(dumphits>1) {
        jout << endl;
        jout << " Tagger ADC row,column are " << dtaggerhits[i]->row << ", " << dtaggerhits[i]->column<< endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for TAGGER ADC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }


      hitCount++;
      nhits=1;
      cscRef cscTDC      = DTaggerTranslationTDC(dtaggerhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscTDC.crate;
      hit[0].slot_id     = cscTDC.slot;
      hit[0].chan_id     = cscTDC.channel;
      hit[0].module_id   = F1TDC32;
      hit[0].module_mode = 0;
      hit[0].nwords      = 1;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = t/F1TDC32tick;
      
      if(dumphits>1) {
        jout << endl;
        jout << " Tagger TDC row,column are " << dtaggerhits[i]->row << ", " << dtaggerhits[i]->column << endl;
        jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot << ", " << cscTDC.channel << endl;
        jout << " hdata is: " << hit[0].hdata[0] << endl;
        jout << " E,t are " << E << ", " << t << endl;
        jout << endl;
      }
      
      if(nomc2coda==0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if(stat!=nhits) {
          jerr << "?error return from mc2codaWrite() for Tagger TDC: " << stat << endl << endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if((dumphits>=1)&&(hc>0)) {
    jout << endl << "Tagger hits: " << hc << endl << endl;
  }



  // close event
  if(nomc2coda==0) {
    int nwords = mc2codaCloseEvent(eventID);
    if(nwords<0) {
      jerr << "?error return from mc2codaCloseEVent(): " << nwords << endl << endl;
      exit(EXIT_FAILURE);
    }
  }


  // write event
  pthread_mutex_lock(&rawMutex);
  chan->write(eventID->evbuf);
  pthread_mutex_unlock(&rawMutex);


  // not needed, using reset instead
  // free event
  // if(nomc2coda==0) {
  //   mc2codaFreeEvent(eventID);
  // }


  // done
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

  jout << "Reading translation table " << translationTableName << endl;

  // create parser and specify element handlers
  XML_Parser xmlParser = XML_ParserCreate(NULL);
  if(xmlParser==NULL) {
    jerr << endl << endl << "readTranslationTable...unable to create parser" << endl << endl;
    exit(EXIT_FAILURE);
  }
  XML_SetElementHandler(xmlParser,startElement,endElement);


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


int type2detID(string &type) {
  if(type=="vmecpu") {
    return(VMECPU);
  } else if (type=="tid") {
    return(TID);
  } else if (type=="fadc250") {
    return(FADC250);
  } else if (type=="fadc125") {
    return(FADC125);
  } else if (type=="f1tdcv2") {
    return(F1TDC32);
  } else if (type=="f1tdcv3") {
    return(F1TDC48);
  } else if (type=="jldisc") {
    return(JLDISC);
  } else if (type=="vx1290a") {
    return(CAEN1290); // (should this be CAEN1290_MODE_TM ??)
  } else {
    return(USERMOD);
  }
}



void JEventProcessor_rawevent::startElement(void *userData, const char *xmlname, const char **atts) {
  
  static int crate=0, slot=0;

  static string type,Type;
  int mc2codaType;
  int channel = 0;
  string Detector;
  string end;
  string row,column,module,sector,layer,chan;
  string ring,straw,plane,bar,gPlane,element;


  // store crate summary info, fill both maps
  if(strcasecmp(xmlname,"halld_online_translation_table")==0) {
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
        type = string(atts[i+1]);
        std::transform(type.begin(), type.end(), type.begin(), (int(*)(int)) tolower);
      }
    }

    mc2codaType = type2detID(type);
    if(mc2codaType!=USERMOD) {
      nModules[nCrate-1]++;
      modules[nCrate-1][slot-1] = mc2codaType;
      if((mc2codaType==VMECPU) || (mc2codaType==TID)) {
        detID[nCrate-1][slot-1]   = 0;
      } else {
        detID[nCrate-1][slot-1]   = 1;
      }
    }


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

	// ignore certain module types
	if(type == "disc") return;
	if(type == "ctp") return;
	if(type == "sd") return;
	if(type == "a1535sn") return;


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
      if(type=="f1tdcv2") {
        s = "bcaltdc::";
      } else if (type=="fadc250") {
        s = "bcaladc::";
      } else {
        s = "unknownBCAL::";
        jerr << endl << endl << "?startElement...illegal type for BCAL: " << Type << " ("<<type<<")" << endl << endl;
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
        
      
    } else if(detector=="st") {
      if(type=="f1tdcv2") {
        s = "sttdc::";
      } else if (type=="fadc250") {
        s = "stadc::";
      } else {
        s = "unknownST::";
        jerr << endl << endl << "?startElement...illegal type for ST: " << Type << endl << endl;
      }
      s += sector;
      cscMap[s] = csc;
    

    } else if(detector=="fdc_cathodes") {
      if(type=="fadc125") {
        s = "fdccathode::";
      } else {
        s = "unknownFDCCathode::";
        jerr << endl << endl << "?startElement...illegal type for FDC Cathode: " << Type << endl << endl;
      }
      s += gPlane + ":" + element;
      cscMap[s] = csc;

      
    } else if(detector=="fdc_wires") {
      if(type=="f1tdcv3") {
        s = "fdcanode::";
      } else {
        s = "unknownFDCAnode::";
        jerr << endl << endl << "?startElement...illegal type for FDC Anode: " << Type << endl << endl;
      }
      s += gPlane + ":" + element;
      cscMap[s] = csc;

      
    } else if(detector=="tof") {
      if(type=="vx1290a") {
        s = "toftdc::";
      } else if (type=="fadc250") {
        s = "tofadc::";
      } else {
        s = "unknownTOF::";
        jerr << endl << endl << "?startElement...illegal type for TOF: " << Type << endl << endl;
      }
      s += plane + ":" + bar + ":" + end;
      cscMap[s] = csc;

     
    } else if(detector=="tagh") {
	  if(type=="f1tdcv2") {
        s = "taghtdc::";
	  }else if (type=="fadc250") {
        s = "taghadc::";
      } else {
        s = "unknownTagger::";
        jerr << endl << endl << "?startElement...illegal type for TAGH: " << Type << endl << endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

    } else if(detector=="tagm") {
      if(type=="f1tdcv2") {
        s = "tagmtdc::";
      } else if (type=="fadc250") {
        s = "tagmadc::";
      } else {
        s = "unknownTagger::";
        jerr << endl << endl << "?startElement...illegal type for TAGM: " << Type << endl << endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

      
    } else if(detector=="psc") {
      if(type=="f1tdcv2") {
        s = "psctdc::";
      } else if (type=="fadc250") {
        s = "pscadc::";
      } else {
        s = "unknownPSC::";
        jerr << endl << endl << "?startElement...illegal type for PSC: " << Type << endl << endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

   } else if(detector=="ps") {
      if (type=="fadc250") {
        s = "psadc::";
      } else {
        s = "unknownPS::";
        jerr << endl << endl << "?startElement...illegal type for PS: " << Type << endl << endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

      
    } else {
      jerr << endl << endl << "?startElement...unknown detector " << Detector << endl << endl;
    }

	if(crate<0 || crate>=MAXDCRATE){
		jerr << " Crate value of "<<crate<<" is not in range 0 <= crate < " << MAXDCRATE << endl;
		exit(-1);
	}

	if(slot<0 || slot>=MAXDSLOT){
		jerr << " Slot value of "<<slot<<" is not in range 0 <= slot < " << MAXDSLOT << endl;
		exit(-1);
	}

	if(channel<0 || channel>=MAXDCHANNEL){
		jerr << " Crate value of "<<channel<<" is not in range 0 <= channel < " << MAXDCHANNEL << endl;
		exit(-1);
	}


    // fill detector map, index is crate,slot,channel
    detectorMap[crate][slot][channel] = s;
    


  } else {
    jerr << endl << endl << "?startElement...unknown xml tag " << xmlname << endl << endl;
  }

}
    

//--------------------------------------------------------------------------


void JEventProcessor_rawevent::endElement(void *userData, const char *xmlname) {
  // nothing to do yet...
}


//--------------------------------------------------------------------------



//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//  aux routines encode hit info into string for inverse lookup table
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTOFRawHitTranslationADC(const DTOFRawHit* hit) const {
  string end;
  if(hit->plane == 0){
    end = (hit->lr==0 ? "UP":"DW");
  }else{
    end = (hit->lr==0 ? "N":"S");
  }
  string s = "tofadc::" + lexical_cast<string>(hit->plane) + ":" + lexical_cast<string>(hit->bar)
    + ":" + end;
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTOFRawHitTranslationTDC(const DTOFRawHit* hit) const {
  string end;
  if(hit->plane == 0){
    end = (hit->lr==0 ? "UP":"DW");
  }else{
    end = (hit->lr==0 ? "N":"S");
  }
  string s = "toftdc::" + lexical_cast<string>(hit->plane) + ":" + lexical_cast<string>(hit->bar)
    + ":" + end;
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DBCALHitTranslationADC(const DBCALHit *hit) const {
  string end = hit->end==0 ? "U":"D";
  string s = "bcaladc::" + lexical_cast<string>(hit->module) + ":" + lexical_cast<string>(hit->sector)
    + ":" + lexical_cast<string>(hit->layer) + ":" + end;
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DBCALHitTranslationTDC(const DBCALHit *hit) const {
  string end = hit->end==0 ? "U":"D";
  string s = "bcaltdc::" + lexical_cast<string>(hit->module) + ":" + lexical_cast<string>(hit->sector)
    + ":" + lexical_cast<string>(hit->layer) + ":" + end;
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFCALHitTranslationADC(const DFCALHit* hit) const {
  string s = "fcaladc::" + lexical_cast<string>(hit->row) + ":" + lexical_cast<string>(hit->column);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFDCAnodeHitTranslation(const DFDCHit* hit) const {
  string s = "fdcanode::"  + lexical_cast<string>(hit->gPlane) + ":" + lexical_cast<string>(hit->element);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFDCCathodeHitTranslation(const DFDCHit* hit) const {
  string s = "fdccathode::"  + lexical_cast<string>(hit->gPlane) + ":" + lexical_cast<string>(hit->element);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DCDCHitTranslationADC(const DCDCHit* hit) const {
  string s = "cdcadc::" + lexical_cast<string>(hit->ring) + ":" + lexical_cast<string>(hit->straw);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DSTHitTranslationADC(const DSCHit* hit) const {
  string s = "stadc::" + lexical_cast<string>(hit->sector);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DSTHitTranslationTDC(const DSCHit* hit) const {
  string s = "sttdc::" + lexical_cast<string>(hit->sector);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTaggerTranslationTDC(const DTagger* hit) const {
  string s = "taggertdc::" + lexical_cast<string>(hit->row) +":" + lexical_cast<string>(hit->column);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTaggerTranslationADC(const DTagger* hit) const {
  string s = "taggeradc::" + lexical_cast<string>(hit->row) +":" + lexical_cast<string>(hit->column);
  if(cscMap.count(s)<=0)jerr << "?unknown map entry " << s << endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//--------------------------------------------------------------------------
