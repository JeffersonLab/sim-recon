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
//      115 ps/count for F1TDC48
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
static TH1F *tagmEnergies;
static TH1F *taghEnergies;

static TH1F *fdcCharges;
static TH1F *cdcCharges;

static TH1F *tofTimes;
static TH1F *fdcTimes;
static TH1F *cdcTimes;
static TH1F *fcalTimes;
static TH1F *bcalTimes;
static TH1F *stTimes;
static TH1F *tagmTimes;
static TH1F *taghTimes;


//  is this thread-safe?
extern "C" {
#include "mc2coda.h"
}


#include<sstream>
#include<iomanip>
#include<algorithm>
#include<fstream>
#include <expat.h>


// Replace use of BOOST lexical_cast with a simple templated function
// so that BOOST is no longer required.   1/24/2014  DL
//#include <boost/lexical_cast.hpp>
//using namespace boost;
template<typename T>
string lexical_cast(T t)
{
   stringstream ss;
   ss << t;
   return ss.str();
}


// to protect writing to output file
static pthread_mutex_t rawMutex = PTHREAD_MUTEX_INITIALIZER;


// for evio output and translation table
static evioFileChannel *chan       = NULL;
static string fileBase             = "rawevent";
static string outputFileName;
static string translationTableName = "tt.xml";


// current run number
static int runNumber;
static unsigned int user_runNumber=0xdeadbeef; // specified via configuration parameter


// csc map converts from detector spec to (crate,slot,channel)
//  key is detector-dependent encoded string (e.g. "cdcadc::2:25" for CDC ring 2 straw 25)
//  value is struct containing (crate,slot,channel)
static map<string,cscVal> cscMap;

// Useful "non-value" for a cscRef
static cscVal CDCBAL_NULL = {-1, -1, -1};
static cscRef CSCREF_NULL = CDCBAL_NULL;

// detector map (inverse of csc map) is 3-dimensional array of strings with indices (crate,slot,channel)
//  content is detector-dependent encoded string
#define MAXDCRATE   78+1
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
static int maxCrateNum = 0;
static int crateID[MAXCRATE];
static int nModules[MAXCRATE];
static int modules[MAXCRATE][MAXSLOT];
static int detID[MAXCRATE][MAXSLOT];


// misc
static uint64_t trigTime    = 32000000;    // in picoseconds
static float tMin           = -100000.;    // minimum hit time in picoseconds

// default time steps
static double trigtick      = 4000;    // in picoseconds
static double CAENTDCtick   = 25.;      // in picoseconds
static double F1TDC32tick   = 60.;     // in picoseconds
static double F1TDC48tick   = 115.;    // in picoseconds
static double FADC250tick   = 62.5;    // in picoseconds
static double FADC125tick   = 800.;    // in picoseconds

static double CDC_ADCscale  = 0.25E5/1.0E6;
static double FDC_ADCscale  = 1.3E5/2.4E4;
static double FCAL_ADCscale = 2.5E5/4.0E1;
static double BCAL_ADCscale = 10.;
static double TOF_ADCscale  = 5.2E5/0.2;
static double SC_ADCscale   = 5.2E-5/2.0E-2 ;

static double CDC_ADCtick = FADC125tick;
static double FDC_ADCtick = FADC125tick;
static double FCAL_ADCtick = FADC250tick;
static double BCAL_ADCtick = FADC250tick;
static double TOF_ADCtick = FADC250tick;
static double SC_ADCtick = FADC250tick;

static double FDC_TDCtick = F1TDC48tick;
static double BCAL_TDCtick = F1TDC32tick;
static double TOF_TDCtick = CAENTDCtick;
static double SC_TDCtick = F1TDC32tick;


// debug
static int dumphits  = 0;
static int nomc2coda = 0;
static int noroot    = 0;
static int dumpmap   = 0;

// Translation table
bool NO_CCDB = false;  // try and read from CCDB by default
string XML_FILENAME = "tt.xml"; // default filename for XML file if CCDB fails or is not used


//----------------------------------------------------------------------------


// required for JANA framework to find and identify this plugin
extern "C"{
void InitPlugin(JApplication *app) {
   InitJANAPlugin(app);
   app->AddProcessor(new JEventProcessor_rawevent());
}
} // "C"


//----------------------------------------------------------------------------

// Comparison operator for testing if two cscRefs are equal
bool operator == (cscRef a, cscRef b) {
   if (a.channel != b.channel)
      return false;
   if (a.slot != b.slot) 
      return false;
   return a.crate == b.crate;
}


// JEventProcessor_rawevent (Constructor) invoked once only
JEventProcessor_rawevent::JEventProcessor_rawevent() {

   // Default is to just read translation table from CCDB. If this fails,
   // then an attempt will be made to read from a file on the local disk.
   // The filename can be specified to be anything, but if the user specifies
   // this, then we assume that they want to use it and skip using the CCDB.
   // They may also specify that they want to skip checking the CCDB via
   // the "TT:NO_CCDB" parameter. This would only be useful if they want to
   // force the use of a local file named "tt.xml".
   gPARMS->SetDefaultParameter("TT:NO_CCDB", NO_CCDB, "Don't try getting"
           " translation table from CCDB and just look for file. Only useful"
           " if you want to force reading tt.xml. This is automatically set"
           " if you specify a different filename via the "
           "TT:XML_FILENAME parameter.");
   JParameter *p = gPARMS->SetDefaultParameter("TT:XML_FILENAME", XML_FILENAME,
           "Fallback filename of translation table XML file."
           " If set to non-default, CCDB will not be checked.");
   if (p->GetDefault() != p->GetValue())
      NO_CCDB = true;


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

  // option to dump map to file
  gPARMS->SetDefaultParameter("RAWEVENT:DUMPMAP",dumpmap,
          "Dump map of translation table map to file (for debugging)");

  // option to set run number
  gPARMS->SetDefaultParameter("RAWEVENT:RUNNUMBER",user_runNumber,
          "Override run number from input file with this one"
          " which will be written to every event in output file");

}


//----------------------------------------------------------------------------


// ~JEventProcessor_rawevent (Destructor) once only
JEventProcessor_rawevent::~JEventProcessor_rawevent() {
  if (nomc2coda == 0) {
    mc2codaFree(expID);
  }
}


//----------------------------------------------------------------------------


// init called once-only at beginning, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::init(void) {


  // read translation table, fill crate id arrays
  readTranslationTable();

  // initialize mc2coda package
  if (nomc2coda == 0) {
    expID=mc2codaInitExp(maxCrateNum+1,expName.c_str());
    if (expID == NULL) {
      jerr << "?NULL return from mc2codaInitExp()" << std::endl;
      exit(EXIT_FAILURE);
    }
  }


  // feed crate-specific info to mc2coda package,
  // note that VMECPU and TID not included in module count
  int stat;
  if (nomc2coda == 0) {
    for (int i=0; i<nCrate; i++) {
     int nMod = nModules[i]-2;
     if (nMod < 1)
        continue; // ignore crates with no digitization modules

      stat = mc2codaSetCrate(expID,crateID[i],nModules[i]-2,modules[i],detID[i]);
      if (stat == -1) {
        jerr << "?JEventProcessor_rawevent...error return from mc2codaSetCrate()"
             << std::endl << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }
  
  // Optionally dump translation table map into file for debugging
  if (dumpmap) {
     ofstream *ofs = new ofstream("cscMap.out");
     if (ofs) {
        jout << "Dumping translation table map to cscMap.out ..." << std::endl;
        map<string,cscVal>::iterator iter = cscMap.begin();
        for (; iter != cscMap.end(); iter++) {
           cscVal &csc = iter->second;
           *ofs << iter->first << " " << csc.crate << " " << csc.slot 
                << " " << csc.channel << std::endl;
        }
        ofs->close();
        delete ofs;
     }
  }

  // root histograms
  if (noroot == 0) {
    tofEnergies  = new TH1F("tofe",  "TOF energies in keV",1000,0.,5000.);
    stEnergies   = new TH1F("ste",   "ST energies in keV",1000,0.,5000.);
    fcalEnergies = new TH1F("fcale", "FCAL energies in keV",1000,0.,500000.);
    bcalEnergies = new TH1F("bcale", "BCAL energies in keV",1000,0.,1000000.);
    tagmEnergies = new TH1F("tagme", "Tagger microscope energies in keV",1000,0.,40000.);
    taghEnergies = new TH1F("taghe", "Tagger hodoscope array energies in keV",1000,0.,40000.);

    fdcCharges   = new TH1F("fdcq",  "FDC charges in fC",1000,0.,1500.);
    cdcCharges   = new TH1F("cdcq",  "CDC charges in fC",1000,0.,25000.);

    tofTimes     = new TH1F("toft","  TOF times in nsec",1000,0.,500.-tMin/1000);
    stTimes      = new TH1F("stt",   "ST times in nsec",1000,0.,250.-tMin/1000);
    fdcTimes     = new TH1F("fdct",  "FDC times in nsec",1000,0.,2000.-tMin/1000);
    cdcTimes     = new TH1F("cdct",  "CDC times in nsec",1000,0.,1000.-tMin/1000);
    bcalTimes    = new TH1F("bcalt", "BCAL times in nsec",1000,0.,200.-tMin/1000);
    fcalTimes    = new TH1F("fcalt", "FCAL times in nsec",1000,0.,200.-tMin/1000);
    tagmTimes    = new TH1F("tagmt", "Tagger microscope times in nsec",1000,0.,250.-tMin/1000);
    taghTimes    = new TH1F("tahft", "Tagger hodoscope times in nsec",1000,0.,250.-tMin/1000);
  }

  return NOERROR;
}

//----------------------------------------------------------------------------


// brun called once-only at beginning of run, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::brun(JEventLoop *eventLoop, int runnumber) {

  runNumber=runnumber;
  jout << std::endl << "   brun called for run " << runNumber << std::endl;
  if (user_runNumber != 0xdeadbeef) {
    jout << "   *** overriding with user-supplied run number: " 
         << user_runNumber << " ***" << std::endl;
    runNumber = runnumber = user_runNumber;
  }
  mc2codaSetRunNumber(runNumber);

  // close old output file
  if (chan != NULL) {
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
  jout << std::endl << "   opening output file:   " << outputFileName << std::endl << std::endl << std::endl;

  
  // load scale factors for converting from physical units into detector hit units
  // the converstion factors are set to reasonable default values when they are defined
  // but try to load them from the CCDB so that we are using the same factors to create
  // the EVIO files as are used to read them in
  jout << "Loading ADC/TDC scale factors..." << endl;

  map<string,double> scale_factors;
  if (eventLoop->GetCalib("/CDC/digi_scales", scale_factors))
	  jout << "Error loading /CDC/digi_scales !" << endl;
  if ( scale_factors.find("CDC_ADC_ASCALE") != scale_factors.end() ) {
	  CDC_ADCscale = 1. / scale_factors["CDC_ADC_ASCALE"];
  } else {
	  jerr << "Unable to get CDC_ADC_ASCALE from /CDC/digi_scales !" << endl;
  }
  if ( scale_factors.find("CDC_ADC_TSCALE") != scale_factors.end() ) {
	  CDC_ADCtick = 1000. * scale_factors["CDC_ADC_TSCALE"];
  } else {
	  jerr << "Unable to get CDC_ADC_TSCALE from /CDC/digi_scales !" << endl;
  }

  if (eventLoop->GetCalib("/FDC/digi_scales", scale_factors))
	  jout << "Error loading /FDC/digi_scales !" << endl;
  if ( scale_factors.find("FDC_ADC_ASCALE") != scale_factors.end() ) {
	  FDC_ADCscale = 1. / scale_factors["FDC_ADC_ASCALE"];
  } else {
	  jerr << "Unable to get FDC_ADC_ASCALE from /FDC/digi_scales !" << endl;
  }
  if ( scale_factors.find("FDC_ADC_TSCALE") != scale_factors.end() ) {
	  FDC_ADCtick = 1000. * scale_factors["FDC_ADC_TSCALE"];
  } else {
	  jerr << "Unable to get FDC_ADC_TSCALE from /FDC/digi_scales !" << endl;
  }
  if ( scale_factors.find("FDC_TDC_SCALE") != scale_factors.end() ) {
	  FDC_TDCtick = 1000. * scale_factors["FDC_TDC_SCALE"];
  } else {
	  jerr << "Unable to get FDC_TDC_SCALE from /FDC/digi_scales !" << endl;
  }

  if (eventLoop->GetCalib("/FCAL/digi_scales", scale_factors))
	  jout << "Error loading /FCAL/digi_scales !" << endl;
  if ( scale_factors.find("FCAL_ADC_ASCALE") != scale_factors.end() ) {
	  FCAL_ADCscale = 1. / scale_factors["FCAL_ADC_ASCALE"];
  } else {
	  jerr << "Unable to get FCAL_ADC_ASCALE from /FCAL/digi_scales !" << endl;
  }
  if ( scale_factors.find("FCAL_ADC_TSCALE") != scale_factors.end() ) {
	  FCAL_ADCtick = 1000. * scale_factors["FCAL_ADC_TSCALE"];
  } else {
	  jerr << "Unable to get FCAL_ADC_TSCALE from /FCAL/digi_scales !" << endl;
  }

  if (eventLoop->GetCalib("/BCAL/digi_scales", scale_factors))
	  jout << "Error loading /BCAL/digi_scales !" << endl;
  if ( scale_factors.find("BCAL_ADC_ASCALE") != scale_factors.end() ) {
	  BCAL_ADCscale = 1. / scale_factors["BCAL_ADC_ASCALE"];
  } else {
	  jerr << "Unable to get BCAL_ADC_ASCALE from /BCAL/digi_scales !" << endl;
  }
  if ( scale_factors.find("BCAL_ADC_TSCALE") != scale_factors.end() ) {
	  BCAL_ADCtick = 1000. * scale_factors["BCAL_ADC_TSCALE"];
  } else {
	  jerr << "Unable to get BCAL_ADC_TSCALE from /BCAL/digi_scales !" << endl;
  }
  if ( scale_factors.find("BCAL_TDC_SCALE") != scale_factors.end() ) {
	  BCAL_TDCtick = 1000. * scale_factors["BCAL_TDC_SCALE"];
  } else {
	  jerr << "Unable to get BCAL_TDC_SCALE from /BCAL/digi_scales !" << endl;
  }

  if (eventLoop->GetCalib("/TOF/digi_scales", scale_factors))
	  jout << "Error loading /TOF/digi_scales !" << endl;
  if ( scale_factors.find("TOF_ADC_ASCALE") != scale_factors.end() ) {
	  TOF_ADCscale = 1. / scale_factors["TOF_ADC_ASCALE"];
  } else {
	  jerr << "Unable to get TOF_ADC_ASCALE from /TOF/digi_scales !" << endl;
  }
  if ( scale_factors.find("TOF_ADC_TSCALE") != scale_factors.end() ) {
	  TOF_ADCtick = 1000. * scale_factors["TOF_ADC_TSCALE"];
  } else {
	  jerr << "Unable to get TOF_ADC_TSCALE from /TOF/digi_scales !" << endl;
  }
  if ( scale_factors.find("TOF_TDC_SCALE") != scale_factors.end() ) {
	  TOF_TDCtick = 1000. * scale_factors["TOF_TDC_SCALE"];
  } else {
	  jerr << "Unable to get TOF_TDC_SCALE from /TOF/digi_scales !" << endl;
  }

  if (eventLoop->GetCalib("/START_COUNTER/digi_scales", scale_factors))
	  jout << "Error loading /START_COUNTER/digi_scales !" << endl;
  if ( scale_factors.find("SC_ADC_ASCALE") != scale_factors.end() ) {
	  SC_ADCscale = 1. / scale_factors["SC_ADC_ASCALE"];
  } else {
	  jerr << "Unable to get SC_ADC_ASCALE from /START_COUNTER/digi_scales !" << endl;
  }
  if ( scale_factors.find("SC_ADC_TSCALE") != scale_factors.end() ) {
	  SC_ADCtick = 1000. * scale_factors["SC_ADC_TSCALE"];
  } else {
	  jerr << "Unable to get SC_ADC_TSCALE from /START_COUNTER/digi_scales !" << endl;
  }
  if ( scale_factors.find("SC_TDC_SCALE") != scale_factors.end() ) {
	  SC_TDCtick = 1000. * scale_factors["SC_TDC_SCALE"];
  } else {
	  jerr << "Unable to get SC_TDC_SCALE from /START_COUNTER/digi_scales !" << endl;
  }


  // add header event if required
  // ...
//  map<string,cscVal>::iterator iter = cscMap.begin();
//  for (; iter != cscMap.end(); iter++) {
//     string key = iter->first;
//   if (key.find("toftdc") == 0) _DBG_ << key << std::endl;
//  }
  
  
  return NOERROR;
}


//----------------------------------------------------------------------------


  static bool compareDTOFHits(const DTOFHit* h1, const DTOFHit* h2) {
    if (h1->plane != h2->plane) {
      return(h1->plane<h2->plane);
    } else if (h1->bar != h2->bar) {
      return(h1->bar<h2->bar);
    } else if (h1->end != h2->end) {
      return(h1->end<h2->end);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDFCALHits(const DFCALHit* h1, const DFCALHit* h2) {
    if (h1->row != h2->row) {
      return(h1->row<h2->row);
    } else if (h1->column != h2->column) {
      return(h1->column<h2->column);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDBCALHits(const DBCALHit* h1, const DBCALHit* h2) {
    if (h1->module != h2->module) {
      return(h1->module<h2->module);
    } else if (h1->sector != h2->sector) {
      return(h1->sector<h2->sector);
    } else if (h1->layer != h2->layer) {
      return(h1->layer<h2->layer);
    } else if (h1->end != h2->end) {
      return(h1->end<h2->end);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDBCALTDCHits(const DBCALTDCHit* h1, const DBCALTDCHit* h2) {
    if (h1->module != h2->module) {
      return(h1->module<h2->module);
    } else if (h1->sector != h2->sector) {
      return(h1->sector<h2->sector);
    } else if (h1->layer != h2->layer) {
      return(h1->layer<h2->layer);
    } else if (h1->end != h2->end) {
      return(h1->end<h2->end);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDFDCHits(const DFDCHit* h1, const DFDCHit* h2) {
    if (h1->gPlane != h2->gPlane) {
      return(h1->gPlane<h2->gPlane);
    } else if (h1->element != h2->element) {
      return(h1->element<h2->element);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDCDCHits(const DCDCHit* h1, const DCDCHit* h2) {
    if (h1->ring != h2->ring) {
      return(h1->ring<h2->ring);
    } else if (h1->straw != h2->straw) {
      return(h1->straw<h2->straw);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDSTHits(const DSCHit* h1, const DSCHit* h2) {
    if (h1->sector != h2->sector) {
      return(h1->sector<h2->sector);
    } else {
      return(h1->t<h2->t);
    }
  }

  static bool compareDTAGMHits(const DTAGMHit* h1, const DTAGMHit* h2) {
    if (h1->column != h2->column) {
      return (h1->column < h2->column);
    }
    else if (h1->row != h2->row) {
      return (h1->row < h2->row);
    }
    else if (h1->t != h2->t) {
      return (h1->t <h2->t);
    }
    else {
      return (h1->time_fadc < h2->time_fadc);
    }
  }

  static bool compareDTAGHHits(const DTAGHHit* h1, const DTAGHHit* h2) {
    if (h1->counter_id != h2->counter_id) {
      return (h1->counter_id < h2->counter_id);
    }
    else if (h1->t != h2->t) {
      return (h1->t < h2->t);
    }
    else {
      return (h1->time_fadc < h2->time_fadc);
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
  if (nomc2coda == 0) {
    if (first_time) {
      first_time=false;
      eventID = mc2codaOpenEvent(expID, (uint64_t)eventnumber, trigTime/trigtick, eventType, MAXEVENTSIZE);
      if (eventID == NULL) {
        jerr << "?NULL return from mc2codaOpenEvent()" << std::endl << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      int stat = mc2codaResetEvent(eventID, (uint64_t)eventnumber, trigTime/trigtick, eventType);
      if (stat != 0) {
        jerr << "?error return from mc2codaResetEvent()" << std::endl << std::endl;
        exit(EXIT_FAILURE);
      }
    }
  }


  // get all raw hit data, sort according to id and time order before feeding to mc2coda package
  // scale data appropriately to minimize loss of precision by using full available bit range

  // FADC125   pulse integral = 17 bits or about 1.3E5
  //           pulse_time     = 14 bits or about 1.6E4
  //
  // FADC250   pulse_integral = 19 bits or about 5.2E5
  //           pulse_time     = 16 bits or about 6.5E4
  //
  // F21TDC    time           = 7.8 us for 115ps resolution
  //                          = 3.9 us for  60ps resolution


  // DCDCHit - FADC125
  vector<const DCDCHit*> dcdchits; 
  eventLoop->Get(dcdchits);
  sort(dcdchits.begin(),dcdchits.end(),compareDCDCHits);

  hc=0;
  for (i=0; i<dcdchits.size(); i++) {
    if ((dcdchits[i]->q > 0) && ((dcdchits[i]->t*1000.) > tMin) &&
        (dcdchits[i]->t*1000. < trigTime))
    {
      //uint32_t q     = dcdchits[i]->q * (1./5.18) * (1.3E5/1.0E6); // q is in femtoCoulombs (max is ~1E6) 
      uint32_t q     = dcdchits[i]->q * CDC_ADCscale; // q is in femtoCoulombs (max is ~1E6) 
      uint32_t t     = dcdchits[i]->t*1000.0 -tMin;    // t is in nanoseconds (max is ~900ns)
      
      if (noroot == 0)
         cdcCharges->Fill(dcdchits[i]->q);
      if (noroot == 0)
         cdcTimes->Fill(dcdchits[i]->t-tMin/1000);

      cscRef cscADC = DCDCHitTranslationADC(dcdchits[i]);
      if (cscADC == CSCREF_NULL)
         continue;
      hc++;
      hitCount++;
      nhits=1;
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC125;
      hit[0].module_mode = FADC125_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = q;  // in fADC counts
      hit[0].hdata[1]    = static_cast<double>(t)/CDC_ADCtick;
      //if (q > 0x7ffff) std::cerr << "q too large for CDC: " << q << std::endl;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " CDC ADC ring,straw are " << dcdchits[i]->ring 
             << ", " << dcdchits[i]->straw << std::endl;
        jout << " c,s,c are " << cscADC.crate << ", " 
             << cscADC.slot << ", " << cscADC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " 
             << hit[0].hdata[1] << std::endl;
        jout << " q,t are " << q << ", " << t << std::endl;
        jout << std::endl;
      }
      
      if (nomc2coda == 0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for CDC ADC: "
               << stat << std::endl << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }      
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "CDC hits: " << hc << std::endl << std::endl;
  }


  // DTOFHit - FADC250 and CAEN TDC (25 ps)
  vector<const DTOFHit*> dtofrawhits; 
  eventLoop->Get(dtofrawhits);
  sort(dtofrawhits.begin(),dtofrawhits.end(),compareDTOFHits);

  hc=0;
  for (i=0; i < dtofrawhits.size(); i++) {
    if ((dtofrawhits[i]->dE > 0) && ((dtofrawhits[i]->t*1000.) > tMin) &&
        (dtofrawhits[i]->t*1000. < trigTime))
    {
	//uint32_t E  = dtofrawhits[i]->dE*(5.2E5/0.2);   // E is GeV (max ~0.1)
      uint32_t E  = dtofrawhits[i]->dE*TOF_ADCscale;   // E is GeV (max ~0.1)
      uint32_t t  = dtofrawhits[i]->t*1000.-tMin;  // in picoseconds
    
      if (noroot == 0)
         tofEnergies->Fill(dtofrawhits[i]->dE*1000000.);
      if (noroot == 0)
         tofTimes->Fill(dtofrawhits[i]->t-tMin/1000);

      hc++;
      hitCount++;
      nhits=1;
      cscRef cscADC      = DTOFHitTranslationADC(dtofrawhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC250;
      hit[0].module_mode = FADC250_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = E;  // in fADC counts
      hit[0].hdata[1]    = static_cast<double>(t)/TOF_ADCtick;
      if (E > 0x7ffff)
         std::cerr << "E too large for TOF: " << E << std::endl;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " TOF ADC plane,bar,lr are " << dtofrawhits[i]->plane << ", " << dtofrawhits[i]->bar 
             << ", " << dtofrawhits[i]->end << std::endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", " << cscADC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }
      
      if (nomc2coda == 0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for TOF ADC: " << stat << std::endl << std::endl;
          exit(EXIT_FAILURE);
        }
      }

      hitCount++;
      nhits=1;
      cscRef cscTDC      = DTOFHitTranslationTDC(dtofrawhits[i]);
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscTDC.crate;
      hit[0].slot_id     = cscTDC.slot;
      hit[0].chan_id     = cscTDC.channel;
      hit[0].module_id   = CAEN1290;
      hit[0].module_mode = 0;
      hit[0].nwords      = 1;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = static_cast<double>(t)/TOF_TDCtick;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " TOF TDC plane,bar,lr are " << dtofrawhits[i]->plane << ", " << dtofrawhits[i]->bar 
             << ", " << dtofrawhits[i]->end << std::endl;
        jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot << ", " << cscTDC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }

      if (nomc2coda == 0) {
       stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for TOF TDC: " << stat << std::endl << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "TOF hits: " << hc << std::endl << std::endl;
  }


  // DBCALHit - FADC250 and F1TDC32 (60 ps)
  vector<const DBCALHit*> dbcalhits;
  eventLoop->Get(dbcalhits);
  sort(dbcalhits.begin(),dbcalhits.end(),compareDBCALHits);

  hc=0;
  for (i=0; i < dbcalhits.size(); i++) {
    if ((dbcalhits[i]->E > 0) && ((dbcalhits[i]->t*1000.) > tMin) &&
        (dbcalhits[i]->t*1000. < trigTime))
    {
      // All calorimeter hits have E in GeV, but BCAL_ADCscale assumes MeV
      // so fix that here by converting E to MeV first, then apply scale factor.
      uint32_t E = dbcalhits[i]->E*1000.*BCAL_ADCscale;  // (each fADC count ~ 100keV) (max ~2.5E4)
      uint32_t t = dbcalhits[i]->t*1000.-tMin;     // in picoseconds

      if (noroot == 0)
         bcalEnergies->Fill(dbcalhits[i]->E*1000.);
      if (noroot == 0)
         bcalTimes->Fill(dbcalhits[i]->t-tMin/1000);

      cscRef cscADC = DBCALHitTranslationADC(dbcalhits[i]);
      if (cscADC == CSCREF_NULL)
         continue;
      hc++;
      hitCount++;
      nhits=1;
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscADC.crate;
      hit[0].slot_id     = cscADC.slot;
      hit[0].chan_id     = cscADC.channel;
      hit[0].module_id   = FADC250;
      hit[0].module_mode = FADC250_MODE_IP;
      hit[0].nwords      = 2;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = E;
      hit[0].hdata[1]    = static_cast<double>(t)/BCAL_ADCtick;
      if (E/10 > 0x7ffff)
         std::cerr << "E too large for BCAL: " << E << std::endl;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " BCAL ADC module,sector,layer,end are " 
             << dbcalhits[i]->module << ", " << dbcalhits[i]->sector
             << ", " << dbcalhits[i]->layer << ", " 
             << dbcalhits[i]->end << std::endl;
        jout << " c,s,c are " << cscADC.crate << ", " 
             << cscADC.slot << ", " << cscADC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " 
             << hit[0].hdata[1] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }

      
      if (nomc2coda == 0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for BCAL ADC: " 
               << stat << std::endl << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      /**
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
      hit[0].hdata[0]    = static_cast<double>(t)/F1TDC32tick;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " BCAL TDC module,sector,layer,end are " << dbcalhits[i]->module << ", " << dbcalhits[i]->sector
             << ", " << dbcalhits[i]->layer << ", " << dbcalhits[i]->end << std::endl;
        jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot << ", " << cscTDC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }
      
      if (nomc2coda == 0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for BCAL TDC: " << stat << std::endl << std::endl;
          exit(EXIT_FAILURE);
        }
      }
      **/
    }
  }
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "BCAL hits: " << hc << std::endl << std::endl;
  }


  // BCAL TDC hits are handled in separate objects
  vector<const DBCALTDCHit*> dbcaltdchits;
  eventLoop->Get(dbcaltdchits);
  sort(dbcaltdchits.begin(),dbcaltdchits.end(),compareDBCALTDCHits);
  
  hc=0;
  for (i=0; i<dbcaltdchits.size(); i++) {
    if (((dbcaltdchits[i]->t*1000.) > tMin) && 
         (dbcaltdchits[i]->t*1000. < trigTime))
    {
      int32_t E     = 0;
      int32_t t     = dbcaltdchits[i]->t*1000.-tMin;  // in picoseconds

      if (noroot == 0)
         bcalTimes->Fill(dbcaltdchits[i]->t-tMin/1000);

      cscRef cscTDC = DBCALHitTranslationTDC(dbcaltdchits[i]);
      if (cscTDC == CSCREF_NULL)
         continue;

      hc++;
      hitCount++;
      nhits=1;
      hit[0].hit_id      = hitCount;
      hit[0].det_id      = detID;
      hit[0].crate_id    = cscTDC.crate;
      hit[0].slot_id     = cscTDC.slot;
      hit[0].chan_id     = cscTDC.channel;
      hit[0].module_id   = F1TDC32;
      hit[0].module_mode = 0;
      hit[0].nwords      = 1;
      hit[0].hdata       = mcData;
      hit[0].hdata[0]    = static_cast<double>(t)/BCAL_TDCtick;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " BCAL TDC module,sector,layer,end are " 
             << dbcaltdchits[i]->module << ", " << dbcaltdchits[i]->sector
             << ", " << dbcaltdchits[i]->layer << ", " 
             << dbcaltdchits[i]->end << std::endl;
        jout << " c,s,c are " << cscTDC.crate << ", " 
             << cscTDC.slot << ", " << cscTDC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }
      
      if (nomc2coda == 0) {
         stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
         if (stat != nhits) {
            jerr << "?error return from mc2codaWrite() for BCAL TDC: " 
                 << stat << std::endl << std::endl;
            exit(EXIT_FAILURE);
         }
      }
    }
  }
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "BCAL TDC hits: " << hc << std::endl << std::endl;
  }


  // DFCALHit - FADC250
  vector<const DFCALHit*> dfcalhits;
  eventLoop->Get(dfcalhits);
  sort(dfcalhits.begin(),dfcalhits.end(),compareDFCALHits);

  hc=0;
  for (i=0; i<dfcalhits.size(); i++) {
    if ((dfcalhits[i]->E > 0) && ((dfcalhits[i]->t*1000.) > tMin) &&
        (dfcalhits[i]->t*1000. < trigTime))
    {
      uint32_t E     = dfcalhits[i]->E*FCAL_ADCscale;  // E in GeV (max ~4)
      uint32_t t     = dfcalhits[i]->t*1000.-tMin;  // in picoseconds
      
      if (noroot == 0)
         fcalEnergies->Fill(dfcalhits[i]->E*1000000.);
      if (noroot == 0)
         fcalTimes->Fill(dfcalhits[i]->t-tMin/1000);

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
      hit[0].hdata[0]    = E; 
      hit[0].hdata[1]    = static_cast<double>(t)/FCAL_ADCtick;
      if (E/10 > 0x7ffff)
         std::cerr << "E too large for FCAL: " << E << std::endl;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " FCAL ADC row,column are " << dfcalhits[i]->row 
             << ", " << dfcalhits[i]->column << std::endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot
             << ", " << cscADC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", "
             << hit[0].hdata[1] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }
      
      if (nomc2coda == 0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for FCAL ADC: " 
               << stat << std::endl << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "FCAL hits: " << hc << std::endl << std::endl;
  }


  // DFDCHit - cathode strips FADC125 or anode wires F1TDC48 (115 ps)
  vector<const DFDCHit*> dfdchits; 
  eventLoop->Get(dfdchits);
  sort(dfdchits.begin(),dfdchits.end(),compareDFDCHits);

  hc=0;
  for (i=0; i < dfdchits.size(); i++) {
    if ((dfdchits[i]->q > 0) && ((dfdchits[i]->t*1000.) > tMin) &&
        (dfdchits[i]->t*1000. < trigTime))
    {
      uint32_t q = dfdchits[i]->q*FDC_ADCscale; // for cathodes
      if (dfdchits[i]->type == 0)
         q = 0.0; // No amplitude read for wires 
      uint32_t t = dfdchits[i]->t*1000.-tMin; // in picoseconds
      
      if (noroot == 0)
         fdcCharges->Fill(dfdchits[i]->q);
      if (noroot == 0)fdcTimes->Fill(dfdchits[i]->t-tMin/1000);

      int type = dfdchits[i]->type;
      // FADC125
      if (type == 1) {
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
        hit[0].hdata[0]    = q; 
        hit[0].hdata[1]    = static_cast<double>(t)/FDC_ADCtick;
        if (q > 0x7ffff)
           std::cerr << "q too large for FDC: " << q << std::endl;
        
        if (dumphits > 2) {
          jout << std::endl;
          jout << " FDC ADC gPlane,element are " << dfdchits[i]->gPlane 
               << ", " << dfdchits[i]->element << std::endl;
          jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot 
               << ", " << cscADC.channel << std::endl;
          jout << " hdata is: " << hit[0].hdata[0] << ", " 
               << hit[0].hdata[1] << std::endl;
          jout << " q,t are " << q << ", " << t << std::endl;
          jout << std::endl;
        }
        
        if (nomc2coda == 0) {
          stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
          if (stat != nhits) {
            jerr << "?error return from mc2codaWrite() for FDC ADC: " 
                 << stat << std::endl << std::endl;
            exit(EXIT_FAILURE);
          }
        }

      // F1TDC48
      }
      else if (type == 0) {
        hc++;
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
        hit[0].hdata[0]    = static_cast<double>(t)/FDC_TDCtick;
        
        if (dumphits > 2) {
          jout << std::endl;
          jout << " FDC TDC gPlane,element are " << dfdchits[i]->gPlane 
               << ", " << dfdchits[i]->element << std::endl;
          jout << " c,s,c are " << cscTDC.crate << ", " 
               << cscTDC.slot << ", " << cscTDC.channel << std::endl;
          jout << " hdata is: " << hit[0].hdata[0] << std::endl;
          jout << " q,t are " << q << ", " << t << std::endl;
          jout << std::endl;
        }
      
        if (nomc2coda == 0) {
          stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
          if (stat != nhits) {
            jerr << "?error return from mc2codaWrite() for Fdc TDC: " 
                 << stat << std::endl << std::endl;
            exit(EXIT_FAILURE);
          }
        }
      }
    } 
  }

  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "FDC hits: " << hc << std::endl << std::endl;
  }


  // DSCHit - FADC250 and F1TDC32 (60 ps)
  vector<const DSCHit*> dsthits;
  eventLoop->Get(dsthits);
  sort(dsthits.begin(),dsthits.end(),compareDSTHits);

  hc=0;
  for (i=0; i < dsthits.size(); i++) {
    if ((dsthits[i]->dE > 0) && ((dsthits[i]->t*1000.) > tMin) &&
        (dsthits[i]->t*1000. < trigTime))
    {
      uint32_t E     = dsthits[i]->dE*SC_ADCscale;   // dE in GeV (max ~2E-2)
      uint32_t t     = dsthits[i]->t*1000.-tMin;  // in picoseconds

      if (noroot == 0)
         stEnergies->Fill(dsthits[i]->dE*1000000.);
      if (noroot == 0)
         stTimes->Fill(dsthits[i]->t-tMin/1000);

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
      hit[0].hdata[0]    = E; 
      hit[0].hdata[1]    = static_cast<double>(t)/SC_ADCtick;
      if (E > 0x7ffff)
         std::cerr << "E too large for ST: " << E << std::endl;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " ST ADC sector is " << dsthits[i]->sector << std::endl;
        jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot
             << ", " << cscADC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << ", " 
             << hit[0].hdata[1] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }
      
      if (nomc2coda == 0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for ST ADC: " 
               << stat << std::endl << std::endl;
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
      hit[0].hdata[0]    = static_cast<double>(t)/SC_TDCtick;
      
      if (dumphits > 1) {
        jout << std::endl;
        jout << " ST TDC sector is " << dsthits[i]->sector << std::endl;
        jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot 
             << ", " << cscTDC.channel << std::endl;
        jout << " hdata is: " << hit[0].hdata[0] << std::endl;
        jout << " E,t are " << E << ", " << t << std::endl;
        jout << std::endl;
      }
      
      if (nomc2coda == 0) {
        stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
        if (stat != nhits) {
          jerr << "?error return from mc2codaWrite() for ST TDC: " 
               << stat << std::endl << std::endl;
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "ST hits: " << hc << std::endl << std::endl;
  }


  // DTAGM - FADC250 and F1TDC32 (60 ps)
  vector<const DTAGMHit*> dtagmhits;
  eventLoop->Get(dtagmhits);
  sort(dtagmhits.begin(),dtagmhits.end(),compareDTAGMHits);

  hc=0;
  for (i=0; i < dtagmhits.size(); i++) {
    if (dtagmhits[i]->npix_fadc > 0 && 
        ((dtagmhits[i]->t*1000 > tMin && dtagmhits[i]->t*1000 < trigTime) ||
         (dtagmhits[i]->time_fadc*1000 > tMin &&
          dtagmhits[i]->time_fadc*1000 < trigTime)) )
    {
      uint32_t pA = dtagmhits[i]->npix_fadc + 2500;       // in SiPM pixels
      uint32_t t = dtagmhits[i]->t*1000 - tMin;           // in ps
      uint32_t tA = dtagmhits[i]->time_fadc*1000 - tMin;  // in ps

      if (noroot == 0)
        tagmEnergies->Fill(dtagmhits[i]->npix_fadc);
      if (noroot == 0)
        tagmTimes->Fill(dtagmhits[i]->time_fadc-tMin/1000);

      cscRef cscADC = DTAGMHitTranslationADC(dtagmhits[i]);
      if (! (cscADC == CSCREF_NULL)) {
        hc++;
        hitCount++;
        nhits=1;
        hit[0].hit_id      = hitCount;
        hit[0].det_id      = detID;
        hit[0].crate_id    = cscADC.crate;
        hit[0].slot_id     = cscADC.slot;
        hit[0].chan_id     = cscADC.channel;
        hit[0].module_id   = FADC250;
        hit[0].module_mode = FADC250_MODE_IP;
        hit[0].nwords      = 2;
        hit[0].hdata       = mcData;
        hit[0].hdata[0]    = pA;   // in SiPM pixels
        hit[0].hdata[1]    = static_cast<double>(tA)/FADC250tick;
        if (pA > 0x7ffff)
          std::cerr << "pA too large for tagger microscope: " 
                    << pA << std::endl;
  
        if (dumphits > 1) {
           jout << std::endl;
           jout << " Tagger microscope ADC row,column are " 
                << dtagmhits[i]->row << ", " << dtagmhits[i]->column << std::endl;
           jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", "
                << cscADC.channel << std::endl;
           jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1]
                << std::endl;
           jout << " pA,tA are " << pA << ", " << tA << std::endl;
           jout << std::endl;
        }
 
        if (nomc2coda == 0) {
           stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
           if (stat != nhits) {
             jerr << "?error return from mc2codaWrite() for TAGGER ADC: " 
                  << stat << std::endl << std::endl;
             exit(EXIT_FAILURE);
           }
        }
      }

      cscRef cscTDC = DTAGMHitTranslationTDC(dtagmhits[i]);
      if (! (cscADC == CSCREF_NULL)) {
        hitCount++;
        nhits=1;
        hit[0].hit_id      = hitCount;
        hit[0].det_id      = detID;
        hit[0].crate_id    = cscTDC.crate;
        hit[0].slot_id     = cscTDC.slot;
        hit[0].chan_id     = cscTDC.channel;
        hit[0].module_id   = F1TDC32;
        hit[0].module_mode = 0;
        hit[0].nwords      = 1;
        hit[0].hdata       = mcData;
        hit[0].hdata[0]    = static_cast<double>(t)/F1TDC32tick;
 
        if (dumphits > 1) {
           jout << std::endl;
           jout << " Tagger microscope TDC row,column are " 
                << dtagmhits[i]->row << ", " << dtagmhits[i]->column << std::endl;
           jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot 
                << ", " << cscTDC.channel << std::endl;
           jout << " hdata is: " << hit[0].hdata[0] << std::endl;
           jout << " t is " << t << std::endl;
           jout << std::endl;
        }
  
        if (nomc2coda == 0) {
           stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
           if (stat != nhits) {
              jerr << "?error return from mc2codaWrite() for Tagger TDC: " 
                   << stat << std::endl << std::endl;
              exit(EXIT_FAILURE);
           }
        }
      }
    }
  }
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "Tagger microscope hits: " << hc << std::endl << std::endl;
  }

  // DTAGH - FADC250 and F1TDC32 (60 ps)
  vector<const DTAGHHit*> dtaghhits;
  eventLoop->Get(dtaghhits);
  sort(dtaghhits.begin(),dtaghhits.end(),compareDTAGHHits);

  hc=0;
  for (i=0; i < dtaghhits.size(); i++) {
    if (dtaghhits[i]->npe_fadc > 0 && 
        ((dtaghhits[i]->t*1000 > tMin && dtaghhits[i]->t*1000 < trigTime) ||
         (dtaghhits[i]->time_fadc*1000 > tMin &&
          dtaghhits[i]->time_fadc*1000 < trigTime)) )
    {
      uint32_t pA = dtaghhits[i]->npe_fadc + 2500;       // in SiPM pixels
      uint32_t t = dtaghhits[i]->t*1000 - tMin;           // in ps
      uint32_t tA = dtaghhits[i]->time_fadc*1000 - tMin;  // in ps

      if (noroot == 0)
        taghEnergies->Fill(dtaghhits[i]->npe_fadc);
      if (noroot == 0)
        taghTimes->Fill(dtaghhits[i]->time_fadc-tMin/1000);

      cscRef cscADC = DTAGHHitTranslationADC(dtaghhits[i]);
      if (! (cscADC == CSCREF_NULL)) {
        hc++;
        hitCount++;
        nhits=1;
        hit[0].hit_id      = hitCount;
        hit[0].det_id      = detID;
        hit[0].crate_id    = cscADC.crate;
        hit[0].slot_id     = cscADC.slot;
        hit[0].chan_id     = cscADC.channel;
        hit[0].module_id   = FADC250;
        hit[0].module_mode = FADC250_MODE_IP;
        hit[0].nwords      = 2;
        hit[0].hdata       = mcData;
        hit[0].hdata[0]    = pA;   // in SiPM pixels
        hit[0].hdata[1]    = static_cast<double>(tA)/FADC250tick;
        if (pA > 0x7ffff)
          std::cerr << "pA too large for tagger hodoscope: " << pA << std::endl;
  
        if (dumphits > 1) {
           jout << std::endl;
           jout << " Tagger hodoscope ADC counter_id is " 
                << dtaghhits[i]->counter_id << std::endl;
           jout << " c,s,c are " << cscADC.crate << ", " << cscADC.slot << ", "
                << cscADC.channel << std::endl;
           jout << " hdata is: " << hit[0].hdata[0] << ", " << hit[0].hdata[1]
                << std::endl;
           jout << " pA,tA are " << pA << ", " << tA << std::endl;
           jout << std::endl;
        }
 
        if (nomc2coda == 0) {
           stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
           if (stat != nhits) {
             jerr << "?error return from mc2codaWrite() for TAGGER ADC: " 
                  << stat << std::endl << std::endl;
             exit(EXIT_FAILURE);
           }
        }
      }

      cscRef cscTDC = DTAGHHitTranslationTDC(dtaghhits[i]);
      if (! (cscADC == CSCREF_NULL)) {
        hitCount++;
        nhits=1;
        hit[0].hit_id      = hitCount;
        hit[0].det_id      = detID;
        hit[0].crate_id    = cscTDC.crate;
        hit[0].slot_id     = cscTDC.slot;
        hit[0].chan_id     = cscTDC.channel;
        hit[0].module_id   = F1TDC32;
        hit[0].module_mode = 0;
        hit[0].nwords      = 1;
        hit[0].hdata       = mcData;
        hit[0].hdata[0]    = static_cast<double>(t)/F1TDC32tick;
 
        if (dumphits > 1) {
           jout << std::endl;
           jout << " Tagger hodoscope TDC counter_id is " 
                << dtaghhits[i]->counter_id << std::endl;
           jout << " c,s,c are " << cscTDC.crate << ", " << cscTDC.slot 
                << ", " << cscTDC.channel << std::endl;
           jout << " hdata is: " << hit[0].hdata[0] << std::endl;
           jout << " t is " << t << std::endl;
           jout << std::endl;
        }
  
        if (nomc2coda == 0) {
           stat = mc2codaWrite(eventID,nhits,(struct coda_hit_info *)&hit[0]);
           if (stat != nhits) {
             jerr << "?error return from mc2codaWrite() for Tagger TDC: " 
                  << stat << std::endl << std::endl;
             exit(EXIT_FAILURE);
           }
        }
      }
    }
  }
  if ((dumphits >= 1) && (hc > 0)) {
    jout << std::endl << "Tagger hodoscope hits: " << hc << std::endl << std::endl;
  }

  // close event
  if (nomc2coda == 0) {
    int nwords = mc2codaCloseEvent(eventID);
    if (nwords < 0) {
      jerr << "?error return from mc2codaCloseEVent(): " << nwords 
           << std::endl << std::endl;
      exit(EXIT_FAILURE);
    }
  }

  // write event
  pthread_mutex_lock(&rawMutex);
  chan->write(eventID->evbuf);
  pthread_mutex_unlock(&rawMutex);

  return NOERROR;
}


//----------------------------------------------------------------------------


// erun called once-only at end of run, independent of the number of processing threads
jerror_t JEventProcessor_rawevent::erun(void) {

  jout << std::endl << "   erun called for run " << runNumber << std::endl << std::endl;


  // add end event if required
  // ...


  // close evio output file and delete channel
  if (chan != NULL) {
    chan->close();
    delete(chan);
    chan=NULL;
    jout << std::endl << "  output file " << outputFileName << " closed" << std::endl << std::endl;
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

  jout << "Reading translation table " << translationTableName << std::endl;


   // Get the calibration object
   JCalibration *jcalib = japp->GetJCalibration(1);
   
   //--------------------------------------------------------------
   // (this block cut and pasted from TTab/DTranslationTable.cc)
   
   // String to hold entire XML translation table
   string tt_xml; 

   // Try getting it from CCDB first
   if (jcalib && !NO_CCDB) {
      map<string,string> tt;
      string namepath = "Translation/DAQ2detector";
      jout << "Reading translation table from calib DB: " << namepath << " ..." << std::endl;
      jcalib->GetCalib(namepath, tt);
      if (tt.size() != 1) {
         jerr << " Error: Unexpected translation table format!" << std::endl;
         jerr << "        tt.size()=" << tt.size() << " (expected 1)" << std::endl;
      }else{
         // Copy table into tt string
         map<string,string>::iterator iter = tt.begin();
         tt_xml = iter->second;
      }
   }
   
   // If getting from CCDB fails, try just reading in local file
   if (tt_xml.size() == 0) {
      if (!NO_CCDB) jout << "Unable to get translation table from CCDB." << std::endl;
      jout << "Will try reading TT from local file: " << XML_FILENAME << std::endl;

      // Open file
      ifstream ifs(XML_FILENAME.c_str());
      if (! ifs.is_open()) {
         jerr << " Error: Cannot open file! Translation table unavailable." << std::endl;
         exit(-1);
      }
      
      // read lines into stringstream object 
      stringstream ss;
      while (ifs.good()) {
         char line[4096];
         ifs.getline(line, 4096);
         ss << line;
      }

      // Close file
      ifs.close();
      
      // Copy from stringstream to tt
      tt_xml = ss.str();
   }

   // create parser and specify element handlers
   XML_Parser xmlParser = XML_ParserCreate(NULL);
   if (xmlParser == NULL) {
      jerr << "readTranslationTable...unable to create parser" << std::endl;
      exit(EXIT_FAILURE);
   }
   XML_SetElementHandler(xmlParser,StartElement,EndElement);

   // Parse XML string
   int status=XML_Parse(xmlParser, tt_xml.c_str(), tt_xml.size(), 1); // "1" indicates this is the final piece of XML
   if (status == 0) {
      jerr << "  ?readTranslationTable...parseXMLFile parse error for " << XML_FILENAME << std::endl;
      jerr << XML_ErrorString(XML_GetErrorCode(xmlParser)) << std::endl;
   }

   //--------------------------------------------------------------


   XML_ParserFree(xmlParser);
}


//----------------------------------------------------------------------------


int type2detID(string &type) {
  if (type == "vmecpu" || type == "cpu") {
    return(VMECPU);
  } else if (type == "tid" || type == "ti") {
    return(TID);
  } else if (type == "fadc250") {
    return(FADC250);
  } else if (type == "fadc125") {
    return(FADC125);
  } else if (type == "f1tdcv2") {
    return(F1TDC32);
  } else if (type == "f1tdcv3") {
    return(F1TDC48);
  } else if (type == "jldisc" || type == "disc") {
    return(JLAB_DISC);
  } else if (type == "vx1290a") {
    return(CAEN1290);
  } else {
    return(UNKNOWN);
  }
}



void JEventProcessor_rawevent::StartElement(void *userData, const char *xmlname, const char **atts) {
  
  static int crate=0, slot=0;

  static string type,Type;
  int mc2codaType;
  int channel = 0;
  string Detector;
  string end;
  string id,row,column,module,sector,layer,chan;
  string ring,straw,plane,bar,gPlane,element;
  string package,chamber,view,strip,wire;


  // store crate summary info, fill both maps
  if (strcasecmp(xmlname,"halld_online_translation_table") == 0) {
    // do nothing


  } else if (strcasecmp(xmlname,"crate") == 0) {
    for (int i=0; atts[i]; i+=2) {
      if (strcasecmp(atts[i],"number") == 0) {
        crate = atoi(atts[i+1]);
        break;
      }
    }
    nCrate++;
    crateID[nCrate-1]=crate;
    nModules[nCrate-1]=0;
   if (crate > maxCrateNum) maxCrateNum = crate;

  } else if (strcasecmp(xmlname,"slot") == 0) {
    for (int i=0; atts[i]; i+=2) {
      if (strcasecmp(atts[i],"number") == 0) {
        slot = atoi(atts[i+1]);
      } else if (strcasecmp(atts[i],"type") == 0) {
        Type = string(atts[i+1]);
        type = string(atts[i+1]);
        std::transform(type.begin(), type.end(), type.begin(), (int(*)(int)) tolower);
      }
    }

   // The detID value set here shows up in the header of the Data Block Bank
   // of the output file. It should be set to one if this crate has JLab
   // made modules that output in the standard format (see document:
   // "VME Data Format Standards for JLAB Modules"). These would include
   // f250ADC, f125ADC, F1TDC, .... Slots containing other types of modules
   // (e.g. CAEN1290) should have their own unique detID. We use detID of
   // zero for non-digitizing modules like CPUs nd TIDs even though potentially,
   // one could read data from these.
    mc2codaType = type2detID(type);
    if (mc2codaType != UNKNOWN) {
      nModules[nCrate-1]++;
      modules[nCrate-1][slot-1] = mc2codaType;
     switch(mc2codaType) {
      case FADC250:
      case FADC125:
      case F1TDC32:
      case F1TDC48:
         detID[nCrate-1][slot-1]   = 1;
         break;
      case CAEN1190:
      case CAEN1290:
         detID[nCrate-1][slot-1]   = 20;
         break;
      default:
         detID[nCrate-1][slot-1]   = 0;
     }
    }


  } else if (strcasecmp(xmlname,"channel") == 0) {
      
    for (int i=0; atts[i]; i+=2) {
      if (strcasecmp(atts[i],"number") == 0) {
        channel = atoi(atts[i+1]);
      } else if (strcasecmp(atts[i],"detector") == 0) {
        Detector = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"id") == 0) {
        id = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"row") == 0) {
        row = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"column") == 0) {
        column = string(atts[i+1]);
       } else if (strcasecmp(atts[i],"col") == 0) {
        column = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"module") == 0) {
        module = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"sector") == 0) {
        sector = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"layer") == 0) {
        layer = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"end") == 0) {
        end = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"chan") == 0) {
        chan = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"ring") == 0) {
        ring = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"straw") == 0) {
        straw = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"gPlane") == 0) {
        gPlane = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"element") == 0) {
        element = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"plane") == 0) {
        plane = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"bar") == 0) {
        bar = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"package") == 0) {
        package = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"chamber") == 0) {
        chamber = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"view") == 0) {
        view = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"strip") == 0) {
        strip = string(atts[i+1]);
      } else if (strcasecmp(atts[i],"wire") == 0) {
        wire = string(atts[i+1]);
      }
    }

   // ignore certain module types
   if (type == "disc")
      return;
   if (type == "ctp")
      return;
   if (type == "sd")
      return;
   if (type == "a1535sn")
      return;

    // fill maps

    cscVal csc = {crate,slot,channel};
    string detector = Detector;
    std::transform(detector.begin(), detector.end(), detector.begin(),
                   (int(*)(int)) tolower);
   
    string s="unknown::";

    if (detector == "fcal") {
      if (type == "fadc250") {
        s = "fcaladc::";
      } else {
        s = "unknownFCAL::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for FCAL: " 
             << Type << std::endl << std::endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;
      

    } else if (detector == "bcal") {
      if (type == "f1tdcv2") {
        s = "bcaltdc::";
      } else if (type == "fadc250") {
        s = "bcaladc::";
      } else {
        s = "unknownBCAL::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for BCAL: " 
             << Type << " (" << type << ")" << std::endl << std::endl;
      }
      s += module + ":" + sector + ":" + layer + ":" + end;
      cscMap[s] = csc;

    } else if (detector == "cdc") {
      if (type == "fadc125") {
        s = "cdcadc::";
      } else {
        s = "unknownCDC::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for CDC: " 
             << Type << std::endl << std::endl;
      }
      s += ring + ":" + straw;
      cscMap[s] = csc;
        
      
    } else if (detector == "st") {
      if (type == "f1tdcv2") {
        s = "sttdc::";
      } else if (type == "fadc250") {
        s = "stadc::";
      } else if (type == "iseg") {
        s = "stiseg::"; // this just here to prevent warning message below
      } else {
        s = "unknownST::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for ST: " 
             << Type << std::endl << std::endl;
      }
      s += sector;
      cscMap[s] = csc;
 
    } else if (detector == "fdc_cathodes") {
      int ipackage = atoi(package.c_str());
      int ichamber = atoi(chamber.c_str());
      int igPlane = 1 + (ipackage-1)*6*3 + (ichamber-1)*3 + (view == "U" ? 0:2);
      stringstream ss;
      ss << igPlane;
      gPlane = ss.str();
      if (type == "fadc125") {
        s = "fdccathode::";
      } else {
        s = "unknownFDCCathode::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for FDC Cathode: " 
             << Type << std::endl << std::endl;
      }
      s += gPlane + ":" + strip;
      cscMap[s] = csc;

      
    } else if (detector == "fdc_wires") {
      int ipackage = atoi(package.c_str());
      int ichamber = atoi(chamber.c_str());
      int igPlane = 2 + (ipackage-1)*6*3 + (ichamber-1)*3;
      stringstream ss;
      ss << igPlane;
      gPlane = ss.str();
      if (type == "f1tdcv3") {
        s = "fdcanode::";
      } else {
        s = "unknownFDCAnode::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for FDC Anode: " 
             << Type << std::endl << std::endl;
      }
      s += gPlane + ":" + wire;
      cscMap[s] = csc;
      
    } else if (detector == "tof") {
      if (type == "vx1290a") {
        s = "toftdc::";
      } else if (type == "fadc250") {
        s = "tofadc::";
      } else {
        s = "unknownTOF::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for TOF: " 
             << Type << std::endl << std::endl;
      }
      s += plane + ":" + bar + ":" + end;
      cscMap[s] = csc;
     
    } else if (detector == "tagh") {
      if (type == "f1tdcv2") {
        s = "taghtdc::";
      } else if (type == "fadc250") {
        s = "taghadc::";
      } else {
        s = "unknownTagger::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for TAGH: " 
             << Type << std::endl << std::endl;
      }
      s += id;
      cscMap[s] = csc;

    } else if (detector == "tagm") {
      if (type == "f1tdcv2") {
        s = "tagmtdc::";
      } else if (type == "fadc250") {
        s = "tagmadc::";
      } else {
        s = "unknownTagger::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for TAGM: " 
             << Type << std::endl << std::endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

    } else if (detector == "psc") {
      if (type == "f1tdcv2") {
        s = "psctdc::" + id;
      } else if (type == "fadc250") {
        s = "pscadc::" + id;
      } else {
        s = "unknownPSC::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for PSC: " 
             << Type << std::endl << std::endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

   } else if (detector == "ps") {
      if (type == "fadc250") {
        s = "psadc::" + id;
      } else {
        s = "unknownPS::";
        jerr << std::endl << std::endl 
             << "?startElement...illegal type for PS: " 
             << Type << std::endl << std::endl;
      }
      s += row + ":" + column;
      cscMap[s] = csc;

    } else {
      jerr << std::endl << std::endl 
           << "?startElement...unknown detector " 
           << Detector << std::endl << std::endl;
    }

   if (crate < 0 || crate >= MAXDCRATE) {
      jerr << " Crate value of " << crate 
           << " is not in range 0 <= crate < " << MAXDCRATE << std::endl;
      exit(-1);
   }

   if (slot < 0 || slot >= MAXDSLOT) {
      jerr << " Slot value of " << slot 
           << " is not in range 0 <= slot < " << MAXDSLOT << std::endl;
      exit(-1);
   }

   if (channel < 0 || channel >= MAXDCHANNEL) {
      jerr << " Crate value of " << channel 
           << " is not in range 0 <= channel < " << MAXDCHANNEL << std::endl;
      exit(-1);
   }


    // fill detector map, index is crate,slot,channel
    detectorMap[crate][slot][channel] = s;

  } else {
    jerr << std::endl << std::endl 
         << "?startElement...unknown xml tag " 
         << xmlname << std::endl << std::endl;
  }

}


//--------------------------------------------------------------------------


void JEventProcessor_rawevent::EndElement(void *userData, const char *xmlname) {
  // nothing to do yet...
}


//--------------------------------------------------------------------------



//--------------------------------------------------------------------------
//--------------------------------------------------------------------------
//  aux routines encode hit info into string for inverse lookup table
//--------------------------------------------------------------------------
//--------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTOFHitTranslationADC(const DTOFHit* hit) const {
  string end;
  if (hit->plane == 0) {
    end = (hit->end == 0 ? "UP":"DW");
  }else{
    end = (hit->end == 0 ? "N":"S");
  }
  string s = "tofadc::" + lexical_cast(hit->plane) + ":"
                        + lexical_cast(hit->bar) + ":" + end;
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTOFHitTranslationTDC(const DTOFHit* hit) const {
  string end;
  if (hit->plane == 0) {
    end = (hit->end == 0 ? "UP":"DW");
  }else{
    end = (hit->end == 0 ? "N":"S");
  }
  string s = "toftdc::" + lexical_cast(hit->plane) + ":"
                        + lexical_cast(hit->bar) + ":" + end;
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DBCALHitTranslationADC(const DBCALHit *hit) const {
  string end = hit->end == 0 ? "U":"D";
  string s = "bcaladc::" + lexical_cast(hit->module) + ":" 
                         + lexical_cast(hit->sector) + ":" 
                         + lexical_cast(hit->layer) + ":" + end;
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


//cscRef JEventProcessor_rawevent::DBCALHitTranslationTDC(const DBCALHit *hit) const {
cscRef JEventProcessor_rawevent::DBCALHitTranslationTDC(const DBCALTDCHit *hit) const {
  // BCAL does not have TDC channels for layer 4, but some older simulation files
  // have this. Ignore those hits here.
  if (hit->layer > 3)
     return CSCREF_NULL;
  string end = hit->end == 0 ? "U":"D";
  string s = "bcaltdc::" + lexical_cast(hit->module) + ":" 
                         + lexical_cast(hit->sector) + ":" 
                         + lexical_cast(hit->layer) + ":" + end;
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFCALHitTranslationADC(const DFCALHit* hit) const {
  string s = "fcaladc::" + lexical_cast(hit->row) + ":" 
                         + lexical_cast(hit->column);
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFDCAnodeHitTranslation(const DFDCHit* hit) const {
  string s = "fdcanode::"  + lexical_cast(hit->gPlane) + ":" 
                           + lexical_cast(hit->element);
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DFDCCathodeHitTranslation(const DFDCHit* hit) const {
  string s = "fdccathode::"  + lexical_cast(hit->gPlane) + ":"
                             + lexical_cast(hit->element);
  if (cscMap.count(s) <= 0)
    jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DCDCHitTranslationADC(const DCDCHit* hit) const {
  string s = "cdcadc::" + lexical_cast(hit->ring) + ":" 
                       + lexical_cast(hit->straw);
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DSTHitTranslationADC(const DSCHit* hit) const {
  string s = "stadc::" + lexical_cast(hit->sector);
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DSTHitTranslationTDC(const DSCHit* hit) const {
  string s = "sttdc::" + lexical_cast(hit->sector);
  if (cscMap.count(s) <= 0)
     jerr << "?unknown map entry " << s << std::endl;
  return(cscMap[s]);
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTAGMHitTranslationTDC(const DTAGMHit* hit) const {
  if ( hit->column > 100)
    return CSCREF_NULL;
  string s = "tagmtdc::" + lexical_cast(hit->row) + ":"
                         + lexical_cast(hit->column);
  if (cscMap.count(s) <= 0)
    jerr << "?unknown map entry " << s << std::endl;
  return cscMap[s];
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTAGMHitTranslationADC(const DTAGMHit* hit) const {
  if ( hit->column > 100)
    return CSCREF_NULL;
  string s = "tagmadc::" + lexical_cast(hit->row) + ":"
                         + lexical_cast(hit->column);
  if (cscMap.count(s) <= 0)
    jerr << "?unknown map entry " << s << std::endl;
  return cscMap[s];
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTAGHHitTranslationTDC(const DTAGHHit* hit) const {
  if ( hit->counter_id > 274)
    return CSCREF_NULL;
  string s = "taghtdc::" + lexical_cast(hit->counter_id);
  if (cscMap.count(s) <= 0)
    jerr << "?unknown map entry " << s << std::endl;
  return cscMap[s];
}


//----------------------------------------------------------------------------


cscRef JEventProcessor_rawevent::DTAGHHitTranslationADC(const DTAGHHit* hit) const {
  if ( hit->counter_id > 274)
    return CSCREF_NULL;
  string s = "taghadc::" + lexical_cast(hit->counter_id);
  if (cscMap.count(s) <= 0)
    jerr << "?unknown map entry " << s << std::endl;
  return cscMap[s];
}


//----------------------------------------------------------------------------
//----------------------------------------------------------------------------

//--------------------------------------------------------------------------
