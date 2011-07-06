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
//
//
// To run:
//
//     $ hd_ana -PPLUGINS=rawevent inputFile.hddm
//
//
//
// still to do:
//    translation table
//    pair spectrometer
//
//
//
//  ejw, 27-jun-2011



#include "rawevent/JEventProcessor_rawevent.h"

#include<sstream>
#include<iomanip>
#include <expat.h>


// to protect writing to output file
static pthread_mutex_t rawMutex = PTHREAD_MUTEX_INITIALIZER;


// for evio output
static evioFileChannel *chan       = NULL;
static int evioBufSize             = 750000;
static string fileBase             = "rawevent";
static string translationTableName = "translation.xml";
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
  gPARMS->SetDefaultParameter("RAWEVENT:FILEBASE",fileBase);

  // get translation table file name from command line 
  gPARMS->SetDefaultParameter("RAWEVENT:TRANSLATION",translationTableName);

  // read translation table
  readTranslationTable();

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



  // DTOFRawHit - FADC250 and F1TDC (60 ps)
  vector<const DTOFRawHit*> dtofrawhits; 
  eventLoop->Get(dtofrawhits);
  for(i=0; i<dtofrawhits.size(); i++) {
    float dE  = dtofrawhits[i]->dE;
    float t   = dtofrawhits[i]->t;

    // translate to crate/slot/channel
    cscVal cscADC  = DTOFRawHitTranslationADC(dtofrawhits[i]);
    cscVal cscTDC  = DTOFRawHitTranslationTDC(dtofrawhits[i]);

    int crateADC     = cscADC.get<0>();
    int slotADC      = cscADC.get<1>();
    int channelADC   = cscADC.get<2>();


    // do something...
  }


  // DBCALHit - FADC250 and F1TDC (60 ps)
  vector<const DBCALHit*> dbcalhits;
  eventLoop->Get(dbcalhits);
  for(i=0; i<dbcalhits.size(); i++) {
    float E      = dbcalhits[i]->E;
    float t      = dbcalhits[i]->t;

    // translate ADC to crate/slot/channel
    cscVal cscADC  = DBCALHitTranslationADC(dbcalhits[i]);
    cscVal cscTDC  = DBCALHitTranslationTDC(dbcalhits[i]);

    // do something...
  }      



  // DFCALHit - FADC250
  vector<const DFCALHit*> dfcalhits;
  eventLoop->Get(dfcalhits);
  for(i=0; i<dfcalhits.size(); i++) {
    float E      = dfcalhits[i]->E;
    float t      = dfcalhits[i]->t;
    
    // translate to crate/slot/channel
    cscVal cscADC  = DFCALHitTranslationADC(dfcalhits[i]);

    // do something...
  }      


  // DFDCHit - cathode strips FADC125 or anode wires F1TDC (115 ps)
  vector<const DFDCHit*> dfdchits; 
  eventLoop->Get(dfdchits);
  for(unsigned int i=0; i<dfdchits.size(); i++) {
    float q      = dfdchits[i]->q;
    float t      = dfdchits[i]->t;

    // translate to crate/slot/channel
    cscVal csc  = DFDCHitTranslation(dfdchits[i]);

    int type = dfdchits[i]->type;
    if(type==0) {           // F1TDC
      // do something...
    } else if(type==1) {    // FADC125
      // do something...
    }

  }


  // DCDCHit - FADC125
  vector<const DCDCHit*> dcdchits; 
  eventLoop->Get(dcdchits);
  for(i=0; i<dcdchits.size(); i++) {
    float dE     = dcdchits[i]->dE;
    float t      = dcdchits[i]->t;

    // translate to crate/slot/channel
    cscVal cscADC  = DCDCHitTranslationADC(dcdchits[i]);

    // do something...
  }      


  // DSCHit - FADC250 and F1TDC (60 ps)
  vector<const DSCHit*> dschits;
  eventLoop->Get(dschits);
  for(unsigned int i=0; i<dschits.size(); i++) {
    float dE     = dschits[i]->dE;
    float t      = dschits[i]->t;

    // translate to crate/slot/channel
    cscVal cscADC  = DSCHitTranslationADC(dschits[i]);
    cscVal cscTDC  = DSCHitTranslationTDC(dschits[i]);

    // do something...
  }      


  // DTagger - FADC250 and F1TDC (60 ps)
  vector<const DTagger*> dtaggerhits;
  eventLoop->Get(dtaggerhits);
  for(i=0; i<dtaggerhits.size(); i++) {
    float E      = dtaggerhits[i]->E;
    float t      = dtaggerhits[i]->t;

    // translate to crate/slot/channel
    cscVal cscADC  = DTaggerTranslationADC(dtaggerhits[i]);
    cscVal cscTDC  = DTaggerTranslationTDC(dtaggerhits[i]);

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


void JEventProcessor_rawevent::readTranslationTable(void) {

  // create parser and specify start element handler
  XML_Parser xmlParser = XML_ParserCreate(NULL);
  XML_SetElementHandler(xmlParser,startElement,NULL);


  // open and parse the file
  FILE *f = fopen(translationTableName.c_str(),"r");
  if(f!=NULL) {
    int status,len;
    bool done;
    const int bufSize = 50000;
    char buf[bufSize];
    do {
      len  = fread(buf,1,bufSize,f);
      done = len!=bufSize;
      status=XML_Parse(xmlParser,buf,len,done);
      if((!done)&&(status==0)) {
        jerr << endl << endl << endl << "  ?readTranslationTable...parseXMLFile parse error for " << translationTableName
             << endl << endl << XML_ErrorString(XML_GetErrorCode(xmlParser))
             << endl << endl << endl;
        fclose(f);
        XML_ParserFree(xmlParser);
        return;
      }
    } while (!done);
    fclose(f);
    jout << endl << endl << "Successfully read translation table " << translationTableName << endl << endl << endl;

  } else {
    jerr << endl << endl << endl << "  ?readTranslationTable...unable to open " << translationTableName
         << endl << endl << endl;
  }
  XML_ParserFree(xmlParser);

}


//----------------------------------------------------------------------------


void JEventProcessor_rawevent::startElement(void *userData, const char *xmlname, const char **atts) {
  
  static int crate=0, slot=0;
  static string slotType=NULL;

  int channel = 0;
  string type,Detector,detector,end;
  string row,column,module,sector,layer,chan;
  string ring,straw,package,strip,wire,plane,bar;


  if(strcasecmp(xmlname,"translation_table")==0) {
    // do nothing

  } else if(strcasecmp(xmlname,"crate")==0) {
    for (int i=0; atts[i]; i+=2) {
      if(strcasecmp(atts[i],"number")==0) {
        crate = atoi(atts[i+1]);
        break;
      }
    }

  } else if(strcasecmp(xmlname,"slot")==0) {
    for (int i=0; atts[i]; i+=2) {
      if(strcasecmp(atts[i],"number")==0) {
        slot = atoi(atts[i+1]);
      } else if(strcasecmp(atts[i],"type")==0) {
        slotType==string(atts[i+1]);
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
      } else if(strcasecmp(atts[i],"package")==0) {
        package = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"strip")==0) {
        strip = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"wire")==0) {
        wire = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"plane")==0) {
        plane = string(atts[i+1]);
      } else if(strcasecmp(atts[i],"bar")==0) {
        bar = string(atts[i+1]);
      }
    }

    // fill maps
    string s;
    std::transform(Detector.begin(), Detector.end(), detector.begin(), (int(*)(int)) tolower);
    if(detector=="fcal") {
      s = row + ":" + column;
      fcalMap[s] = make_tuple<int,int,int>(crate,slot,channel);

    } else if(detector=="bcal") {
      s = module + ":" + sector + ":" + layer + ":" + end + ":" + chan;
      bcalMap[s] = make_tuple<int,int,int>(crate,slot,channel);
      
    } else if(detector=="cdc") {
      s = ring + ":" + straw;
      cdcMap[s] = make_tuple<int,int,int>(crate,slot,channel);
      
    } else if(detector=="sc") {
      scMap[sector] = make_tuple<int,int,int>(crate,slot,channel);
    
    } else if(detector=="hodoscope") {
      hodoscopeMap[chan] = make_tuple<int,int,int>(crate,slot,channel);
      
    } else if(detector=="fdccathode") {
      s = package + ":" + layer + ":" + strip;
      fdcCathodeMap[s] = make_tuple<int,int,int>(crate,slot,channel);
      
    } else if(detector=="fdcanode") {
      s = package + ":" + layer + ":" + wire;
      fdcAnodeMap[s] = make_tuple<int,int,int>(crate,slot,channel);
      
    } else if(detector=="tof") {
      s = plane + ":" + bar + ":" + end;
      tofMap[s] = make_tuple<int,int,int>(crate,slot,channel);
     
    } else if(detector=="microscope") {
      microscopeMap[chan] =  make_tuple<int,int,int>(crate,slot,channel);

    } else {
      jerr << endl << endl << "?startElement...unknown detector " << Detector << endl << endl;
    }


  } else {
    jerr << endl << endl << "?startElement...unknown xml tag " << xmlname << endl << endl;
  }

}
    

//--------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DTOFRawHitTranslationADC(const DTOFRawHit* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DTOFRawHitTranslationTDC(const DTOFRawHit* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DBCALHitTranslationADC(const DBCALHit *hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DBCALHitTranslationTDC(const DBCALHit *hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DFCALHitTranslationADC(const DFCALHit* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DFDCHitTranslation(const DFDCHit* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DCDCHitTranslationADC(const DCDCHit* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DSCHitTranslationADC(const DSCHit* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DSCHitTranslationTDC(const DSCHit* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DTaggerTranslationADC(const DTagger* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


cscVal JEventProcessor_rawevent::DTaggerTranslationTDC(const DTagger* hit) {
  return(make_tuple(1,2,3));
}


//----------------------------------------------------------------------------


