
#include <stdint.h>
#include <vector>

#include "JEventProcessor_TPOL_online.h"
#include <JANA/JApplication.h>

using namespace std;
using namespace jana;

#include "TTAB/DTTabUtilities.h"
#include "DAQ/Df250WindowRawData.h"
#include <TPOL/DTPOLSectorDigiHit.h>
#include <TPOL/DTPOLHit.h>
#include <TPOL/DTPOLHit_factory.h>

#include <TDirectory.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

// root hist pointers
static TH1I *tpol_num_events;
static TH1I *tpol_hitMultiplicity;
// static TH2F *hnHitBySector;
//
//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_TPOL_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_TPOL_online::JEventProcessor_TPOL_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_TPOL_online::~JEventProcessor_TPOL_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_TPOL_online::init(void) {
  
  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!
  
  // create root folder for TPOL and cd to it, store main dir
  TDirectory *mainDir = gDirectory;
  TDirectory *tpolDir = gDirectory->mkdir("TPOL");
  tpolDir->cd();

  // book hists
  tpol_num_events = new TH1I("tpol_num_events","TPOL number of events",1, 0.5, 1.5);
  tpol_hitMultiplicity = new TH1I("tpol_hitMultiplicity",";multiplicity of TPOL hits/event;",DTPOLHit_factory::NSECTORS,0.5,DTPOLHit_factory::NSECTORS+0.5);

  // hnHitBySector = new TH2F("hnHitBySector","TPOL hit multiplicity by sector;sector;# hits",DTPOLHit_factory::NSECTORS,-0.5,DTPOLHit_factory::NSECTORS-0.5,10,-0.5,9.5);
  // back to main dir
  mainDir->cd();

  japp->RootUnLock(); //RELEASE ROOT LOCK!!
  
  return NOERROR;
}

//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TPOL_online::brun(JEventLoop *eventLoop, int runnumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TPOL_online::evnt(JEventLoop *eventLoop, int eventnumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.

  // This part is for future use, for the spring 2015 data
  // we cannot get the TPOL hits directly due to the inverted
  // polarity of the detector. 
  // vector<const DTPOLHit*> hits;
  // eventLoop->Get(hits);
  // vector<const DTPOLHit_factory*> sectordigihits;
  // eventLoop->Get(sectordigihits);

  // For the spring 2015 data, we need to grab the Df250WindowRawData objects
  // directly due to the inverted polarity.

  // Each detector's hits
  vector<const Df250WindowRawData*> windowrawdata;
    
  // Get Df250PulseIntegral objects
  eventLoop->Get(windowrawdata);

  // fill hists
  japp->RootWriteLock();
  {

    // Fill number of TPOL events for this file.
    Bool_t hadTPOL = false;
    Int_t hitMultiplicity = 0;
    for(unsigned int i=0;i<windowrawdata.size();i++){
      if(windowrawdata[i]->rocid == 84
	 &&
	 (windowrawdata[i]->slot == 13 || windowrawdata[i]->slot == 14)
	 ){

	// If a waveform is found with rocid 84 and slot 13 or 14,
	// that came from the TPOL
	hadTPOL = true;
	hitMultiplicity++;

	// Get slot and channel of hit	
	// UInt_t sl = windowrawdata[i]->slot;
	// UInt_t ch = windowrawdata[i]->channel;
      }
    }

    if(hadTPOL){
      tpol_num_events->Fill(1);
      tpol_hitMultiplicity->Fill(hitMultiplicity);
    }

    // if((sectordigihits.size()>0))
    //   tpol_num_events->Fill(1);

    // UInt_t hits[DTPOLHit_factory::NSECTORS];
    // // initialize counts for each sector
    // for(Int_t sector=0;sector<DTPOLHit_factory::NSECTORS;sector++){
    //   hits[sector] = 0;
    // }

    // // count hits for each sector
    // for(Int_t n=0;n<sectordigihits.size();n++){
    //   Int_t sector = sectordigihits[n]->sector;
    //   hits[sector-1]++;
    // }

    // // fill hit numbers for each sector
    // for(Int_t sector=1;sector<=DTPOLHit_factory::NSECTORS;sector++){
    //   hnHitBySector->Fill(sector,hits[sector-1]);
    // }

  } // end of lock
  japp->RootUnLock();
    
  return NOERROR;
}
//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TPOL_online::erun(void) {
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
  }


  //----------------------------------------------------------------------------------


  jerror_t JEventProcessor_TPOL_online::fini(void) {
    // Called before program exit after event processing is finished.
    return NOERROR;
  }


  //----------------------------------------------------------------------------------
  //----------------------------------------------------------------------------------
