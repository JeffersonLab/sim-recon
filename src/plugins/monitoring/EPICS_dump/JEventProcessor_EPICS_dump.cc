// $Id$
//
//    File: JEventProcessor_FCAL_online.cc
// Created: Fri Nov  9 11:58:09 EST 2012
// Creator: wolin (on Linux stan.jlab.org 2.6.32-279.11.1.el6.x86_64 x86_64)


#include <stdint.h>
#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>

#include "JEventProcessor_EPICS_dump.h"
#include <JANA/JApplication.h>

#include "DLorentzVector.h"
#include "TMatrixD.h"

#include "BCAL/DBCALShower.h"
#include "BCAL/DBCALCluster.h"
#include "BCAL/DBCALPoint.h"
#include "BCAL/DBCALHit.h"
#include "FCAL/DFCALShower.h"
#include "FCAL/DFCALCluster.h"
#include "FCAL/DFCALHit.h"
#include "TRACKING/DMCThrown.h"
#include "ANALYSIS/DAnalysisUtilities.h"
#include "TRIGGER/DL1Trigger.h"
#include <DANA/DStatusBits.h>
#include "FCAL/DFCALGeometry.h"


using namespace std;
using namespace jana;

// #include "TRIG/DTRIG.h"

#include <TDirectory.h>
#include <TH1.h>


// root hist pointers

     static TH1I* h1trig_trgbits = NULL;


//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *locApplication){
    InitJANAPlugin(locApplication);
    locApplication->AddProcessor(new JEventProcessor_TRIG_online());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_TRIG_online::JEventProcessor_TRIG_online() {
}


//----------------------------------------------------------------------------------


JEventProcessor_TRIG_online::~JEventProcessor_TRIG_online() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_TRIG_online::init(void) {

  // lock all root operations
  japp->RootWriteLock(); //ACQUIRE ROOT LOCK!!


 // First thread to get here makes all histograms. If one pointer is
 // already not NULL, assume all histograms are defined and return now
	if(h1trig_trgbits != NULL){
		japp->RootUnLock();
		return NOERROR;
	}

  // create root folder for trig and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("epics")->cd();


  // book hist
        int const nbins=20;
	
	h1trig_trgbits = new TH1I("h1trig_trgbits", "Trig Trgbits",nbins,0,20);
	h1trig_trgbits->SetXTitle("trig_mask || (10+fp_trig_mask)");
	h1trig_trgbits->SetYTitle("counts");

  // back to main dir
  main->cd();


  // unlock
  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TRIG_online::brun(jana::JEventLoop* locEventLoop, int locRunNumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TRIG_online::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber) {
  // This is called for every event. Use of common resources like writing
  // to a file or filling a histogram should be mutex protected. Using
  // loop-Get(...) to get reconstructed objects (and thereby activating the
  // reconstruction algorithm) should be done outside of any mutex lock
  // since multiple threads may call this method at the same time.


	vector<const DFCALShower*> locFCALShowers;
	vector<const DBCALPoint*> bcalpoints;
	vector<const DFCALHit*> fcalhits;
	vector<const DFCALCluster*> locFCALClusters;
	//const DDetectorMatches* locDetectorMatches = NULL;
	//locEventLoop->GetSingle(locDetectorMatches);
	locEventLoop->Get(locFCALShowers);
	locEventLoop->Get(bcalpoints);
	locEventLoop->Get(fcalhits);
	locEventLoop->Get(locFCALClusters);
	DFCALGeometry fcalgeom;

	bool isPhysics = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_PHYSICS_EVENT);
	bool isEPICS = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_EPICS_EVENT);

	japp->RootWriteLock();

	if (isPhysics) {
	  // first get trigger bits

	  const DL1Trigger *trig_words = NULL;
	  uint32_t trig_mask, fp_trig_mask;
	  try {
	    locEventLoop->GetSingle(trig_words);
	  } catch(...) {};
	  if (trig_words) {
	    trig_mask = trig_words->trig_mask;
	    fp_trig_mask = trig_words->fp_trig_mask;
	  }
	  else {
	    trig_mask = 0;
	    fp_trig_mask = 0;
	  }

	  int trig_bits = fp_trig_mask > 0? 10 + fp_trig_mask: trig_mask;
	  // printf (" Event=%d trig_bits=%d trig_mask=%X fp_trig_mask=%X\n",(int)locEventNumber,trig_bits,trig_mask,fp_trig_mask);

	  /* fp_trig_mask & 0x100 - upstream LED
	   fp_trig_mask & 0x200 - downstream LED
	   trig_mask & 0x1 - cosmic trigger*/

	  h1trig_trgbits->Fill(trig_bits);
	}
	if (isEPICS) {
	  // process EPICS records
	  printf (" Event=%d is an EPICS record\n",(int)locEventNumber);
	}

        //UnlockState();	
	japp->RootUnLock();

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TRIG_online::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_TRIG_online::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
