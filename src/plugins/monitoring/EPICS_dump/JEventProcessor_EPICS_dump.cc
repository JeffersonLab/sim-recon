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
#include "DAQ/DEPICSvalue.h"


using namespace std;
using namespace jana;

// #include "TRIG/DTRIG.h"

#include <TDirectory.h>
#include <TH1.h>


// root hist pointers

     static TH1I* h1epics_trgbits = NULL;
     static TH1I* h1epics_AD00 = NULL;
     static TH2I* h2epics_pos_inner = NULL;
     static TH2I* h2epics_pos_outer = NULL;


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
	if(h1epics_trgbits != NULL){
		japp->RootUnLock();
		return NOERROR;
	}

  // create root folder for trig and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("EPICS_dump")->cd();


  // book hist
        int const nbins=100;
	
	h1epics_trgbits = new TH1I("h1epics_trgbits", "Trig Trgbits",30,0,30);
	h1epics_trgbits->SetXTitle("trig_mask || (10+fp_trig_mask)");
	h1epics_trgbits->SetYTitle("counts");
	h1epics_trgbits = new TH1I("h1epics_trgbits", "Trig Trgbits",30,0,30);
	h1epics_trgbits->SetXTitle("trig_mask || (10+fp_trig_mask)");
	h1epics_trgbits->SetYTitle("counts");



	h1epics_AD00 = new TH1I("h1epics_AD00", "Current AD00",nbins,0,500);
	h1epics_AD00->SetXTitle("Current AD00 (nA)");
	h1epics_AD00->SetYTitle("counts");

	h2epics_pos_inner = new TH2I("h1epics_pos_inner", "Position AC inner",nbins,-50,50,nbins,-50,50);
	h2epics_pos_inner->SetXTitle("Position AC inner x (mm)");
	h2epics_pos_inner->SetYTitle("Position AC inner y (mm)");
	h2epics_pos_outer = new TH2I("h1epics_pos_outer", "Position AC outer",nbins,-50,50,nbins,-50,50);
	h2epics_pos_outer->SetXTitle("Position AC outer x (mm)");
	h2epics_pos_outer->SetYTitle("Position AC outer y (mm)");

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
	vector<const DEPICSvalue*> epicsvalues;
	//const DDetectorMatches* locDetectorMatches = NULL;
	//locEventLoop->GetSingle(locDetectorMatches);
	locEventLoop->Get(locFCALShowers);
	locEventLoop->Get(bcalpoints);
	locEventLoop->Get(fcalhits);
	locEventLoop->Get(locFCALClusters);
	locEventLoop->Get(epicsvalues);	
	DFCALGeometry fcalgeom;

	bool isPhysics = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_PHYSICS_EVENT);
	bool isEPICS = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_EPICS_EVENT);

	japp->RootWriteLock();

	if (isPhysics) {
	  // first get trigger bits

	  const DL1Trigger *trig_words = NULL;
	  uint32_t trig_mask, fp_trig_mask;
	  // uint32_t nsync; /* sync event number */
	  // uint32_t trig_number;
	  uint32_t livetime; /* accumulated livetime */
	  uint32_t busytime; /* accumulated busy time */
	  // uint32_t live_inst; /* instantaneous livetime */
	  uint32_t timestamp;   /*unix time*/

	  try {
	    locEventLoop->GetSingle(trig_words);
	  } catch(...) {};
	  if (trig_words) {
	    trig_mask = trig_words->trig_mask;
	    fp_trig_mask = trig_words->fp_trig_mask;
	    livetime = trig_words->livetime;
	    busytime = trig_words->busytime;
	  }
	  else {
	    trig_mask = 0;
	    fp_trig_mask = 0;
	    livetime = 0;
	    busytime = 0;
	  }

	  int trig_bits = fp_trig_mask > 0? 20 + fp_trig_mask/256: trig_mask;
	  if (fp_trig_mask>0) printf (" Event=%d trig_bits=%d trig_mask=%X fp_trig_mask=%X\n",(int)locEventNumber,trig_bits,trig_mask,fp_trig_mask);

	  /* fp_trig_mask & 0x100 - upstream LED
	   fp_trig_mask & 0x200 - downstream LED
	   trig_mask & 0x1 - cosmic trigger*/

	  h1epics_trgbits->Fill(trig_bits);
	}
	else if (isEPICS) {
	  // else if (TEST) {
	  // process EPICS records
	  printf (" Event=%d is an EPICS record\n",(int)locEventNumber);


	  // read in whatever epics values are in this event

	  // save their values
	  float pos_default=-1000.;
	  float xpos_inner = pos_default;
	  float ypos_inner = pos_default;;
	  float xpos_outer = pos_default;
	  float ypos_outer = pos_default;
	  for(vector<const DEPICSvalue*>::const_iterator val_itr = epicsvalues.begin();
	      val_itr != epicsvalues.end(); val_itr++) {
		const DEPICSvalue* epics_val = *val_itr;
		cout << "EPICS:  " << epics_val->name << " = " << epics_val->sval << endl;
		float fconv = atof(epics_val->sval.c_str());
		bool isDigit = epics_val->name.length()> 12 && isdigit(epics_val->name[12]);
		// cout << "isDigit=" << isDigit << " string=" << epics_val->name << endl;
		if ((epics_val->name.substr(0,11) == "BCAL:pulser") & isDigit) {
		  double freq = 1.e8/fconv;  // cover to s: period is in units 10 ns
		    cout << "BCAL:pulser=" << epics_val->name.substr(0,11) << epics_val->fval <<  " freq=" << freq <<endl;
		}     
		else if (epics_val->name == "IBCAD00CRCUR6") {
		  h1epics_AD00->Fill(fconv);
	          cout << "IBCAD00CRCUR6 " << epics_val->name << " fconv=" << fconv << endl; 
		}
		else if (epics_val->name == "AC:inner:position:x") {
		  xpos_inner = fconv;
	        }  
		else if (epics_val->name == "AC:inner:position:y") {
		  ypos_inner = fconv;
	        }  
		else if (epics_val->name == "AC:outer:position:x") {
		  xpos_outer = fconv;
	        }  
		else if (epics_val->name == "AC:outer:position:y") {
		  ypos_outer = fconv;
	        }  
	  }
	  if (xpos_inner> pos_default && ypos_inner > pos_default) { 
	    h2epics_pos_inner->Fill(xpos_inner,ypos_inner);
	  }
	  if (xpos_outer> pos_default && ypos_outer > pos_default) { 
	    h2epics_pos_outer->Fill(xpos_outer,ypos_outer);
	  }

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
