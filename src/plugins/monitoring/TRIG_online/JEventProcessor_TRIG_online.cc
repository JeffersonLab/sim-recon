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

#include "JEventProcessor_TRIG_online.h"
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

using namespace std;
using namespace jana;

// #include "TRIG/DTRIG.h"

#include <TDirectory.h>
#include <TH1.h>


// root hist pointers

     static TH1I* h1trig_trgbits = NULL; 
     static TH1I* h1trig_fcal = NULL;
     static TH1I* h1trig_fcalN = NULL;
     static TH1I* h1trig_bcal = NULL;
     static TH1I* h1trig_bcalN = NULL;
     static TH1I* h1trig_tot = NULL;
     static TH2I* h2trig_fcalVSbcal = NULL;


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
	if(h1trig_fcal != NULL){
		japp->RootUnLock();
		return NOERROR;
	}

  // create root folder for trig and cd to it, store main dir
  TDirectory *main = gDirectory;
  gDirectory->mkdir("trig")->cd();


  // book hist
        int const nbins=100;
	
	h1trig_trgbits = new TH1I("h1trig_trgbits", "Trig Trgbits",nbins,0,100);
	h1trig_trgbits->SetXTitle("trig_mask + (50+fp_trig_mask)");
	h1trig_trgbits->SetYTitle("counts");

	h1trig_fcal = new TH1I("h1trig_fcal", "Trig Fcal energy (GeV)",nbins,0,4);
	h1trig_fcal->SetXTitle("Fcal sum energy (GeV)");
	h1trig_fcal->SetYTitle("counts");
	h1trig_fcalN = new TH1I("h1trig_fcalN", "Trig FcalN hits",nbins,0,100);
	h1trig_fcalN->SetXTitle("FcalN hits");
	h1trig_fcalN->SetYTitle("counts");
	h1trig_bcal = new TH1I("h1trig_bcal", "Trig Bcal energy (GeV)",nbins,0,4);
	h1trig_bcal->SetXTitle("Bcal sum energy (GeV)");
	h1trig_bcal->SetYTitle("counts");
	h1trig_bcalN = new TH1I("h1trig_bcalN", "Trig BcalN hits",nbins,0,100);
	h1trig_bcalN->SetXTitle("BcalN hits");
	h1trig_bcalN->SetYTitle("counts");

	h1trig_tot = new TH1I("h1trig_tot", "Trig Tot energy (GeV)",nbins,0,4);
	h1trig_tot->SetXTitle("Total energy (GeV)");
	h1trig_tot->SetYTitle("counts");
	h2trig_fcalVSbcal= new TH2I("h2trig_fcalVSbcal", "E fcal vs E bcal (GeV)",nbins,0,4,nbins,0,4);
	h2trig_fcalVSbcal->SetXTitle("Bcal Energy (GeV)");
	h2trig_fcalVSbcal->SetYTitle("Fcal Energy (GeV)");


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

	japp->RootWriteLock();

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

	int trig_bits = fp_trig_mask > 0? trig_mask + 50 + fp_trig_mask: trig_mask;
	printf (" Event=%d trig_bits=%d trig_mask=%X fp_trig_mask=%X\n",locEventNumber,trig_bits,trig_mask,fp_trig_mask);

	/* fp_trig_mask & 0x100 - upstream LED
	   fp_trig_mask & 0x200 - downstream LED
	   trig_mask & 0x1 - cosmic trigger*/

        // Compute total energy sums for fcal and bcal (as in the trigger)
        
	// loop over all points in FCAL

	float fcal_energy = 0;
	for (unsigned int jj=0; jj<fcalhits.size(); jj++) {
	    fcal_energy += fcalhits[jj]->E;
	    }

	// loop over all points in BCAL

	float bcal_energy = 0;
	for (unsigned int jj=0; jj<bcalpoints.size(); jj++) {
	    bcal_energy += bcalpoints[jj]->E();
	    }

	h1trig_trgbits->Fill(trig_bits);
        h1trig_fcal->Fill(fcal_energy);
        h1trig_fcalN->Fill(fcalhits.size());
        h1trig_bcal->Fill(bcal_energy);
        h1trig_bcalN->Fill(bcalpoints.size());
        h2trig_fcalVSbcal->Fill(bcal_energy,fcal_energy);
        h1trig_tot->Fill(bcal_energy+fcal_energy);
        

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
