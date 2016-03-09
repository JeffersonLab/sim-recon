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
#include <DANA/DStatusBits.h>
#include "FCAL/DFCALGeometry.h"


using namespace std;
using namespace jana;

// #include "TRIG/DTRIG.h"

#include <TDirectory.h>
#include <TH1.h>


// root hist pointers

     static TH1I* h1trig_epics = NULL; 
     static TH1I* h1trig_trgbits = NULL; 
     static TH1I* h1trig_fcal = NULL;
     static TH1I* h1trig_fcalN = NULL;
     static TH1I* h1trig_bcal = NULL;
     static TH1I* h1trig_bcalN = NULL;
     static TH1I* h1trig_tot = NULL;
     static TH2I* h2trig_fcalVSbcal = NULL;
     static TH1I* h1trig_fcal_time = NULL;
     static TH1I* h1trig_bcal_time = NULL;
     static TH2I* h2trig_tfcalVStbcal = NULL;
     static TH2I* h2trig_tfcalVSfcal = NULL;
     static TH2I* h2trig_tbcalVSbcal = NULL;
 
     static TH1I* h1trig1_fcal = NULL;
     static TH1I* h1trig1_fcalN = NULL;
     static TH1I* h1trig1_bcal = NULL;
     static TH1I* h1trig1_bcalN = NULL;
     static TH2I* h2trig1_fcalVSbcal = NULL;
     static TH1I* h1trig1_fcal_time = NULL;
     static TH1I* h1trig1_bcal_time = NULL;
     static TH2I* h2trig1_tfcalVStbcal = NULL;
     static TH2I* h2trig1_tfcalVSfcal = NULL;
     static TH2I* h2trig1_tbcalVSbcal = NULL;
 
     static TH1I* h1trig3_fcal = NULL;
     static TH1I* h1trig3_fcalN = NULL;
     static TH1I* h1trig3_bcal = NULL;
     static TH1I* h1trig3_bcalN = NULL;
     static TH2I* h2trig3_fcalVSbcal = NULL;
     static TH1I* h1trig3_fcal_time = NULL;
     static TH1I* h1trig3_bcal_time = NULL;
     static TH2I* h2trig3_tfcalVStbcal = NULL;
     static TH2I* h2trig3_tfcalVSfcal = NULL;
     static TH2I* h2trig3_tbcalVSbcal = NULL;
 
     static TH1I* h1trig5_fcal = NULL;
     static TH1I* h1trig5_fcalN = NULL;
     static TH1I* h1trig5_bcal = NULL;
     static TH1I* h1trig5_bcalN = NULL;
     static TH2I* h2trig5_fcalVSbcal = NULL;
     static TH1I* h1trig5_fcal_time = NULL;
     static TH1I* h1trig5_bcal_time = NULL;
     static TH2I* h2trig5_tfcalVStbcal = NULL;
     static TH2I* h2trig5_tfcalVSfcal = NULL;
     static TH2I* h2trig5_tbcalVSbcal = NULL;
 
     static TH1I* h1trig7_fcal = NULL;
     static TH1I* h1trig7_fcalN = NULL;
     static TH1I* h1trig7_bcal = NULL;
     static TH1I* h1trig7_bcalN = NULL;
     static TH2I* h2trig7_fcalVSbcal = NULL;
     static TH1I* h1trig7_fcal_time = NULL;
     static TH1I* h1trig7_bcal_time = NULL;
     static TH2I* h2trig7_tfcalVStbcal = NULL;
     static TH2I* h2trig7_tfcalVSfcal = NULL;
     static TH2I* h2trig7_tbcalVSbcal = NULL;


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
	
	h1trig_epics = new TH1I("h1trig_epics", "Epics triggers",20,0,20);
	h1trig_epics->SetXTitle("Epics triggers");
	h1trig_epics->SetYTitle("counts");

	h1trig_trgbits = new TH1I("h1trig_trgbits", "Trig Trgbits",30,0,30);
	h1trig_trgbits->SetXTitle("trig_mask || (20+fp_trig_mask/256)");
	h1trig_trgbits->SetYTitle("counts");

	h1trig_fcal = new TH1I("h1trig_fcal", "Trig Fcal energy (GeV)",nbins,0,2);
	h1trig_fcal->SetXTitle("Fcal sum energy (GeV)");
	h1trig_fcal->SetYTitle("counts");
	h1trig_fcalN = new TH1I("h1trig_fcalN", "Trig FcalN hits",nbins,0,100);
	h1trig_fcalN->SetXTitle("FcalN hits");
	h1trig_fcalN->SetYTitle("counts");
	h1trig_bcal = new TH1I("h1trig_bcal", "Trig Bcal energy (GeV)",nbins,0,1);
	h1trig_bcal->SetXTitle("Bcal sum energy (GeV)");
	h1trig_bcal->SetYTitle("counts");
	h1trig_bcalN = new TH1I("h1trig_bcalN", "Trig BcalN hits",nbins,0,100);
	h1trig_bcalN->SetXTitle("BcalN hits");
	h1trig_bcalN->SetYTitle("counts");

	h1trig_tot = new TH1I("h1trig_tot", "Trig Tot energy (GeV)",nbins,0,2);
	h1trig_tot->SetXTitle("Total energy (GeV)");
	h1trig_tot->SetYTitle("counts");
	h2trig_fcalVSbcal= new TH2I("h2trig_fcalVSbcal", "E fcal vs E bcal (GeV)",nbins,0,1,nbins,0,2);
	h2trig_fcalVSbcal->SetXTitle("Bcal Energy (GeV)");
	h2trig_fcalVSbcal->SetYTitle("Fcal Energy (GeV)");
	h1trig_fcal_time = new TH1I("h1trig_fcal_time", "Trig Fcal_Time (ns)",nbins,-200,200);
	h1trig_fcal_time->SetXTitle("Fcal time (ns)");
	h1trig_fcal_time->SetYTitle("counts");
	h1trig_bcal_time = new TH1I("h1trig_bcal_time", "Trig Bcal_Time (ns)",nbins,-200,200);
	h1trig_bcal_time->SetXTitle("Bcal time (ns)");
	h1trig_bcal_time->SetYTitle("counts");
	h2trig_tfcalVStbcal= new TH2I("h2trig_tfcalVStbcal", "T fcal vs T bcal (ns)",nbins,-200,200,nbins,-200,200);
	h2trig_tfcalVStbcal->SetXTitle("Bcal time (ns)");
	h2trig_tfcalVStbcal->SetYTitle("Fcal time (ns)");
	h2trig_tfcalVSfcal= new TH2I("h2trig_tfcalVSfcal", "T fcal vs E fcal",nbins,0,2,nbins,-200,200);
	h2trig_tfcalVSfcal->SetXTitle("Fcal energ (GeV)");
	h2trig_tfcalVSfcal->SetYTitle("Fcal time (ns)");
	h2trig_tbcalVSbcal= new TH2I("h2trig_tbcalVSbcal", "T bcal vs E bcal",nbins,0,2,nbins,-200,200);
	h2trig_tbcalVSbcal->SetXTitle("Bcal energ (GeV)");
	h2trig_tbcalVSbcal->SetYTitle("Bcal time (ns)");

	h1trig1_fcal = new TH1I("h1trig1_fcal", "Trig 1 Fcal energy (GeV)",nbins,0,2);
	h1trig1_fcal->SetXTitle("Fcal sum energy (GeV)");
	h1trig1_fcal->SetYTitle("counts");
	h1trig1_fcalN = new TH1I("h1trig1_fcalN", "Trig 1 FcalN hits",nbins,0,100);
	h1trig1_fcalN->SetXTitle("FcalN hits");
	h1trig1_fcalN->SetYTitle("counts");
	h1trig1_bcal = new TH1I("h1trig1_bcal", "Trig 1 Bcal energy (GeV)",nbins,0,1);
	h1trig1_bcal->SetXTitle("Bcal sum energy (GeV)");
	h1trig1_bcal->SetYTitle("counts");
	h1trig1_bcalN = new TH1I("h1trig1_bcalN", "Trig 1 BcalN hits",nbins,0,100);
	h1trig1_bcalN->SetXTitle("BcalN hits");
	h1trig1_bcalN->SetYTitle("counts");

	h2trig1_fcalVSbcal= new TH2I("h2trig1_fcalVSbcal", "Trig 1 E fcal vs E bcal (GeV)",nbins,0,1,nbins,0,2);
	h2trig1_fcalVSbcal->SetXTitle("Bcal Energy (GeV)");
	h2trig1_fcalVSbcal->SetYTitle("Fcal Energy (GeV)");
	h1trig1_fcal_time = new TH1I("h1trig1_fcal_time", "Trig1 Fcal_Time (ns)",nbins,-200,200);
	h1trig1_fcal_time->SetXTitle("Fcal time (ns)");
	h1trig1_fcal_time->SetYTitle("counts");
	h1trig1_bcal_time = new TH1I("h1trig1_bcal_time", "Trig1 Bcal_Time (ns)",nbins,-200,200);
	h1trig1_bcal_time->SetXTitle("Bcal time (ns)");
	h1trig1_bcal_time->SetYTitle("counts");
	h2trig1_tfcalVStbcal= new TH2I("h2trig1_tfcalVStbcal", "T fcal vs T bcal (ns)",nbins,-200,200,nbins,-200,200);
	h2trig1_tfcalVStbcal->SetXTitle("Bcal time (ns)");
	h2trig1_tfcalVStbcal->SetYTitle("Fcal time (ns)");
	h2trig1_tfcalVSfcal= new TH2I("h2trig1_tfcalVSfcal", "T fcal vs E fcal",nbins,0,2,nbins,-200,200);
	h2trig1_tfcalVSfcal->SetXTitle("Fcal energ (GeV)");
	h2trig1_tfcalVSfcal->SetYTitle("Fcal time (ns)");
	h2trig1_tbcalVSbcal= new TH2I("h2trig1_tbcalVSbcal", "T bcal vs E bcal",nbins,0,2,nbins,-200,200);
	h2trig1_tbcalVSbcal->SetXTitle("Bcal energ (GeV)");
	h2trig1_tbcalVSbcal->SetYTitle("Bcal time (ns)");

	h1trig3_fcal = new TH1I("h1trig3_fcal", "Trig 3 Fcal energy (GeV)",nbins,0,2);
	h1trig3_fcal->SetXTitle("Fcal sum energy (GeV)");
	h1trig3_fcal->SetYTitle("counts");
	h1trig3_fcalN = new TH1I("h1trig3_fcalN", "Trig 3 FcalN hits",nbins,0,100);
	h1trig3_fcalN->SetXTitle("FcalN hits");
	h1trig3_fcalN->SetYTitle("counts");
	h1trig3_bcal = new TH1I("h1trig3_bcal", "Trig 3 Bcal energy (GeV)",nbins,0,1);
	h1trig3_bcal->SetXTitle("Bcal sum energy (GeV)");
	h1trig3_bcal->SetYTitle("counts");
	h1trig3_bcalN = new TH1I("h1trig3_bcalN", "Trig 3 BcalN hits",nbins,0,100);
	h1trig3_bcalN->SetXTitle("BcalN hits");
	h1trig3_bcalN->SetYTitle("counts");

	h2trig3_fcalVSbcal= new TH2I("h2trig3_fcalVSbcal", "Trig 3 E fcal vs E bcal (GeV)",nbins,0,1,nbins,0,2);
	h2trig3_fcalVSbcal->SetXTitle("Bcal Energy (GeV)");
	h2trig3_fcalVSbcal->SetYTitle("Fcal Energy (GeV)");
	h1trig3_fcal_time = new TH1I("h1trig3_fcal_time", "Trig3 Fcal_Time (ns)",nbins,-200,200);
	h1trig3_fcal_time->SetXTitle("Fcal time (ns)");
	h1trig3_fcal_time->SetYTitle("counts");
	h1trig3_bcal_time = new TH1I("h1trig3_bcal_time", "Trig3 Bcal_Time (ns)",nbins,-200,200);
	h1trig3_bcal_time->SetXTitle("Bcal time (ns)");
	h1trig3_bcal_time->SetYTitle("counts");
	h2trig3_tfcalVStbcal= new TH2I("h2trig3_tfcalVStbcal", "T fcal vs T bcal (ns)",nbins,-200,200,nbins,-200,200);
	h2trig3_tfcalVStbcal->SetXTitle("Bcal time (ns)");
	h2trig3_tfcalVStbcal->SetYTitle("Fcal time (ns)");
	h2trig3_tfcalVSfcal= new TH2I("h2trig3_tfcalVSfcal", "T fcal vs E fcal",nbins,0,2,nbins,-200,200);
	h2trig3_tfcalVSfcal->SetXTitle("Fcal energ (GeV)");
	h2trig3_tfcalVSfcal->SetYTitle("Fcal time (ns)");
	h2trig3_tbcalVSbcal= new TH2I("h2trig3_tbcalVSbcal", "T bcal vs E bcal",nbins,0,2,nbins,-200,200);
	h2trig3_tbcalVSbcal->SetXTitle("Bcal energ (GeV)");
	h2trig3_tbcalVSbcal->SetYTitle("Bcal time (ns)");

	h1trig5_fcal = new TH1I("h1trig5_fcal", "Trig 5 Fcal energy (GeV)",nbins,0,2);
	h1trig5_fcal->SetXTitle("Fcal sum energy (GeV)");
	h1trig5_fcal->SetYTitle("counts");
	h1trig5_fcalN = new TH1I("h1trig5_fcalN", "Trig 5 FcalN hits",nbins,0,100);
	h1trig5_fcalN->SetXTitle("FcalN hits");
	h1trig5_fcalN->SetYTitle("counts");
	h1trig5_bcal = new TH1I("h1trig5_bcal", "Trig 5 Bcal energy (GeV)",nbins,0,1);
	h1trig5_bcal->SetXTitle("Bcal sum energy (GeV)");
	h1trig5_bcal->SetYTitle("counts");
	h1trig5_bcalN = new TH1I("h1trig5_bcalN", "Trig 5 BcalN hits",nbins,0,100);
	h1trig5_bcalN->SetXTitle("BcalN hits");
	h1trig5_bcalN->SetYTitle("counts");

	h2trig5_fcalVSbcal= new TH2I("h2trig5_fcalVSbcal", "Trig 5 E fcal vs E bcal (GeV)",nbins,0,1,nbins,0,2);
	h2trig5_fcalVSbcal->SetXTitle("Bcal Energy (GeV)");
	h2trig5_fcalVSbcal->SetYTitle("Fcal Energy (GeV)");
	h1trig5_fcal_time = new TH1I("h1trig5_fcal_time", "Trig5 Fcal_Time (ns)",nbins,-200,200);
	h1trig5_fcal_time->SetXTitle("Fcal time (ns)");
	h1trig5_fcal_time->SetYTitle("counts");
	h1trig5_bcal_time = new TH1I("h1trig5_bcal_time", "Trig5 Bcal_Time (ns)",nbins,-200,200);
	h1trig5_bcal_time->SetXTitle("Bcal time (ns)");
	h1trig5_bcal_time->SetYTitle("counts");
	h2trig5_tfcalVStbcal= new TH2I("h2trig5_tfcalVStbcal", "T fcal vs T bcal (ns)",nbins,-200,200,nbins,-200,200);
	h2trig5_tfcalVStbcal->SetXTitle("Bcal time (ns)");
	h2trig5_tfcalVStbcal->SetYTitle("Fcal time (ns)");
	h2trig5_tfcalVSfcal= new TH2I("h2trig5_tfcalVSfcal", "T fcal vs E fcal",nbins,0,2,nbins,-200,200);
	h2trig5_tfcalVSfcal->SetXTitle("Fcal energ (GeV)");
	h2trig5_tfcalVSfcal->SetYTitle("Fcal time (ns)");
	h2trig5_tbcalVSbcal= new TH2I("h2trig5_tbcalVSbcal", "T bcal vs E bcal",nbins,0,2,nbins,-200,200);
	h2trig5_tbcalVSbcal->SetXTitle("Bcal energ (GeV)");
	h2trig5_tbcalVSbcal->SetYTitle("Bcal time (ns)");

	h1trig7_fcal = new TH1I("h1trig7_fcal", "Trig 7 Fcal energy (GeV)",nbins,0,2);
	h1trig7_fcal->SetXTitle("Fcal sum energy (GeV)");
	h1trig7_fcal->SetYTitle("counts");
	h1trig7_fcalN = new TH1I("h1trig7_fcalN", "Trig 7 FcalN hits",nbins,0,100);
	h1trig7_fcalN->SetXTitle("FcalN hits");
	h1trig7_fcalN->SetYTitle("counts");
	h1trig7_bcal = new TH1I("h1trig7_bcal", "Trig 7 Bcal energy (GeV)",nbins,0,1);
	h1trig7_bcal->SetXTitle("Bcal sum energy (GeV)");
	h1trig7_bcal->SetYTitle("counts");
	h1trig7_bcalN = new TH1I("h1trig7_bcalN", "Trig 7 BcalN hits",nbins,0,100);
	h1trig7_bcalN->SetXTitle("BcalN hits");
	h1trig7_bcalN->SetYTitle("counts");

	h2trig7_fcalVSbcal= new TH2I("h2trig7_fcalVSbcal", "Trig 7 E fcal vs E bcal (GeV)",nbins,0,1,nbins,0,2);
	h2trig7_fcalVSbcal->SetXTitle("Bcal Energy (GeV)");
	h2trig7_fcalVSbcal->SetYTitle("Fcal Energy (GeV)");
	h1trig7_fcal_time = new TH1I("h1trig7_fcal_time", "Trig7 Fcal_Time (ns)",nbins,-200,200);
	h1trig7_fcal_time->SetXTitle("Fcal time (ns)");
	h1trig7_fcal_time->SetYTitle("counts");
	h1trig7_bcal_time = new TH1I("h1trig7_bcal_time", "Trig7 Bcal_Time (ns)",nbins,-200,200);
	h1trig7_bcal_time->SetXTitle("Bcal time (ns)");
	h1trig7_bcal_time->SetYTitle("counts");
	h2trig7_tfcalVStbcal= new TH2I("h2trig7_tfcalVStbcal", "T fcal vs T bcal (ns)",nbins,-200,200,nbins,-200,200);
	h2trig7_tfcalVStbcal->SetXTitle("Bcal time (ns)");
	h2trig7_tfcalVStbcal->SetYTitle("Fcal time (ns)");
	h2trig7_tfcalVSfcal= new TH2I("h2trig7_tfcalVSfcal", "T fcal vs E fcal",nbins,0,2,nbins,-200,200);
	h2trig7_tfcalVSfcal->SetXTitle("Fcal energ (GeV)");
	h2trig7_tfcalVSfcal->SetYTitle("Fcal time (ns)");
	h2trig7_tbcalVSbcal= new TH2I("h2trig7_tbcalVSbcal", "T bcal vs E bcal",nbins,0,2,nbins,-200,200);
	h2trig7_tbcalVSbcal->SetXTitle("Bcal energ (GeV)");
	h2trig7_tbcalVSbcal->SetYTitle("Bcal time (ns)");


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
	if(! isPhysics) {
	  printf ("Non-physics Event=%d\n",(int)locEventNumber);
	  h1trig_epics->Fill(0.);
	  return NOERROR;
	}
	h1trig_epics->Fill(1.);

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

	h1trig_epics->Fill(2.);
	
	int trig_bits = fp_trig_mask > 0? 20 + fp_trig_mask/256: trig_mask;
	if (fp_trig_mask>0) printf (" Event=%d trig_bits=%d trig_mask=%X fp_trig_mask=%X\n",(int)locEventNumber,trig_bits,trig_mask,fp_trig_mask);

	/* fp_trig_mask & 0x100 - upstream LED
	   fp_trig_mask & 0x200 - downstream LED
	   trig_mask & 0x1 - cosmic trigger*/

        // Compute total energy sums for fcal and bcal (as in the trigger)
        
	// loop over all points in FCAL

	float fcal_energy = 0;
	float fcal_time = 0;
	float rmin = 4*4*sqrt(2);    // 4 layers x 4 cm  on the diagonal.
	for (unsigned int jj=0; jj<fcalhits.size(); jj++) {
	  DVector2 pos = fcalgeom.positionOnFace(jj);
	  double r = sqrt(pos.X()*pos.X() + pos.Y()*pos.Y());
	  if (r <= rmin) continue;    // keep only hits that are outside a minimum radius

	     // require trigger threshold in sum
	     // if (fcalhits[jj]->E > 65*0.27*7.5/1000) {
	  if (fcalhits[jj]->E*7.5/fcalhits[jj]->intOverPeak > 1.0*(65*0.27*7.5/1000)) {
	       fcal_energy += fcalhits[jj]->E*7.5/fcalhits[jj]->intOverPeak;
	       fcal_time += fcalhits[jj]->t*fcalhits[jj]->E*7.5/fcalhits[jj]->intOverPeak;    // calculate energy weighted time
	     }
	}
	fcal_time = fcal_energy > 0 ? fcal_time/fcal_energy : -200; 

	// loop over all points in BCAL

	float bcal_energy = 0;
	float bcal_time = 0;
	for (unsigned int jj=0; jj<bcalpoints.size(); jj++) {
	  if (bcalpoints[jj]->E() > (20*0.045*15/2)/1000) { 
	    bcal_energy += bcalpoints[jj]->E();
	    bcal_time += bcalpoints[jj]->t()*bcalpoints[jj]->E();
	  }
        }
	bcal_time = bcal_energy > 0 ? bcal_time/bcal_energy: -200; 

	h1trig_trgbits->Fill(trig_bits);
        h1trig_fcal->Fill(fcal_energy);
        h1trig_fcalN->Fill(fcalhits.size());
        h1trig_bcal->Fill(bcal_energy);
        h1trig_bcalN->Fill(bcalpoints.size());
        h2trig_fcalVSbcal->Fill(bcal_energy,fcal_energy);
        h1trig_tot->Fill(bcal_energy+fcal_energy);
        h1trig_fcal_time->Fill(fcal_time);
        h1trig_bcal_time->Fill(bcal_time);
        h2trig_tfcalVStbcal->Fill(bcal_time,fcal_time);
        h2trig_tfcalVSfcal->Fill(fcal_energy,fcal_time);
        h2trig_tbcalVSbcal->Fill(bcal_energy,bcal_time);

	if (trig_bits == 1) {
	  h1trig1_fcal->Fill(fcal_energy);
	  h1trig1_fcalN->Fill(fcalhits.size());
	  h1trig1_bcal->Fill(bcal_energy);
	  h1trig1_bcalN->Fill(bcalpoints.size());
	  h2trig1_fcalVSbcal->Fill(bcal_energy,fcal_energy);
	  h1trig1_fcal_time->Fill(fcal_time);
	  h1trig1_bcal_time->Fill(bcal_time);
	  h2trig1_tfcalVStbcal->Fill(bcal_time,fcal_time);
          h2trig1_tfcalVSfcal->Fill(fcal_energy,fcal_time);
          h2trig1_tbcalVSbcal->Fill(bcal_energy,bcal_time);
	}
	else if(trig_bits == 3) {
	  h1trig3_fcal->Fill(fcal_energy);
	  h1trig3_fcalN->Fill(fcalhits.size());
	  h1trig3_bcal->Fill(bcal_energy);
	  h1trig3_bcalN->Fill(bcalpoints.size());
	  h2trig3_fcalVSbcal->Fill(bcal_energy,fcal_energy);
	  h1trig3_fcal_time->Fill(fcal_time);
	  h1trig3_bcal_time->Fill(bcal_time);
	  h2trig3_tfcalVStbcal->Fill(bcal_time,fcal_time);
          h2trig3_tfcalVSfcal->Fill(fcal_energy,fcal_time);
          h2trig3_tbcalVSbcal->Fill(bcal_energy,bcal_time);
	}
	else if (trig_bits == 5) {
	  h1trig5_fcal->Fill(fcal_energy);
	  h1trig5_fcalN->Fill(fcalhits.size());
	  h1trig5_bcal->Fill(bcal_energy);
	  h1trig5_bcalN->Fill(bcalpoints.size());
	  h2trig5_fcalVSbcal->Fill(bcal_energy,fcal_energy);
	  h1trig5_fcal_time->Fill(fcal_time);
	  h1trig5_bcal_time->Fill(bcal_time);
	  h2trig5_tfcalVStbcal->Fill(bcal_time,fcal_time);
          h2trig5_tfcalVSfcal->Fill(fcal_energy,fcal_time);
          h2trig5_tbcalVSbcal->Fill(bcal_energy,bcal_time);
	}
	else if (trig_bits == 7) {
	  h1trig7_fcal->Fill(fcal_energy);
	  h1trig7_fcalN->Fill(fcalhits.size());
	  h1trig7_bcal->Fill(bcal_energy);
	  h1trig7_bcalN->Fill(bcalpoints.size());
	  h2trig7_fcalVSbcal->Fill(bcal_energy,fcal_energy);
	  h1trig7_fcal_time->Fill(fcal_time);
	  h1trig7_bcal_time->Fill(bcal_time);
	  h2trig7_tfcalVStbcal->Fill(bcal_time,fcal_time);
          h2trig7_tfcalVSfcal->Fill(fcal_energy,fcal_time);
          h2trig7_tbcalVSbcal->Fill(bcal_energy,bcal_time);
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
