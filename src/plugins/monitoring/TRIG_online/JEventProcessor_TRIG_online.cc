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

#include "BCAL/DBCALPoint.h"
#include "FCAL/DFCALHit.h"
#include "TRIGGER/DL1Trigger.h"
#include <DANA/DStatusBits.h>
#include "FCAL/DFCALGeometry.h"


using namespace std;
using namespace jana;

#include <TDirectory.h>
#include <TH1.h>

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

	timing = 0;
	gPARMS->SetDefaultParameter("TRIG_ONLINE:TIMING", timing, "Fill trigger timing histograms: default = false");
	
	// create root folder for trig and cd to it, store main dir
	TDirectory *main = gDirectory;
	gDirectory->mkdir("trig")->cd();
	
        int const nbins=100;

	h1trig_trgbits = new TH1I("h1trig_trgbits", "Trig Trgbits",150,0,150);
	h1trig_trgbits->SetXTitle("trig_mask || (128+fp_trig_mask/256)");
	h1trig_trgbits->SetYTitle("counts");
	//h2trig_fcalVSbcal= new TH2I("h2trig_fcalVSbcal", "E fcal vs E bcal (GeV); Bcal Energy (GeV); Fcal Energy (GeV)",nbins,0,1,nbins,0,2);
	h2trig_fcalVSbcal= new TH2I("h2trig_fcalVSbcal", "E fcal vs E bcal (GeV); Bcal Energy (GeV); Fcal Energy (GeV)",nbins,0,15,nbins,0,15);

	// include timing histograms only if flag set
	if(timing) {	
		h2trig_tfcalVStbcal = new TH2I("h2trig_tfcalVStbcal", "T fcal vs T bcal (ns); Bcal time (ns); Fcal time (ns)",nbins,-200,200,nbins,-200,200);
		h2trig_tfcalVSfcal = new TH2I("h2trig_tfcalVSfcal", "T fcal vs E fcal; Fcal energy (GeV); Fcal time (ns)",nbins,0,2,nbins,-200,200);
		h2trig_tbcalVSbcal = new TH2I("h2trig_tbcalVSbcal", "T bcal vs E bcal; Bcal energ (GeV); Bcal time (ns)",nbins,0,2,nbins,-200,200);
	}

	// monitor some trigger bits separately
	dTrigBits.push_back(1);   // FCAL-BCAL
	dTrigBits.push_back(3);   // FCAL-BCAL && FCAL 
	dTrigBits.push_back(4);   // BCAL only 
	dTrigBits.push_back(5);   // FCAL-BCAL && FCAL
	dTrigBits.push_back(7);   // FCAL-BCAL && FCAL && BCAL
	dTrigBits.push_back(32);  //
	dTrigBits.push_back(33);  //
	dTrigBits.push_back(37);  //
	dTrigBits.push_back(64);  //
	dTrigBits.push_back(65);  //
	dTrigBits.push_back(97);  //
	dTrigBits.push_back(101);  //

	for(size_t loc_i = 0; loc_i < dTrigBits.size(); ++loc_i) {
               
		h2trigbits_fcalVSbcal[dTrigBits[loc_i]] = new TH2I(Form("h2trigbit%d_fcalVSbcal", dTrigBits[loc_i]), Form("Trig %d: E fcal vs E bcal (GeV); Bcal Energy (GeV); Fcal Energy (GeV)",dTrigBits[loc_i]),nbins,0,2,nbins,0,4);

		// include timing histograms only if flag set
		if(timing) {
			h2trigbits_tfcalVStbcal[dTrigBits[loc_i]] = new TH2I(Form("h2trigbit%d_tfcalVStbcal", dTrigBits[loc_i]), Form("Trig %d: T fcal vs T bcal (ns); Bcal time (ns); Fcal time (ns)", dTrigBits[loc_i]),nbins,-200,200,nbins,-200,200);
			h2trigbits_tfcalVSfcal[dTrigBits[loc_i]] = new TH2I(Form("h2trigbit%d_tfcalVSfcal", dTrigBits[loc_i]), Form("Trig %d: T fcal vs E fcal; Fcal energy (GeV); Fcal time (ns)", dTrigBits[loc_i]),nbins,0,2,nbins,-200,200);
			h2trigbits_tbcalVSbcal[dTrigBits[loc_i]] = new TH2I(Form("h2trigbit%d_tbcalVSbcal", dTrigBits[loc_i]), Form("Trig %d: T bcal vs E bcal; Bcal energy (GeV); Bcal time (ns)", dTrigBits[loc_i]),nbins,0,2,nbins,-200,200);
		}
	}

	// back to main dir
	main->cd();

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

	vector<const DBCALPoint*> bcalpoints;
	vector<const DFCALHit*> fcalhits;
	locEventLoop->Get(bcalpoints);
	locEventLoop->Get(fcalhits);
	DFCALGeometry fcalgeom;

	bool isPhysics = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_PHYSICS_EVENT);
	if(!isPhysics)
	  return NOERROR;

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

	int trig_bits = fp_trig_mask > 0? 128 + fp_trig_mask/256: trig_mask;
	//if (fp_trig_mask>0) printf (" Event=%d trig_bits=%d trig_mask=%X fp_trig_mask=%X\n",(int)locEventNumber,trig_bits,trig_mask,fp_trig_mask);

	/* fp_trig_mask & 0x100 - upstream LED
	   fp_trig_mask & 0x200 - downstream LED
	   trig_mask & 0x1 - cosmic trigger*/

        // Compute total energy sums for fcal and bcal (as in the trigger)
        
	// loop over all points in FCAL
	float fcal_energy = 0;
	float fcal_time = 0;
	double fcal_rings_masked = 2;
	if(locEventLoop->GetJEvent().GetRunNumber() < 11127) // ugly run dependence for now
		fcal_rings_masked = 4;
	float rmin = fcal_rings_masked*4*sqrt(2);    // N rings x 4 cm  on the diagonal.
	for (unsigned int jj=0; jj<fcalhits.size(); jj++) {
	  int rowhit = fcalhits[jj]->row;
	  int columnhit = fcalhits[jj]->column;
	  // printf (" Event=%d, jj=%d, rowhit=%d, columnhit=%d\n",(int)locEventNumber,jj,rowhit,columnhit);
	  DVector2 pos = fcalgeom.positionOnFace(rowhit,columnhit);
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

	// FILL HISTOGRAMS
	// Since we are filling histograms local to this plugin, it will not interfere with other ROOT operations: can use plugin-wide ROOT fill lock
	japp->RootFillLock(this); //ACQUIRE ROOT FILL LOCK

	h1trig_trgbits->Fill(trig_bits);
        h2trig_fcalVSbcal->Fill(bcal_energy,fcal_energy);
	if(timing) {
		h2trig_tfcalVStbcal->Fill(bcal_time,fcal_time);
		h2trig_tfcalVSfcal->Fill(fcal_energy,fcal_time);
		h2trig_tbcalVSbcal->Fill(bcal_energy,bcal_time);
	}

	for(size_t loc_i = 0; loc_i < dTrigBits.size(); ++loc_i) {
		if(trig_bits != (int)dTrigBits[loc_i])
			continue;

		h2trigbits_fcalVSbcal[dTrigBits[loc_i]]->Fill(bcal_energy,fcal_energy);
		if(timing) {
			h2trigbits_tfcalVStbcal[dTrigBits[loc_i]]->Fill(bcal_time,fcal_time);
			h2trigbits_tfcalVSfcal[dTrigBits[loc_i]]->Fill(fcal_energy,fcal_time);
			h2trigbits_tbcalVSbcal[dTrigBits[loc_i]]->Fill(bcal_energy,bcal_time);
		}
	}

	japp->RootFillUnLock(this); //RELEASE ROOT FILL LOCK

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
