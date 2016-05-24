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
#include "DAQ/DTSscalers.h"


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
     Int_t const nscalers=32;
     static TH1I* h1_trig_rates[nscalers]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                             NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                             NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                             NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
     static TH1I* h1_trig_livetimes[nscalers]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                             NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                             NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,
                                             NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL};
     static TH1I* h1epics_liveinst = NULL;
     static TH1F* h1epics_AD00_VSevent = NULL;
     static TH1F* h1epics_entries1_VSevent = NULL;;
     static TH1F* h1epics_entries2_VSevent = NULL;
     static TH1F* h1epics_liveinst_VSevent = NULL;

     static uint64_t save_ntrig[nscalers];
     static uint64_t save_ntrig0[nscalers];
     static bool init_ntrig = false;


//----------------------------------------------------------------------------------


// Routine used to create our JEventProcessor
extern "C"{
  void InitPlugin(JApplication *locApplication){
    InitJANAPlugin(locApplication);
    locApplication->AddProcessor(new JEventProcessor_EPICS_dump());
  }
}


//----------------------------------------------------------------------------------


JEventProcessor_EPICS_dump::JEventProcessor_EPICS_dump() {
}


//----------------------------------------------------------------------------------


JEventProcessor_EPICS_dump::~JEventProcessor_EPICS_dump() {
}


//----------------------------------------------------------------------------------

jerror_t JEventProcessor_EPICS_dump::init(void) {

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
        Int_t const nscalers=32;
	char string[132];
	
	h1epics_trgbits = new TH1I("h1epics_trgbits", "Trig Trgbits",80,0,80);
	h1epics_trgbits->SetXTitle("trig_mask || (70+fp_trig_mask/256)");
	h1epics_trgbits->SetYTitle("counts");
	h1epics_trgbits = new TH1I("h1epics_trgbits", "Trig Trgbits",80,0,80);
	h1epics_trgbits->SetXTitle("trig_mask || (70+fp_trig_mask/256)");
	h1epics_trgbits->SetYTitle("counts");
    
	for (Int_t j=0; j<nscalers;j++) {
	  sprintf (string,"Rates%d(kHz)",j);
	  h1_trig_rates[j] = new TH1I(string,string,nbins,0,100);
	  h1_trig_rates[j]->SetTitle(string);
	  sprintf (string,"Livetimes%d",j);
	  h1_trig_livetimes[j] = new TH1I(string,string,nbins,0,1);
	  h1_trig_livetimes[j]->SetTitle(string);    
	  }

	h1epics_liveinst = new TH1I("h1epics_liveinst", "Liveinst",nbins,0,1);
	h1epics_liveinst->SetXTitle("Instantaneous Live time");
	h1epics_liveinst->SetYTitle("counts");	

	h1epics_AD00 = new TH1I("h1epics_AD00", "Current AD00",nbins,0,500);
	h1epics_AD00->SetXTitle("Current AD00 (nA)");
	h1epics_AD00->SetYTitle("counts");

	h2epics_pos_inner = new TH2I("h1epics_pos_inner", "Position AC inner",nbins,-10,10,nbins,-10,10);
	h2epics_pos_inner->SetXTitle("Position AC inner x (mm)");
	h2epics_pos_inner->SetYTitle("Position AC inner y (mm)");
	h2epics_pos_outer = new TH2I("h1epics_pos_outer", "Position AC outer",nbins,-20,20,nbins,-20,20);
	h2epics_pos_outer->SetXTitle("Position AC outer x (mm)");
	h2epics_pos_outer->SetYTitle("Position AC outer y (mm)");


        h1epics_AD00_VSevent = new TH1F("h1epics_AD00_VSevent", "Current AD00 (nA)",200,0,2e6);
	h1epics_AD00_VSevent->SetXTitle("Event number");
	h1epics_AD00_VSevent->SetYTitle("Current (nA)");
        h1epics_entries1_VSevent = new TH1F("h1epics_entries1_VSevent", "No synch events/bin",200,0,200e6);
	h1epics_entries1_VSevent->SetXTitle("Event number");
	h1epics_entries1_VSevent->SetYTitle("No of synch events");
        h1epics_entries2_VSevent = new TH1F("h1epics_entries2_VSevent", "No epics events/bin",200,0,2e6);
	h1epics_entries2_VSevent->SetXTitle("Event number");
	h1epics_entries2_VSevent->SetYTitle("No of synch events");
        h1epics_liveinst_VSevent = new TH1F("h1epics_liveinst_VSevent", "",200,0,200e6);
	h1epics_liveinst_VSevent->SetXTitle("Event number");
	h1epics_liveinst_VSevent->SetYTitle("Live Time");



  // back to main dir
  main->cd();

  // initialize live times
  for (Int_t j=0; j<nscalers; j++) {
    save_ntrig[j] = 0;
    save_ntrig0[j] = 0;
    }


  // unlock
  japp->RootUnLock(); //RELEASE ROOT LOCK!!

  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_EPICS_dump::brun(jana::JEventLoop* locEventLoop, int locRunNumber) {
  // This is called whenever the run number changes
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_EPICS_dump::evnt(jana::JEventLoop* locEventLoop, uint64_t locEventNumber) {
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
        Int_t const nscalers=32;

	bool isPhysics = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_PHYSICS_EVENT);
	bool isEPICS = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_EPICS_EVENT);
	// bool isSynch = locEventLoop->GetJEvent().GetStatusBit(kSTATUS_SYNCH_EVENT);

	japp->RootWriteLock();

	uint32_t trig_mask=0, fp_trig_mask=0;
	uint32_t temp_mask=0;

	if (isPhysics) {
	  // first get trigger bits

	  const DL1Trigger *trig_words = NULL;

	  try {
	    locEventLoop->GetSingle(trig_words);
	  } catch(...) {};
	  if (trig_words) {
	    trig_mask = trig_words->trig_mask;
	    fp_trig_mask = trig_words->fp_trig_mask; 
	    for (int j=0; j<nscalers; j++) {
	      temp_mask = trig_mask & 1<<j;
	      if (temp_mask && init_ntrig) save_ntrig[j] += 1;
	    } 
	  }
	  else {
	    trig_mask = 0;
	    fp_trig_mask = 0;
	  }

	  int trig_bits = fp_trig_mask > 0? 70 + fp_trig_mask/256: trig_mask;
	  if (fp_trig_mask>0) printf (" Event=%d trig_bits=%d trig_mask=%X fp_trig_mask=%X\n",(int)locEventNumber,trig_bits,trig_mask,fp_trig_mask);

	  /* fp_trig_mask & 0x100 - upstream LED
	   fp_trig_mask & 0x200 - downstream LED
	   trig_mask & 0x1 - cosmic trigger*/

	  h1epics_trgbits->Fill(trig_bits);


	  // now get scalers

	  const DTSscalers *ts_scalers = NULL;
	  uint32_t livetime; /* accumulated livetime */
	  uint32_t busytime; /* accumulated busy time */
	  uint32_t live_inst; /* instantaneous livetime */
	  uint32_t timestamp;   /*unix time*/

	  uint32_t gtp_sc[nscalers]; /* number of input triggers from GTP for 32 lanes (32 trigger bits) */
	  // uint32_t fp_sc[nscalers1]; /* number of TS front pannel triggers for 16 fron pannel lines (16 trigger bits) */
	  uint32_t gtp_rate[nscalers]; /* instant. rate of GTP triggers */
	  // uint32_t fp_rate[nscalers1]; /* instant. rate of FP triggers */

	  try {
	    locEventLoop->GetSingle(ts_scalers);
	  } catch(...) {};
	  if (ts_scalers) {
	    livetime = ts_scalers->live_time;
	    busytime = ts_scalers->busy_time;
	    live_inst = ts_scalers->inst_livetime;
	    timestamp = ts_scalers->time;
	    h1epics_liveinst->Fill((float)live_inst/1000.);
	    h1epics_liveinst_VSevent->Fill((float)locEventNumber,(float)live_inst/1000.);
	    h1epics_entries1_VSevent->Fill((float)locEventNumber);
	    // printf ("Event=%d livetime=%d busytime=%d time=%d live_inst=%d\n",(int)locEventNumber,livetime,busytime,(int)timestamp,live_inst);
	    if (! init_ntrig) {
	      for (int j=0; j<nscalers; j++) {
		save_ntrig0[j] = ts_scalers->gtp_scalers[j];
		init_ntrig = true;
	      }
	    }
	    for (int j=0; j<nscalers; j++) {
	      gtp_sc[j] = ts_scalers->gtp_scalers[j] - save_ntrig0[j];
	      gtp_rate[j] = ts_scalers->gtp_rate[j];
	      // printf ("TSscalers: Event=%d j=%d gtp_sc=%d gtp_rate=%d\n",(int)locEventNumber,j,gtp_sc[j],gtp_rate[j]);
	      // printf ("TSscalers: Event=%d trig_mask=%X temp_mask=%X save_ntrig=%d\n",(int)locEventNumber,trig_mask,temp_mask,(int)save_ntrig[j]);
	      h1_trig_rates[j]->Fill(gtp_rate[j]/1000);     // plot in kHz
	      if (gtp_sc[j] >0) {
		h1_trig_livetimes[j]->Fill((float)save_ntrig[j]/(float)gtp_sc[j]);
		save_ntrig0[j] = ts_scalers->gtp_scalers[j];
		save_ntrig[j] = 0;
	      }
	    }
	  }
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
		unsigned int iconv = atoi(epics_val->sval.c_str());
		bool isDigit = epics_val->name.length()> 12 && isdigit(epics_val->name[12]);
		// cout << "isDigit=" << isDigit << " string=" << epics_val->name << endl;
		if ((epics_val->name.substr(0,11) == "BCAL:pulser") & isdigit(epics_val->name[11])) {
		  val_itr++;   // increment counter and get low order word
		  epics_val = *val_itr;
		  cout << "EPICS:  " << epics_val->name << " = " << epics_val->sval << endl;
		  unsigned int iconv_low = atoi(epics_val->sval.c_str());
		  // cout << "BCAL:pulser status=" << epics_val->name.substr(0,11) << endl;
		  // printf ("BCAL:pulser iconv=%d, %0X, iconv_low=%d %0X\n",iconv,iconv,iconv_low,iconv_low);
		  iconv_low = iconv_low >> 31;
		  iconv  = ((iconv <<1) & 0XFFFF) + iconv_low;
		  // cout << "BCAL:pulser status=" << epics_val->name.substr(0,11) << endl;
		  printf ("BCAL:pulser status=%d, %0X, iconv_low=%0X\n",iconv,iconv,iconv_low);
		}     
		else if ((epics_val->name.substr(0,11) == "BCAL:pulser") & isDigit) {
		  double freq = 1.e8/fconv;  // cover to s: period is in units 10 ns
		  cout << "BCAL:pulser=" << epics_val->name.substr(0,11) << epics_val->fval <<  " freq=" << freq <<endl;
		  }
		else if (epics_val->name == "IBCAD00CRCUR6") {
		  h1epics_AD00->Fill(fconv);
		  h1epics_AD00_VSevent->Fill((float)locEventNumber,fconv);
		  h1epics_entries2_VSevent->Fill((float)locEventNumber);
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


jerror_t JEventProcessor_EPICS_dump::erun(void) {
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}


//----------------------------------------------------------------------------------


jerror_t JEventProcessor_EPICS_dump::fini(void) {
  // Called before program exit after event processing is finished.
  return NOERROR;
}


//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
