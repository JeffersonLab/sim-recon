// $Id$
//
//    File: JEventProcessor_event_size.cc
// Created: Tue Jun  7 10:07:36 EDT 2011
// Creator: davidl (on Darwin eleanor.jlab.org 10.7.0 i386)
//

#include "JEventProcessor_event_size.h"
using namespace jana;

// From level1_trigger plugin
#include <level1_trigger/DTrigger.h>

#include <PID/DBeamPhoton.h>
#include <BCAL/DBCALHit.h>
#include <FCAL/DFCALHit.h>
#include <CCAL/DCCALHit.h>
#include <CDC/DCDCHit.h>
#include <FDC/DFDCHit.h>
#include <TOF/DTOFRawHit.h>
#include <START_COUNTER/DSCHit.h>
#include <TAGGER/DTagger.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_event_size());
}
} // "C"


//------------------
// JEventProcessor_event_size (Constructor)
//------------------
JEventProcessor_event_size::JEventProcessor_event_size()
{
	pthread_mutex_init(&mutex, NULL);
}

//------------------
// ~JEventProcessor_event_size (Destructor)
//------------------
JEventProcessor_event_size::~JEventProcessor_event_size()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_event_size::init(void)
{
	evt_tree = new TTree("event","Event Size info");
	evt = new Event();
	evt_tree->Branch("A", &evt);
	
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_event_size::brun(JEventLoop *eventLoop, int runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_event_size::evnt(JEventLoop *loop, int eventnumber)
{
	// Ignore events not passing Level 1 trigger
	const DTrigger *trig=NULL;
	try{
		loop->GetSingle(trig);
		if(!trig->L1fired)return NOERROR;
	}catch(...){
		return NOERROR;
	};

	// Get hit data objects
	vector<const DBeamPhoton*> beamphotons;
	vector<const DBCALHit*> bcalhits;
	vector<const DFCALHit*> fcalhits;
	vector<const DCCALHit*> ccalhits;
	vector<const DCDCHit*> cdchits;
	vector<const DFDCHit*> fdchits;
	vector<const DTOFRawHit*> tofhits;
	vector<const DSCHit*> schits;
	vector<const DTagger*> taggerhits;
	
	loop->Get(beamphotons);
	loop->Get(bcalhits);
	loop->Get(fcalhits);
	loop->Get(ccalhits);
	loop->Get(cdchits);
	loop->Get(fdchits);
	loop->Get(tofhits);
	loop->Get(schits);
	loop->Get(taggerhits);
	
	// Count inner and outer BCAL hits
	unsigned int Nbcalhits_inner = 0;
	unsigned int Nbcalhits_outer = 0;
	for(unsigned int i=0; i<bcalhits.size(); i++){
		if(bcalhits[i]->layer < DBCALGeometry::BCALMID){
			Nbcalhits_inner++;
		}else{
			Nbcalhits_outer++;
		}
	}
	
	// Count FDC wire and FDC strip hits
	unsigned int Nfdchits_anode = 0;
	unsigned int Nfdchits_cathode = 0;
	for(unsigned int i=0; i<fdchits.size(); i++){
		if(fdchits[i]->type == 1){
			Nfdchits_cathode++;
		}else{
			Nfdchits_anode++;
		}
	}

	// Lock mutex while we fill in Event tree
	pthread_mutex_lock(&mutex);
	
	evt->event = eventnumber;
	evt->Egamma = beamphotons.size()>0 ? beamphotons[0]->momentum().Mag():0.0;
	evt->Nbcalhits_inner = Nbcalhits_inner;
	evt->Nbcalhits_outer = Nbcalhits_outer;
	evt->Nfcalhits = fcalhits.size();
	evt->Nccalhits = ccalhits.size();
	evt->Ncdchits = cdchits.size();
	evt->Ncdchits = cdchits.size();
	evt->Nfdchits_anode = Nfdchits_anode;
	evt->Nfdchits_cathode = Nfdchits_cathode;
	evt->Ntofhits = tofhits.size();
	evt->Nschits = schits.size();
	evt->Ntaggerhits = taggerhits.size();

	evt->Ndigitized_values = 0;
		evt->Ndigitized_values += 3*evt->Nbcalhits_inner;    // time and energy from fADC plus time from TDC
		evt->Ndigitized_values += 2*evt->Nbcalhits_outer;    // time and energy from fADC only
		evt->Ndigitized_values += 2*evt->Nfcalhits;          // fADC only
		evt->Ndigitized_values += 2*evt->Nccalhits;          // fADC only
		evt->Ndigitized_values += 2*evt->Ncdchits;           // fADC only
		evt->Ndigitized_values += 1*evt->Nfdchits_anode;     // TDC only
		evt->Ndigitized_values += 2*evt->Nfdchits_cathode;   // fADC only
		evt->Ndigitized_values += 3*evt->Ntofhits;           // fADC plus TDC
		evt->Ndigitized_values += 3*evt->Nschits;            // fADC plus TDC
		evt->Ndigitized_values += 3*evt->Ntaggerhits;        // fADC plus TDC

	// Fill event tree
	evt_tree->Fill();
	
	// Unlock mutex
	pthread_mutex_unlock(&mutex);
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_event_size::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_event_size::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

