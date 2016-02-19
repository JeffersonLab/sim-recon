// $Id$
//
//    File: JEventProcessor_TAGGER_online.cc
// Created: Thu Feb 18 07:45:18 EST 2016
// Creator: jrsteven (on Linux gluon110.jlab.org 2.6.32-358.23.2.el6.x86_64 x86_64)
//

#include "JEventProcessor_TAGGER_online.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_TAGGER_online());
}
} // "C"


//------------------
// JEventProcessor_TAGGER_online (Constructor)
//------------------
JEventProcessor_TAGGER_online::JEventProcessor_TAGGER_online()
{

}

//------------------
// ~JEventProcessor_TAGGER_online (Destructor)
//------------------
JEventProcessor_TAGGER_online::~JEventProcessor_TAGGER_online()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_TAGGER_online::init(void)
{
        japp->RootWriteLock();
	gDirectory->Cd("/");
	new TDirectoryFile("TAGGER_ONLINE", "TAGGER_ONLINE");
	gDirectory->cd("TAGGER_ONLINE");

	dTaggerEnergy_DeltaTSC = new TH2D("TaggerEnergy_DeltaTSC", "Tagger Energy vs. #Delta t (TAG-SC); #Delta t (TAG-SC); Tagger Energy", 200, -100, 100, 240, 0., 12.); 

	gDirectory->cd("..");
        japp->RootUnLock();

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_TAGGER_online::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_TAGGER_online::evnt(JEventLoop *loop, uint64_t eventnumber)
{
        vector<const DBeamPhoton*> locBeamPhotons;
	loop->Get(locBeamPhotons);

	vector<const DSCHit*> locSCHits;
	loop->Get(locSCHits);

	for(size_t loc_i = 0; loc_i < locBeamPhotons.size(); loc_i++) {
	  const DTAGMHit* locTAGMHit;	  
	  locBeamPhotons[loc_i]->GetSingle(locTAGMHit);
	  if(locTAGMHit != NULL && locTAGMHit->integral < 500.) continue;

	  for(size_t loc_j = 0; loc_j < locSCHits.size(); loc_j++) {
	    Double_t locDeltaT = locBeamPhotons[loc_i]->time() - locSCHits[loc_j]->t;
	    japp->RootWriteLock();
	    dTaggerEnergy_DeltaTSC->Fill(locDeltaT, locBeamPhotons[loc_i]->momentum().Mag());
	    japp->RootUnLock();

	  }
	}

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_TAGGER_online::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_TAGGER_online::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

