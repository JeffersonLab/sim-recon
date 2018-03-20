//
// JEventProcessor_randomtrigger_skim.cc
//
// JANA event processor plugin to skim random triggers to HDDM files
//
// Author: Sean Dobbs

#include "JEventProcessor_randomtrigger_skim.h"
#include "TRIGGER/DL1Trigger.h"
#include "BCAL/DBCALHit.h"
#include "DAQ/DL1Info.h"
#include <HDDM/DEventWriterHDDM.h>
#include <DAQ/DBeamCurrent.h>

// for initializing plugins
extern "C" {
   void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_randomtrigger_skim(), true);
   }
} // "extern C"


// variables to control which triggers get read out

//-------------------------------
// init
//-------------------------------
jerror_t JEventProcessor_randomtrigger_skim::init(void)
{
    return NOERROR;
}

//-------------------------------
// brun
//-------------------------------
jerror_t JEventProcessor_randomtrigger_skim::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
    dBeamCurrentFactory = new DBeamCurrent_factory();
    dBeamCurrentFactory->init();
    dBeamCurrentFactory->brun(locEventLoop, runnumber);
   
    return NOERROR;
}

//-------------------------------
// evnt
//-------------------------------
jerror_t JEventProcessor_randomtrigger_skim::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
    // Get HDDM writer
    vector<const DEventWriterHDDM*> locEventWriterHDDMVector;
    locEventLoop->Get(locEventWriterHDDMVector);

    // beam current and fiducial definition
    vector<const DBeamCurrent*> beamCurrent;
    locEventLoop->Get(beamCurrent);

    //bool is_cosmic_trigger = false;
    bool is_random_trigger = false;

	const DL1Trigger *trig = NULL;
	try {
		locEventLoop->GetSingle(trig);
	} catch (...) {}

    // parse the triggers we want to save
	if (trig) {

		if (trig->fp_trig_mask & 0x800) {  // Trigger front-panel bit 11
			// Periodic pulser trigger fired
			is_random_trigger = true;
		}

	} 

    // make sure this is a random trigger event
    if(!is_random_trigger)
        return NOERROR;
    
    // make sure we can perform a fiducial beam current cut
    if(beamCurrent.empty())
        return NOERROR;

    // make sure the beam is on
    if(!beamCurrent[0]->is_fiducial)
        return NOERROR;
        

    // Save events to skim file
    locEventWriterHDDMVector[0]->Write_HDDMEvent(locEventLoop, "random"); 

    return NOERROR;
}

//-------------------------------
// erun
//-------------------------------
jerror_t JEventProcessor_randomtrigger_skim::erun(void)
{
   return NOERROR;
}

//-------------------------------
// fini
//-------------------------------
jerror_t JEventProcessor_randomtrigger_skim::fini(void)
{
   return NOERROR;
}

