//
// JEventProcessor_trigger_skims.cc
//
// JANA event processor plugin to skim various trigger types to EVIO files
//
// Author: Sean Dobbs
// Adapted from Ahmed Foda, 13-May-2016 copied from Paul Mattione's 2trackskim

#include "JEventProcessor_trigger_skims.h"
#include "TRIGGER/DL1Trigger.h"
#include "BCAL/DBCALHit.h"

// for initializing plugins
extern "C" {
   void InitPlugin(JApplication *app)
	{
		InitJANAPlugin(app);
		app->AddProcessor(new JEventProcessor_trigger_skims(), true);
   }
} // "extern C"


// variables to control which triggers get read out
bool write_out_bcal_led = true;
bool write_out_fcal_led = true;
bool write_out_random = true;


//-------------------------------
// init
//-------------------------------
jerror_t JEventProcessor_trigger_skims::init(void)
{
    int bcal_led_writeout_toggle = 1;
    int fcal_led_writeout_toggle = 1;
    int random_writeout_toggle = 1;

    gPARMS->SetDefaultParameter("TRIGSKIM:WRITEBCALLED", bcal_led_writeout_toggle, "Write out BCAL LED events");
    gPARMS->SetDefaultParameter("TRIGSKIM:WRITEFCALLED", fcal_led_writeout_toggle, "Write out FCAL LED events");
    gPARMS->SetDefaultParameter("TRIGSKIM:WRITERANDOM", random_writeout_toggle, "Write out random pulser events");

    if(bcal_led_writeout_toggle == 0)
        write_out_bcal_led = false;
    if(fcal_led_writeout_toggle == 0)
        write_out_fcal_led = false;
    if(random_writeout_toggle == 0)
        write_out_random = false;

    return NOERROR;
}

//-------------------------------
// brun
//-------------------------------
jerror_t JEventProcessor_trigger_skims::brun(JEventLoop *locEventLoop, int32_t runnumber)
{
   return NOERROR;
}

//-------------------------------
// evnt
//-------------------------------
jerror_t JEventProcessor_trigger_skims::evnt(JEventLoop *locEventLoop, uint64_t eventnumber)
{
    // Get EVIO writer 
	const DEventWriterEVIO* locEventWriterEVIO = NULL;
	locEventLoop->GetSingle(locEventWriterEVIO);

    // Save BOR events
    if(locEventLoop->GetJEvent().GetStatusBit(kSTATUS_BOR_EVENT)) {
        locEventWriterEVIO->Write_EVIOEvent( locEventLoop, "BCAL-LED" );
        locEventWriterEVIO->Write_EVIOEvent( locEventLoop, "FCAL-LED" );
        locEventWriterEVIO->Write_EVIOEvent( locEventLoop, "random" );
        return NOERROR;
    }

	// Save EPICS events
	vector<const DEPICSvalue*> locEPICSValues;
	locEventLoop->Get(locEPICSValues);
	if(!locEPICSValues.empty()) {
        if (write_out_bcal_led)
            locEventWriterEVIO->Write_EVIOEvent(locEventLoop, "BCAL-LED");
        if (write_out_fcal_led)
            locEventWriterEVIO->Write_EVIOEvent(locEventLoop, "FCAL-LED");
        if (write_out_random)
            locEventWriterEVIO->Write_EVIOEvent(locEventLoop, "random");
		return NOERROR;
	}

    bool cosmic_trigger = false;
	bool BCAL_LED_US_trigger = false;
    bool BCAL_LED_DS_trigger = false;
    bool FCAL_LED_trigger = false;
    bool random_trigger = false;

	const DL1Trigger *trig = NULL;
	try {
		locEventLoop->GetSingle(trig);
	} catch (...) {}

    // parse the triggers we want to save
	if (trig) {
		//printf("%5i  %5i | %5i  %5i  %5i | %i\n",
		//	   trig->trig_mask,trig->trig_mask & 0x1,
		//	   trig->fp_trig_mask, trig->fp_trig_mask & 0x100,trig->fp_trig_mask & 0x200,
		//	   trig->trig_mask && trig->fp_trig_mask);

		if (trig->trig_mask & 0x1) {
			// Cosmic trigger fired
			cosmic_trigger = true;
		}
        
        // Select triggers based on front panel inputs
        // Trigger bits start counting from 0
		if (trig->fp_trig_mask & 0x100) {   // Trigger front-panel bit 8
			// Upstream BCAL LED trigger fired
			BCAL_LED_US_trigger = true;
		}
		if (trig->fp_trig_mask & 0x200) {   // Trigger front-panel bit 9
			// Downstream BCAL LED trigger fired
			BCAL_LED_DS_trigger = true;
		}
		if (trig->fp_trig_mask & 0x800) {  // Trigger front-panel bit 11
			// Periodic pulser trigger fired
			random_trigger = true;
		}
		if (trig->fp_trig_mask & 0x004) {   // Trigger front-panel bit 2
			// FCAL LED trigger fired
			FCAL_LED_trigger = true;
		}
	} 

    // Do some backup calculations for runs in which the BCAL LED trigger did not latch correctly
    vector<const DBCALHit *> bcal_hits;
    double total_bcal_energy = 0.;
    if(write_out_bcal_led) {
        for(unsigned int i=0; i<bcal_hits.size(); i++) {
            total_bcal_energy += bcal_hits[i]->E;
        }
    }

    // Save events to skim file

    // Save BCAL trigger if:
    // 1. Trigger front-panel bits 8 or 9
    // 2. Total energy in BCAL > 12 GeV
    // 3. Number of hits in BCAL > 200
    bool save_BCAL_LED_event = BCAL_LED_US_trigger || BCAL_LED_DS_trigger
        || (bcal_hits.size() >= 200) || (total_bcal_energy > 12.);
	if (write_out_bcal_led && save_BCAL_LED_event) {
        locEventWriterEVIO->Write_EVIOEvent(locEventLoop, "BCAL-LED");
    }
	if (write_out_fcal_led && FCAL_LED_trigger) {
        locEventWriterEVIO->Write_EVIOEvent(locEventLoop, "FCAL-LED");
    }
	if (write_out_random && random_trigger) {
        locEventWriterEVIO->Write_EVIOEvent(locEventLoop, "random");
    }

   return NOERROR;
}

//-------------------------------
// erun
//-------------------------------
jerror_t JEventProcessor_trigger_skims::erun(void)
{
   return NOERROR;
}

//-------------------------------
// fini
//-------------------------------
jerror_t JEventProcessor_trigger_skims::fini(void)
{
   return NOERROR;
}

