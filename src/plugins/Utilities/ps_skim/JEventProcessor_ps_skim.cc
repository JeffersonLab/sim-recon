// $Id$
//
//    File: JEventProcessor_ps_skim.cc
// Created: Mon May 18 09:52:08 EDT 2015
// Creator: nsparks (on Linux cua2.jlab.org 2.6.32-431.5.1.el6.x86_64 x86_64)
//

#include "JEventProcessor_ps_skim.h"
using namespace jana;

#include "evio_writer/DEventWriterEVIO.h"
#include <TRIGGER/DL1Trigger.h>

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
    void InitPlugin(JApplication *app){
        InitJANAPlugin(app);
        app->AddProcessor(new JEventProcessor_ps_skim());
    }
} // "C"


//------------------
// JEventProcessor_ps_skim (Constructor)
//------------------
JEventProcessor_ps_skim::JEventProcessor_ps_skim()
{

}

//------------------
// ~JEventProcessor_ps_skim (Destructor)
//------------------
JEventProcessor_ps_skim::~JEventProcessor_ps_skim()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_ps_skim::init(void)
{

    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_ps_skim::brun(JEventLoop *eventLoop, int32_t runnumber)
{
    // This is called whenever the run number changes
    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_ps_skim::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    const DEventWriterEVIO* locEventWriterEVIO = NULL;
    loop->GetSingle(locEventWriterEVIO);
    // write out BOR events
    if(loop->GetJEvent().GetStatusBit(kSTATUS_BOR_EVENT)) {
        locEventWriterEVIO->Write_EVIOEvent(loop, "ps");
        return NOERROR;
    }
    // write out EPICS events
    if(loop->GetJEvent().GetStatusBit(kSTATUS_EPICS_EVENT)) {
        locEventWriterEVIO->Write_EVIOEvent(loop, "ps");
        return NOERROR;
    }
    // get trigger types
    const DL1Trigger *trig_words = NULL;
    uint32_t trig_mask, fp_trig_mask;
    try {
        loop->GetSingle(trig_words);
    } catch(...) {};
    if (trig_words) {
        trig_mask = trig_words->trig_mask;
        fp_trig_mask = trig_words->fp_trig_mask;
    }
    else {
        trig_mask = 0;
        fp_trig_mask = 0;
    }
    int trig_bits = fp_trig_mask > 0 ? 10 + fp_trig_mask:trig_mask;
    // skim PS triggers
    if (trig_bits==8) {
        locEventWriterEVIO->Write_EVIOEvent(loop, "ps");
        return NOERROR;
    }
    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_ps_skim::erun(void)
{
    // This is called whenever the run number changes, before it is
    // changed to give you a chance to clean up before processing
    // events from the next run number.
    return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_ps_skim::fini(void)
{
    // Called before program exit after event processing is finished.
    return NOERROR;
}

