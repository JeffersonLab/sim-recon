// $Id$
//
//    File: JEventProcessor_merge_rawevents.cc
//

#include <math.h>

#include "JEventProcessor_merge_rawevents.h"
using namespace jana;

// Routine used to create our JEventProcessor
#include "JANA/JApplication.h"
#include "GlueX.h"
#include <vector>
#include <deque>
#include <string>
#include <iostream>
#include <algorithm>
#include <stdio.h>
#include <stdlib.h>
#include "DAQ/JEventSource_EVIO.h"

extern "C"{
  void InitPlugin(JApplication *app){
    InitJANAPlugin(app);
    app->AddProcessor(new JEventProcessor_merge_rawevents());
  }
} // "C"


//------------------
// JEventProcessor_merge_rawevents (Constructor)
//------------------
JEventProcessor_merge_rawevents::JEventProcessor_merge_rawevents() : dEventWriterEVIO(NULL)
{

}

//------------------
// ~JEventProcessor_merge_rawevents (Destructor)
//------------------
JEventProcessor_merge_rawevents::~JEventProcessor_merge_rawevents()
{
}

//------------------
// init
//------------------
jerror_t JEventProcessor_merge_rawevents::init(void)
{
    return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_merge_rawevents::brun(JEventLoop *loop, int32_t runnumber)
{
    // create a new EVIO writer if it doesn't exist
    if(dEventWriterEVIO == NULL) {
        dEventWriterEVIO = new DEventWriterEVIO(loop);
        dEventWriterEVIO->Set_MergeFiles(true);
    }

    return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_merge_rawevents::evnt(JEventLoop *loop, uint64_t eventnumber)
{
    // get EVIO information associated with the event
    JEvent& the_event = loop->GetJEvent();
    void* the_event_ref = the_event.GetRef();
    uint32_t* output_buffer = JEventSource_EVIO::GetEVIOBufferFromRef(the_event_ref);
    uint32_t  output_buffer_size = JEventSource_EVIO::GetEVIOBufferSizeFromRef(the_event_ref);

    //cout << "Writing out event " << eventnumber << " buffer size = " << (output_buffer_size/4) << " words"  << endl;

    // write the buffer out
    // WARNING: this will work for non-entangled events, but hasn't been tested for entagled EVIO events
    dEventWriterEVIO->Write_EVIOBuffer( loop, output_buffer, output_buffer_size, "merged" );
    
    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_merge_rawevents::erun(void)
{
  // This is called whenever the run number changes, before it is
  // changed to give you a chance to clean up before processing
  // events from the next run number.
  return NOERROR;
}

//------------------
// Fin
//------------------
jerror_t JEventProcessor_merge_rawevents::fini(void)
{
    // Called before program exit after event processing is finished.
    delete dEventWriterEVIO;
    return NOERROR;
}


