// Author: David Lawrence  June 24, 2004
//
//
// DEventLoop
//
/// Hall-D DANA package:
/// The DEventLoop class implements the main event processing loop
/// for a Hall-D analysis package. Mainly, it instatiates the
/// proper DEventSource variant (DEventSourceET, DEventSourceHDDM,...)
/// and uses it to repeatedly read in events. It also keeps a list
/// of DEventProcessor objects to call for each event as it is
/// read in. The idea is that the user will define their own
/// class based on DEventProcessor. A pointer to this object will
/// be passed to DEventLoop to add to its standard list of objects
/// to call as it loops over events.
///

#ifndef _DEVENT_LOOP_H_
#define _DEVENT_LOOP_H_

#include "DEventSource.h"
#include "DEventProcessor.h"
#include "derror.h"

#define MAX_EVENT_PROCESSORS 256

class DEventLoop
{
	enum EVENT_PROCESSORS{
		EVENT_PROCESSOR_NONE			= 0x0000,
		EVENT_PROCESSOR_TAGGER		= 0x0001,
		EVENT_PROCESSOR_UPV			= 0x0002,
		EVENT_PROCESSOR_VERTEX		= 0x0004,
		EVENT_PROCESSOR_BCAL			= 0x0008,
		EVENT_PROCESSOR_CDC			= 0x0010,
		EVENT_PROCESSOR_FDC			= 0x0020,
		EVENT_PROCESSOR_CHERENKOV	= 0x0040,
		EVENT_PROCESSOR_TOF			= 0x0080,
		EVENT_PROCESSOR_FCAL			= 0x0100,
		EVENT_PROCESSOR_TRIGGER		= 0x0200,
		EVENT_PROCESSOR_TRACKING	= 0x0400,
		EVENT_PROCESSOR_PID			= 0x0800,
		EVENT_PROCESSOR_ALL			= 0xFFFF
	};

	public:
		/// Constructor for a DEventLoop object. This should normally be
		/// called with the same arguments passed to main() on program start up.
		DEventLoop(int narg, char *argv[]);
		~DEventLoop();

		derror_t AddProcessor(DEventProcessor *processor);
		derror_t SetStandardProcessors(EVENT_PROCESSORS processmask);
		derror_t Run(void);
		derror_t Run(DEventProcessor *proc);
		derror_t Run(EVENT_PROCESSORS mask);
		derror_t Run(DEventProcessor *proc, EVENT_PROCESSORS mask);
		
		EVENT_PROCESSORS proc_mask;
	
	private:
		int Nprocessors;
		DEventProcessor* processors[MAX_EVENT_PROCESSORS];
		DEventSource *source;
};


#endif //_DEVENT_LOOP_H_
