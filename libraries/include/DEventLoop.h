// Author: David Lawrence  June 24, 2004
// $Id$
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

class DEventLoop;
class DEventSource;
class DEventProcessor;

#include "DEvent.h"
#include "derror.h"

#define MAX_EVENT_PROCESSORS 256

class DEventLoop:public DEvent{
	public:
		/// Constructor for a DEventLoop object. This should normally be
		/// called with the same arguments passed to main() on program start up.
		DEventLoop(int narg, char *argv[]);
		~DEventLoop();

		derror_t AddProcessor(DEventProcessor *processor);
		derror_t Init(void);
		void GotoEvent(int eventno){goto_event = eventno;}
		derror_t OneEvent(void);
		derror_t Fini(void);
		derror_t Run(void);
		derror_t Run(DEventProcessor *proc);
		void Quit(void){quit=1;};
		float GetRate(void);
		derror_t PrintRate(void);
		
	private:
		int Nprocessors;
		DEventProcessor* processors[MAX_EVENT_PROCESSORS];
		DEventSource *source;
		time_t last_print_rate_time;
		int quit;
		int goto_event;
};


#endif //_DEVENT_LOOP_H_
