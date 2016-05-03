// $Id$
//
//    File: DEventStoreEvent.h
// Creator: sdobbs 
//


#ifndef _DEventStoreEvent_
#define _DEventStoreEvent_

#include <JANA/JEventSource.h>
#include <JANA/JObject.h>
#include <JANA/JEvent.h>
using namespace jana;

#include <vector>
#include <set>
#include <string>
#include <iostream>
using namespace std;


class DEventStoreEvent {
  public:
    DEventStoreEvent() { event_source = NULL; data = NULL; }
    ~DEventStoreEvent() {}
  
	// data public for now, wrap it later...
	JEventSource *event_source;   // the event source that made the event
	set<string> skims;    // what skims the event belongs to
	void *data;           // the Ref for the associated event

};


#endif  // _DEventStoreEvent_