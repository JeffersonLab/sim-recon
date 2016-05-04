// $Id$
//
//    File: DEventStoreEvent.h
// Creator: sdobbs 
//

#ifndef _DEventStoreEvent_
#define _DEventStoreEvent_

#include <set>
#include <string>

#include <JANA/JEventSource.h>

using namespace std;
using namespace jana;

class DEventStoreEvent
{
	public:
		DEventStoreEvent() : dEventSource(NULL), dSourceRef(NULL){}
		~DEventStoreEvent(){}

		//SOURCES
		void Set_EventSource(JEventSource* locEventSource){dEventSource = locEventSource;}
		JEventSource* Get_EventSource(void) const{return dEventSource;}

		//SKIMS
		void Add_Skim(string locSkim){dSkims.insert(locSkim);}
		bool Is_Skim(string locSkim) const{return (dSkims.find(locSkim) != dSkims.end());}
		const set<string>& Get_Skims(void) const{return dSkims;}

		//REFS
		void Set_SourceRef(void* locSourceRef){dSourceRef = locSourceRef;}
		void* Get_SourceRef(void) const{return dSourceRef;}

	private:
		JEventSource* dEventSource;   // the event source that made the event (e.g. REST)
		void* dSourceRef;           // the Ref for the associated event (e.g. hddm record)
		set<string> dSkims;    // what skims the event belongs to
};

#endif  // _DEventStoreEvent_

