// $Id$
//
//    File: DEventSourceEventStore.h
// Creator: sdobbs
//

#ifndef _DEventSourceEventStore_
#define _DEventSourceEventStore_

#include <JANA/jerror.h>
#include <JANA/JEventSource.h>
#include <JANA/JEvent.h>

#include <vector>

using namespace jana;
using namespace std;


class DEventSourceEventStore:public JEventSource{
	public:
		DEventSourceEventStore(const char* source_name);
		virtual ~DEventSourceEventStore();

		jerror_t GetEvent(JEvent &event);
		void FreeEvent(JEvent &event);
		jerror_t GetObjects(JEvent &event, JFactory_base *factory);
		
		
	protected:
		
	
	private:
	
		JEventSource *event_source;    //  the source we are actually reading from
		string esdb_connection;        //  connection string for database

		int min_run, max_run;
};


#endif // _DEventSourceEventStore_

