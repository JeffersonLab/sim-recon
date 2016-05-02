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
#include <string>

#include "DESDBProvider.h"
#include "DESDBProviderMySQL.h"

using namespace jana;
using namespace std;


class DEventSourceEventStore : public JEventSource {
	public:
		DEventSourceEventStore(const char* source_name);
		virtual ~DEventSourceEventStore();
		const char* className(void){return "DEventSourceEventStore";}

		jerror_t GetEvent(JEvent &event);
		void FreeEvent(JEvent &event);
		jerror_t GetObjects(JEvent &event, JFactory_base *factory);
		
		
	protected:
		void PrintGrades();
		void PrintSkims(string timestamp, string grade);
		//void PrintActualDate();
				 
	private:
	
		JEventSource *event_source;    //  the source we are actually reading from
		string esdb_connection;        //  connection string for database
		DESDBProvider *esdb;           //  the database connection
		
		// We tag which skims JEvents belong to using JEvent::SetStatusBit()
		// We can get away with this now, since no one else is using fields above 16 yet
		// Probably the scheme needs to change or we need our own fields
		int BASE_SKIM_INDEX;           // the first status bit that we use for EventStore
		int MAX_SKIM_INDEX;            // we can store 64 bits in the JEvent, so 64 - MAX_SKIM_INDEX

		vector<string> skim_list;

		int min_run, max_run;
};


#endif // _DEventSourceEventStore_

