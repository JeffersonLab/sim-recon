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
#include <map>
#include <utility>

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
		// reporting functions
		void PrintGrades();
		void PrintRunPeriods();
		void PrintSkims(string timestamp, string grade);
		//void PrintActualDate();
				 
		// utility functions
		jerror_t MoveToNextEvent();
		jerror_t OpenNextFile();

		
	private:
	
		JEventSource *event_source;    //  the source we are actually reading from
		string esdb_connection;        //  connection string for database
		DESDBProvider *esdb;           //  the database connection
		bool es_data_loaded;
		
		// We tag which skims JEvents belong to using JEvent::SetStatusBit()
		// We can get away with this now, since no one else is using fields above 16 yet
		// Probably the scheme needs to change or we need our own fields
		int BASE_SKIM_INDEX;           // the first status bit that we use for EventStore
		int MAX_SKIM_INDEX;            // we can store 64 bits in the JEvent, so 64 - MAX_SKIM_INDEX

		bool load_all_skims;
		vector<string> skim_list;

		string timestamp;
		string grade;
		int min_run, max_run;
		bool run_range_set;
		bool run_period_set;
		
		map< string, pair<int,int> > run_period_map;

		// file information - deprecated
		vector<string> data_files;   // store list of file names
		vector<string>::iterator current_file_itr;

		// store list of runs - EventStore information is keyed off runs
		vector<int> run_numbers;
		vector<int>::iterator current_run_itr;
		
};


#endif // _DEventSourceEventStore_

