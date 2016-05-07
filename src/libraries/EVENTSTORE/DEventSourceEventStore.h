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
#include "DESDBProviderSQLite.h"
#include "DEventStoreDefs.h"

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
		
		bool load_all_skims;
		vector<string> skim_list;

		string timestamp;
		string grade;
		int min_run, max_run;
		bool run_range_set;
		bool run_period_set;
		

		// file information - deprecated
		vector<string> data_files;   // store list of file names
		vector<string>::iterator current_file_itr;

		// Store lists of runs - EventStore information is keyed off runs
		// The master DB query returns a list of run ranges that satisfies
		// some given requirements
		vector<EventStore::RunRange> run_ranges;              
		vector<EventStore::RunRange>::const_iterator current_run_range;
		// the list of actual run numbers to process in the current run range
		vector<int32_t> current_run_numbers;
		vector<int32_t>::const_iterator current_run_itr;
		
		
		// some more metadata
		map< string, pair<int,int> > run_period_map;

};


#endif // _DEventSourceEventStore_

