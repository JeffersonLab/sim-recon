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
		jerror_t LoadNextRunData();
		//jerror_t OpenNextFile();
		jerror_t LoadESData();

	protected:
		jerror_t LoadNextVersionRunRange();
		
	private:
	
		JEventSource *event_source;    //  the source we are actually reading from
		string esdb_connection;        //  connection string for database
		DESDBProvider *esdb;           //  the database connection
		bool es_data_loaded;
		int VERBOSE;
		
		bool load_all_skims;
		vector<string> skim_list;

		string timestamp;
		string grade;
		int min_run, max_run;
		int min_file, max_file;   // optionally specify file ranges
		bool run_range_set;
		bool run_period_set;
		
		// Store lists of runs - EventStore information is keyed off runs
		// The master DB query returns a list of run ranges that satisfies
		// some given requirements
		//vector<EventStore::RunRange> run_ranges;              
		//vector<EventStore::RunRange>::const_iterator current_run_range;
		EventStore::DataVersionList data_versions;              
		EventStore::DataVersionList::iterator current_data_version_itr;
		bool first_data_version;
		// the list of actual run numbers to process in the current run range
		vector<int32_t> current_run_numbers;
		vector<int32_t>::const_iterator current_run_itr;
		int current_graphid;
		bool first_run_in_range;

		// data file information
		map<int,string> data_file_map;   // store map of fid -> file name
		map<int,string> data_type_map;   // store map of fid -> data type
		int current_fid;
		
		// index data
		vector<EventStore::DESEventIndexData> event_index;
		int event_index_pos;
		
		// some more metadata
		map< string, pair<int,int> > run_period_map;

};


#endif // _DEventSourceEventStore_

