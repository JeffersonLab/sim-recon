// $Id$
//
//    File: DEventSourceEventStore.cc
// Creator: sdobbs 
//

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <limits>

using namespace std;

#include <DANA/DApplication.h>
#include <DANA/DStatusBits.h>

#include <TRandom3.h>

#include <HDDM/DEventSourceREST.h>
#include <HDDM/DEventSourceHDDM.h>

#include "DEventSourceEventStore.h"
#include "DEventStoreEvent.h"
#include "DESSkimData.h"

static string EventstoreQueryHelp() {
	string str = "For more information, go to https://github.com/JeffersonLab/HDEventStore/wiki";
	return str;
}

// File-scope variables
static bool ES_TEST_MODE = false;

//---------------------------------
// AcceptRunRange
//---------------------------------
static bool AcceptRunRange(EventStore::RunRange &the_run_range, int min_run, int max_run)
{
	// reject the run range if it's not in the range we're interested in
	// otherwise, modify it if needed to make sure the range is within what we're interested in
	if(the_run_range.second < min_run || the_run_range.first > max_run) {
		// if the run range lies outside of what we want, reject it
		return false;
	} else if(the_run_range.first > min_run && the_run_range.second < max_run) {
		// if the run range is completely contained in our range, keep it
		return true;
	} else {
		// otherwise trim as needed
		if(the_run_range.first < min_run) {
			the_run_range.first = min_run;
		}
		if(the_run_range.second > max_run) {
			the_run_range.second = max_run;
		}
	} 
}


//---------------------------------
// DEventSourceEventStore    (Constructor)
//---------------------------------
DEventSourceEventStore::DEventSourceEventStore(const char* source_name):JEventSource(source_name)
{
	// initialize data members
	VERBOSE = 0;
	es_data_loaded = false;
	event_source = NULL;
	min_run = 0;
	max_run = numeric_limits<int>::max();   // default to something ridiculously large
	esdb_connection = "mysql://es_user@hallddb.jlab.org/EventStoreTMP";    // default to main JLab ES server
	
	// eventstore parameters
	// for more information, see ...
	//grade = "physics";        // REST physics data
	//grade = "recon";        // REST physics data
	grade = "recon-unchecked";        // DEBUG!!
	timestamp = "21000501";   // use the latest data - this should never happen, anyway...
	load_all_skims = true;    // default to processing all skims
	run_period_set = false;
	run_range_set = false;
	current_graphid = -1;
	current_fid = -1;
	
	// load run period mapping
	// hardcode for now, this info should move to RCDB...
	run_period_map["2015-03"] = pair<int,int>(2607,3385);
	run_period_map["2016-02"] = pair<int,int>(10000,20000);
		
	// read in configurations
	// priority:  JANA command line -> environment variable -> default
	if(getenv("EVENTSTORE_CONNECTION") != NULL)
		esdb_connection = getenv("EVENTSTORE_CONNECTION");
	gPARMS->SetDefaultParameter("ESDB:DB_CONNECTION", esdb_connection,
								"Specification of EventStore DB connection.");
	gPARMS->SetDefaultParameter("ESDB:VERBOSE", VERBOSE,
								"Level of verbose output.");

	int test_mode_flag = 0;
	gPARMS->SetDefaultParameter("ESDB:TEST_MODE", test_mode_flag,
								"Toggle test mode features");
	if(test_mode_flag != 0) {
		ES_TEST_MODE = true;
		if(gRandom == NULL)
			gRandom = new TRandom3(0);
	}
	
	// First, parse the eventsource query
	// For details of the query format, see: <...>
	
	// Tokenize the query string
	string es_query(source_name);
	istringstream iss(es_query);
	vector<string> tokens;
	copy(istream_iterator<string>(iss),
    	 istream_iterator<string>(),
     	 back_inserter(tokens));
	
	////////////////////////////////////////////////////////////

	// initialize database connection
	if(esdb_connection.substr(0,8) == "mysql://") {
		cout << "Connecting to MySQL..." << endl;
		esdb = static_cast<DESDBProvider *>(new DESDBProviderMySQL(esdb_connection));
	} else if(esdb_connection.substr(0,8) == "sqlite://") {
		cout << "Connecting to SQLite..." << endl;	
		esdb = static_cast<DESDBProvider *>(new DESDBProviderSQLite(esdb_connection));
	} 
	
	// Connect to database
	esdb->Open();
	
	////////////////////////////////////////////////////////////
	// parse the ES command
	//   format: eventstore in YYYYMMDD <grade> [more optional args]
	// Check query header
	if( (tokens[0] != "eventstore") || tokens.size() < 3)
		throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
		
	if( (tokens[1] == "in") ) {   // read data from the data base
		timestamp = tokens[2];    // require a timestamp
		
		if(tokens.size() > 3) {
			// parse the rest
			size_t token_ind = 3;
			
			// the next argument is the data grade to use
			grade = tokens[token_ind++];
			
			while(token_ind < tokens.size()) {
				if(tokens[token_ind] == "runs") {
					// format: runs [MAX RUN] [MIN RUN]
					if(run_period_set) 
						throw JException("Cannot set run range and run period in the same command!");
					
					// make sure there's enough args
					if(tokens.size() - token_ind < 3)
						throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
						
					min_run = atoi(tokens[token_ind+1].c_str());   // ERROR CHECK!!
					max_run = atoi(tokens[token_ind+2].c_str());   // ERROR CHECK!!
	
					// sanity check
					if(max_run < min_run) {
						throw JException("Maximum run must be larger than minimum run!");
					}
					token_ind += 3;
				} else if(tokens[token_ind] == "files") {
					// format: files [MAX FILES] [MIN FILES]
					// make sure there's enough args
					if(tokens.size() - token_ind < 3)
						throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
						
					min_file = atoi(tokens[token_ind+1].c_str());   // ERROR CHECK!!
					max_file = atoi(tokens[token_ind+2].c_str());   // ERROR CHECK!!
	
					// sanity check
					if(max_run < min_run) {
						throw JException("Maximum run must be larger than minimum run!");
					}
					token_ind += 3;
				} else if(tokens[token_ind] == "run_File") {
					// run_File [RUN] [FILE]
					// make sure there's enough args
					if(tokens.size() - token_ind < 3)
						throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
						
					min_run = atoi(tokens[token_ind+1].c_str());   // ERROR CHECK!!
					min_file = atoi(tokens[token_ind+2].c_str());   // ERROR CHECK!!
					max_run = min_run;
					max_file = min_file;
	
					// sanity check
					if(max_run < min_run) {
						throw JException("Maximum run must be larger than minimum run!");
					}
					token_ind += 3;
				} else if(tokens[token_ind] == "run_period") {
					// format: run_periods [RUN PERIOD STRING]
					if(run_range_set) 
						throw JException("Cannot set run range and run period in the same command!");

					// make sure there's enough args
					if(tokens.size() - token_ind < 2)
						throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
						
					map< string, pair<int,int> >::iterator run_period_itr = run_period_map.find(tokens[token_ind+1]);
					if(run_period_itr == run_period_map.end()) {
						// a bad run period was specified...
						PrintRunPeriods();
						throw JException("Invalid ES query = " + es_query + "\n");
					}
					min_run = run_period_itr->second.first;
					max_run = run_period_itr->second.second;
					token_ind += 2;
				} else if(tokens[token_ind] == "skims") {
					load_all_skims = false;
					// for the skims command, assume the rest of the arguments are skim names
					while(token_ind++ < tokens.size()) {
						skim_list.push_back(tokens[token_ind]);
					}
					
				} else {
					// require a valid command
					throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
				}
			}
			
			// sanity check - make sure the grade exists!
			vector<string> grades_in_db = esdb->GetGrades();
			
			vector<string>::iterator grade_itr = find(grades_in_db.begin(), grades_in_db.end(), grade);
			if(grade_itr == grades_in_db.end()) {
				jerr << "Could not find grade \'" << grade << "\' in DB!" << endl;
				PrintGrades();
				throw JException("Invalid ES query = " + es_query + "\n");
			}
		}
		
		
		// debugging stuff
		if(ES_TEST_MODE) {
			if(skim_list.size() == 0) {
				skim_list.push_back("pi0");
				skim_list.push_back("eta");
				skim_list.push_back("rho");
				skim_list.push_back("omega");				
			}
		}
	} else if( (tokens[1] == "info") ) {   // query information from the DB
		// runs <datestamp> [<grade> [<skim>] ] [runs <min> [<max>]] 
    	//    [<run info query>]    prints available runs in DB  that match criteria
        //PrintRuns();
                           
  		// grades                  prints available grades in DB
		if(tokens[2] == "grades") {
			PrintGrades();
                           
  		// run_periods                  prints available run periods
		} else if(tokens[2] == "run_periods") {
			PrintRunPeriods();

  		// skims <datestamp> <grade>  print skims available for the grade
        } else if(tokens[2] == "skims") {
        	if(tokens.size() < 5)
        		throw JException("Invalid query:  skims <datestamp> <grade>" );

        	PrintSkims(tokens[3], tokens[4]);    
		}
  		// actualDate <datestamp> <grade>    print actual internal date used
        //} else if(tokens[2] == "actualDate") {
        //	if(tokens.size() < 5)
        //		throw JException("Invalid query = actualDate <datestamp> <grade>" );
		//
        //	PrintActualDate();     
		//}
  		// versions   <datestamp> <grade> [-verbose]
       	//                   print run ranges versus versions
        //                use -verbose option for more details
	
		// We're just querying the database, so we can quit here
		exit(0);
	} else {
		throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
	}

	////////////////////////////////////////////////////////////
	LoadESData();
	es_data_loaded = true;
		
	/*
	if(ES_TEST_MODE)   // if we're testing, don't make any more checks 
		return;
		
	// load some data here
		
	// make sure we've found any files
	*/
}

//---------------------------------
// ~DEventSourceEventStore    (Destructor)
//---------------------------------
DEventSourceEventStore::~DEventSourceEventStore()
{
	if (event_source != NULL)
		delete event_source;
}

//---------------------------------
// LoadESData
// 
// Load EventStore data from database
// This should generally only be called once during program execution
//---------------------------------
jerror_t DEventSourceEventStore::LoadESData()
{
	if(es_data_loaded) {
		throw JException("Tried to load EventStore data twice! This should not happen! Exiting...");
	}

	if(VERBOSE>0) {
		jout << "EventStore: Loading data versions..." << endl;
	}

	if(load_all_skims) {
		// if the user didn't ask for specific skims, then load all of them
		skim_list = esdb->GetSkims(timestamp, grade);
	}
	
	// Now, load the data versions and corresponding run ranges
	//EventStore::RunVersionList in_run_ranges = esdb->RunVersionList(timestamp, grade);
	data_versions = esdb->GetDataVersions(timestamp, grade);
	if(data_versions.size() == 0) {
		jerr << "No runs exist for grade = " << grade 
			 << " and timestamp = " << timestamp << endl;
		throw JException("Problems loading EventStore data!");
	}
	current_data_version_itr = data_versions.begin();
	first_data_version = true;

	if(VERBOSE>1) {
		jout << "EventStore: Loaded versions:" << endl;
		for (auto the_data_version : data_versions) {
			jout << "    " << the_data_version.second << "  runs = ( "
				 << the_data_version.first.first << ", " << the_data_version.first.second << ")" << endl;
		}
	}

	// Set everything up, and we should be ready to go
	if(LoadNextVersionRunRange() != NOERROR) {
		throw JException("Problems loading EventStore data!");
	}
}

//---------------------------------
// LoadNextVersionRunRange
//---------------------------------
jerror_t DEventSourceEventStore::LoadNextVersionRunRange()
{
	if(VERBOSE>0) {
		jout << "EventStore: Loading next run range..." << endl;
	}

	// generally, we want to move to the next data version, but we need
	// to make we analyze the first
	if(!first_data_version)
		current_data_version_itr++;
	else
		first_data_version = false;

	// unpack the data version info
	EventStore::RunRange the_run_range = current_data_version_itr->first;
	int the_graphid = current_data_version_itr->second;

	// make sure this run range is within the range of runs we are interested
	while(!AcceptRunRange(the_run_range, min_run, max_run)) {
		current_data_version_itr++;
		if(current_data_version_itr == data_versions.end())
			break;
	}
	
	if(current_data_version_itr == data_versions.end())
		return NO_MORE_EVENTS_IN_SOURCE;
		
	// load run number list for this range
	string all_skim("all");
	current_run_numbers = esdb->GetRunList(the_run_range, the_graphid, all_skim);
	
	if(VERBOSE>1) {
		jout << "EventStore: Loaded runs:  ";
		for (auto run_number : current_run_numbers) {
			cout << " " << run_number;
		}
		cout << endl;
	}

	if(current_run_numbers.size() ==  0) // sanity check
		return LoadNextVersionRunRange();
	
	current_graphid = the_graphid;
	current_run_itr = current_run_numbers.begin();
	first_run_in_range = true;
	return LoadNextRunData();
}

//---------------------------------
// LoadNextRunData
//---------------------------------
jerror_t DEventSourceEventStore::LoadNextRunData()
{
	// generally, we want to move to the next run in the range, but we need
	// to make we analyze the first
	if(!first_run_in_range)
		current_run_itr++;
	else
		first_run_in_range = false;

	if(VERBOSE>0) {
		jout << "EventStore: Loading data for run " << (*current_run_itr) << "..." << endl;
	}

	// first load the primary index for the run
	string all_skim("all");
	string index_file_name = esdb->GetKeyFileName(current_graphid, all_skim, *current_run_itr);

	// LOAD INDEX

	if(VERBOSE>1) {
		jout << "EventStore: Loaded master index" << endl;
	}
	
	// then load the indices for the skims
	for(auto skim_name : skim_list) {
		string skim_index_file_name = esdb->GetKeyFileName(current_graphid, skim_name, *current_run_itr);

		// LOAD INDEX
		
		if(VERBOSE>1) {
			jout << "EventStore: Loaded skim " << skim_name << endl;
		}

	}
	
	// then load information on the data files for this run
	vector< pair<string,string> > data_file_type_vec = esdb->GetDataFileNameTypePairs(current_graphid, all_skim, *current_run_itr);
	// index information by database file id
	data_file_map.clear();
	for(auto data_file_type : data_file_type_vec) {
		int the_fid = esdb->GetFID(data_file_type.first);
		data_file_map[the_fid] = data_file_type.first;
		data_type_map[the_fid] = data_file_type.second;
	}

	if(VERBOSE>1) {
		jout << "EventStore: Data files = ";
		for(auto data_file_type : data_file_type_vec) 
			cout << data_file_type.first << " ";
		cout << endl;
	}
	
	// start from the beginning of the run
	event_index_pos = 0;
}

//---------------------------------
// GetEvent
//---------------------------------
jerror_t DEventSourceEventStore::GetEvent(JEvent &event)
{
	// make sure that all of the EventStore metadata has been loaded
	if(!es_data_loaded) {
		LoadESData();
		es_data_loaded = true;   
	}
	
	// FOR DEBUGGING - EMIT EVENTS FOREVER
	if(ES_TEST_MODE) {
		// output some fake event with skim information
    	event.SetEventNumber(1);
    	event.SetRunNumber(10000);
    	event.SetJEventSource(static_cast<JEventSource *>(this));
   		//event.SetRef(NULL);
    	event.SetStatusBit(kSTATUS_FROM_FILE);
    	event.SetStatusBit(kSTATUS_PHYSICS_EVENT);

		DEventStoreEvent *the_es_event = new DEventStoreEvent();
		event.SetRef(the_es_event);
		for(int i=0; i<4; i++)
			if(gRandom->Uniform() < 0.5)
				the_es_event->Add_Skim(skim_list[i]);

		return NOERROR;
	}
		
	// skip to next event
	jerror_t retval;
	if( (retval = MoveToNextEvent()) != NOERROR)
		return retval;   // if we can't get to another event, then we're done
		
	// make sure the file is open
	if(event_source == NULL)
		return NO_MORE_EVENTS_IN_SOURCE;

	// read the next event in
	retval = event_source->GetEvent(event);
	if(retval == NOERROR) {
		// To store the skim and other EventStore information for the event
		// we wrap the actual event data and store our information in the wrapper
		DEventStoreEvent *the_es_event = new DEventStoreEvent();
		the_es_event->Set_EventSource(event_source);
		the_es_event->Set_SourceRef(event.GetRef());    // save the actual event data
		event.SetRef(the_es_event);
	    event.SetStatusBit(kSTATUS_FROM_FILE);
	    event.SetStatusBit(kSTATUS_PHYSICS_EVENT);

		// tag event with skims
		;
	} else if(retval == NO_MORE_EVENTS_IN_SOURCE) {   
		// if the source is empty, close the current one, then move to the next
		delete event_source;
		event_source = NULL;
	} else {   // if there'a another error, then pass it on...
		return retval;
	}
	
	return NOERROR;
}

//---------------------------------
// FreeEvent
//---------------------------------
void DEventSourceEventStore::FreeEvent(JEvent &event)
{
	if(event_source != NULL)
		event_source->FreeEvent(event);
}

//---------------------------------
// GetObjects
//---------------------------------
jerror_t DEventSourceEventStore::GetObjects(JEvent &event, JFactory_base *factory)
{
	/// This gets called through the virtual method of the
	/// JEventSource base class. It creates the objects of the type
	/// on which factory is based.

	// We must have a factory to hold the data
	if(!factory) throw RESOURCE_UNAVAILABLE;

	// return meta-EventStore information
    string dataClassName = factory->GetDataClassName();
	
	if (dataClassName =="DESSkimData") {
		JFactory<DESSkimData> *essd_factory = dynamic_cast<JFactory<DESSkimData>*>(factory);
		
		DEventStoreEvent *the_es_event = static_cast<DEventStoreEvent *>(event.GetRef());
		DESSkimData *skim_data = new DESSkimData(the_es_event->Get_Skims(), skim_list);

		vector<DESSkimData*> skim_data_vec(1, skim_data);
		essd_factory->CopyTo(skim_data_vec);

		return NOERROR;
	}
	
	if(!event_source) throw RESOURCE_UNAVAILABLE;
	
	// See GetEvent() for the motivation for this
	// Unwrap the event...
	DEventStoreEvent *the_es_event = static_cast<DEventStoreEvent *>(event.GetRef());
	event.SetRef(the_es_event->Get_SourceRef());
	// ..now grab the objects...
	jerror_t retval = event_source->GetObjects(event, factory);
	// ...and wrap it back up
	event.SetRef(the_es_event);

	return retval;
}

//---------------------------------
// MoveToNextEvent
//---------------------------------
jerror_t DEventSourceEventStore::MoveToNextEvent()
{
	//
	
	
	// move to next event
	event_index_pos++;

	// if we're loading all of the skims, then we don't need to skip any events
	if(load_all_skims)
		return NOERROR;
		
	// see if we need to change files
	if(current_fid != event_index[event_index_pos].fid) {
		delete event_source;
		current_fid = event_index[event_index_pos].fid;
		// create new file
		if(data_type_map[current_fid] == "rest") {
			event_source = static_cast<JEventSource*>(new DEventSourceREST(data_file_map[current_fid].c_str()));
		} else if(data_type_map[current_fid] == "sim") {
			event_source = static_cast<JEventSource*>(new DEventSourceHDDM(data_file_map[current_fid].c_str()));
		} 
		//else if(data_type_map[fid] == "evio") {   // FUTURE
		//event_source = JEventSourceEVIO(data_file_map[fid]);
		//}
	}
		
	return NOERROR;
}

//---------------------------------
// OpenNextFile
//---------------------------------
/*
jerror_t DEventSourceEventStore::OpenNextFile()
{

	
	// if there's a current file open, close it so that we don't leak memory
	if(event_source != NULL) {
		delete event_source;
		event_source = NULL;
	}

	//Get generators
	vector<JEventSourceGenerator*> locEventSourceGenerators = japp->GetEventSourceGenerators();
	
	//Get event source
	while( (event_source == NULL) && (current_file_itr != data_files.end()) ) {

		// Loop over JEventSourceGenerator objects and find the one
		// (if any) that has the highest chance of being able to read
		// this source. The return value of 
		// JEventSourceGenerator::CheckOpenable(source) is a liklihood that
		// the named source can be read by the JEventSource objects
		// created by the generator. In most cases, the liklihood will
		// be either 0.0 or 1.0. In the case that 2 or more generators return
		// equal liklihoods, the first one in the list will be used.
		JEventSourceGenerator* locEventSourceGenerator = NULL;
		double liklihood = 0.0;
		string locFileName = current_file_itr->c_str();
		for(unsigned int i=0; i<locEventSourceGenerators.size(); i++)
		{
			double my_liklihood = locEventSourceGenerators[i]->CheckOpenable(locFileName);
			if(my_liklihood > liklihood)
			{
				liklihood = my_liklihood;
				locEventSourceGenerator = locEventSourceGenerators[i];
			}
		}

		if(locEventSourceGenerator != NULL)
		{
			jout<<"Opening source \""<<locFileName<<"\" of type: "<<locEventSourceGenerator->Description()<<endl;
			event_source = locEventSourceGenerator->MakeJEventSource(locFileName);
		}

		if(event_source == NULL){
			jerr<<endl;
			jerr<<"  xxxxxxxxxxxx  Unable to open event source \""<<locFileName<<"\"!  xxxxxxxxxxxx"<<endl;
		}

		current_file_itr++;
	}
	
	// error check
	if(event_source == NULL)
		return NO_MORE_EVENT_SOURCES;
	else
		return NOERROR;
}
*/


//---------------------------------
// PrintGrades
//---------------------------------
void DEventSourceEventStore::PrintGrades()
{
	vector<string> grades = esdb->GetGrades();
	
	// print out information
	cout << endl << "Available grades:" << endl;
	for(vector<string>::iterator it = grades.begin();
		it != grades.end(); it++)
			cout << "  " << *it << endl;
	cout << endl;
}

//---------------------------------
// PrintRunPeriods
//---------------------------------
void DEventSourceEventStore::PrintRunPeriods()
{
	vector<string> grades = esdb->GetGrades();
	
	// print out information
	cout << endl << "Available Run Periods:" << endl;
	for(map< string, pair<int,int> >::iterator it = run_period_map.begin();
		it != run_period_map.end(); it++) {
			pair<int,int> &the_run_range = it->second;
			cout << "  " << it->first << ":  " 
				 << the_run_range.first << " - " << the_run_range.second << endl;
	}
	cout << endl;
}

//---------------------------------
// PrintSkims
//---------------------------------
void DEventSourceEventStore::PrintSkims(string datestamp, string grade)
{
	vector<string> skims = esdb->GetSkims(datestamp, grade);
	
	// print out information
	cout << endl << "Available skims for grade " << grade << ":" << endl;
	for(vector<string>::iterator it = skims.begin();
		it != skims.end(); it++)
			cout << "  " << *it << endl;
	cout << endl;
}

