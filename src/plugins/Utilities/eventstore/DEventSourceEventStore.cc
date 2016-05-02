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
using namespace std;

#include <DANA/DApplication.h>
#include <DANA/DStatusBits.h>

#include <TRandom3.h>

#include "DEventSourceEventStore.h"
#include "DESSkimData.h"

static string EventstoreQueryHelp() {
	string str = "this is a help string";
	return str;
}


// Various variables
static bool TEST_MODE = false;


//---------------------------------
// DEventSourceEventStore    (Constructor)
//---------------------------------
DEventSourceEventStore::DEventSourceEventStore(const char* source_name):JEventSource(source_name)
{
	// initialize data members
	event_source = NULL;
	min_run = 0;
	max_run = INT_MAX;   // default to something ridiculously large
	esdb_connection = "mysql://es_user@hallddb.jlab.org/EventStoreTMP";    // default to main JLab ES server
	BASE_SKIM_INDEX = 20;
	MAX_SKIM_INDEX = 64 - BASE_SKIM_INDEX;
		
	// read in configurations
	// priority:  JANA command line -> environment variable -> default
	if(getenv("EVENTSTORE_CONNECTION") != NULL)
		esdb_connection = getenv("EVENTSTORE_CONNECTION");
	gPARMS->SetDefaultParameter("ESDB:CONNECTION", esdb_connection,
								"Specification of EventStore DB connection.");
	
	int test_mode_flag = 0;
	gPARMS->SetDefaultParameter("ESDB:TEST_MODE", test_mode_flag,
								"Toggle test mode features");
	if(test_mode_flag != 0) {
		TEST_MODE = true;
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
		cout << "MySQL connection" << endl;
		esdb = static_cast<DESDBProvider *>(new DESDBProviderMySQL(esdb_connection));
	} else if(esdb_connection.substr(0,8) == "sqlite://") {
		cout << "SQLite connection" << endl;	
	} 
	
	// Connect to database
	esdb->Open();
	
	////////////////////////////////////////////////////////////
	// parse the ES command
	// Check query header
	if( (tokens[0] != "eventstore") || tokens.size() < 3)
		throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
		
	if( (tokens[1] == "in") ) {   // read data from the data base
		;
		
		// debugging stuff
		if(TEST_MODE) {
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
// GetEvent
//---------------------------------
jerror_t DEventSourceEventStore::GetEvent(JEvent &event)
{
	// open the next file if available
	//if(!event_source)
	//	event_source = OpenNextFile();

	if(TEST_MODE) {
		// FOR DEBUGGING 
		// output some fake event with skim information
    	event.SetEventNumber(1);
    	event.SetRunNumber(10000);
    	event.SetJEventSource(this);
   		event.SetRef(NULL);
    	event.SetStatusBit(kSTATUS_FROM_FILE);
    	event.SetStatusBit(kSTATUS_PHYSICS_EVENT);

		for(int i=0; i<4; i++)
			if(gRandom->Uniform() < 0.5)
				event.SetStatusBit(BASE_SKIM_INDEX+i);

		return NOERROR;
	}
	
	// make sure the event source is open
	// if not, we're out of events
	if(!event_source) //throw RESOURCE_UNAVAILABLE;
		return EVENT_SOURCE_NOT_OPEN;
	
	// skip to next event
	// MoveToNextEvent();

	jerror_t retval = event_source->GetEvent(event);
	if(retval == NOERROR) {
		// tag event
	}
	
	return retval;
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
	JEventLoop* locEventLoop = event.GetJEventLoop();
    string dataClassName = factory->GetDataClassName();
	
	if (dataClassName =="DESSkimData") {
		JFactory<DESSkimData> *essd_factory = dynamic_cast<JFactory<DESSkimData>*>(factory);
		
		LockRead();			// LOCK class data
		vector<DESSkimData*> skim_data_vec;
		DESSkimData *skim_data = new DESSkimData(event, skim_list, BASE_SKIM_INDEX);
		skim_data_vec.push_back(skim_data);
		UnlockRead();       // UNLOCK

		essd_factory->CopyTo(skim_data_vec);
		return NOERROR;
	}
	
	if(!event_source) throw RESOURCE_UNAVAILABLE;
	
	return event_source->GetObjects(event, factory);
}

//---------------------------------
// PrintGrades
//---------------------------------
void DEventSourceEventStore::PrintGrades()
{
	vector<string> grades;
	esdb->GetGrades(grades);
	
	// print out information
	cout << endl << "Available grades:" << endl;
	for(vector<string>::iterator it = grades.begin();
		it != grades.end(); it++)
			cout << "  " << *it << endl;
	cout << endl;
}

//---------------------------------
// PrintSkims
//---------------------------------
void DEventSourceEventStore::PrintSkims(string datestamp, string grade)
{
	vector<string> skims;
	esdb->GetSkims(skims, datestamp, grade);
	
	// print out information
	cout << endl << "Available skims for grade " << grade << ":" << endl;
	for(vector<string>::iterator it = skims.begin();
		it != skims.end(); it++)
			cout << "  " << *it << endl;
	cout << endl;
}

