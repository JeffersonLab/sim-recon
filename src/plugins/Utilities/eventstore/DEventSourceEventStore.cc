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

#include "DEventSourceEventStore.h"

static string EventstoreQueryHelp() {
	string str = "this is a help string";
	return str;
}

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
	// skip to next event
	if(!event_source) throw RESOURCE_UNAVAILABLE;

	return event_source->GetEvent(event);
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

