// $Id$
//
//    File: DEventSourceEventStore.cc
// Creator: sdobbs 
//

#include <iostream>
#include <string>
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
	esdb_connection = "mysql://eventstore@hallddb.jlab.org/esdb";    // default to main JLab ES server
	
	// First, parse the eventsource query
	// For details of the query format, see: <...>
	
	// Tokenize the query string
	string es_query(source_name);
	istringstream iss(es_query);
	vector<string> tokens;
	copy(istream_iterator<string>(iss),
    	istream_iterator<string>(),
     	back_inserter(tokens));
	
	// Check query header
	if( (tokens[0] != "eventstore") || (tokens[1] != "in") )
		throw JException("Invalid ES query = " + es_query + "\n\n" + EventstoreQueryHelp() );
		


	// initialize database connection
	if(esdb_connection.substr(0,8) == "mysql://") {
		cout << "MySQL connection" << endl;
	} else if(esdb_connection.substr(0,8) == "sqlite://") {
		cout << "SQLite connection" << endl;	
	} 
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

