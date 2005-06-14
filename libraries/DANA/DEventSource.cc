// $Id$
//
//    File: DEventSource.cc
// Created: Wed Jun  8 12:31:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DEventSource.h"
#include "DEventSourceHDDM.h"

//---------------------------------
// DEventSource    (Constructor)
//---------------------------------
DEventSource::DEventSource(const char *source_name)
{
	this->source_name = source_name;
	source_is_open = 0;
	pthread_mutex_init(&read_mutex, NULL);
	event_buffer_size = 10;
	Nevents_read = 0;
}

//---------------------------------
// ~DEventSource    (Destructor)
//---------------------------------
DEventSource::~DEventSource()
{

}

//----------------
// GuessSourceType
//----------------
const char* DEventSource::GuessSourceType(const char* source_name)
{
	/// Guesses type of event source
	/// Only hddm is supported for now.
	return DEventSourceHDDM::static_className();
}

//----------------
// OpenSource
//----------------
DEventSource* DEventSource::OpenSource(const char* name)
{
	// Get the name of the subclass of DEventSource for this source
	string type = GuessSourceType(name);

	// Create a new object of the appropriate class and return a pointer to it
	if(type == DEventSourceHDDM::static_className())	return new DEventSourceHDDM(name);

	// Mmmm... Don't seem to know this source type.
	cerr<<__FILE__<<":"<<__LINE__<<" Unknown source type for \""<<name<<"\"."<<endl;
	return NULL;
}

//----------------
// GetObjects
//----------------
derror_t DEventSource::GetObjects(const char *name, vector<void*> &v, const char* tag, void* ref, DFactory_base *factory)
{
	/// This only gets called if the subclass doesn't define it. In that
	/// case, the subclass must not support objects.
	return OBJECT_NOT_AVAILABLE;
}

