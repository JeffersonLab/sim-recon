// $Id$
//
//    File: JEventSource_FileEVIO.h
// Created: Tue Nov 27 08:12:29 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2  x86_64)
//

#include "JEventSource_FileEVIO.h"
using namespace jana;

#include <evioFileChannel.hxx>
using namespace evio;

//----------------
// Constructor
//----------------
JEventSource_FileEVIO::JEventSource_FileEVIO(const char* source_name):JEventSource_EVIO(source_name)
{
	// open event source here
	chan = new evioFileChannel(source_name);
	if(chan)chan->open();
}

//----------------
// ReadEVIOEvent
//----------------
jerror_t JEventSource_FileEVIO::ReadEVIOEvent(void)
{
	if(!chan->read())return NO_MORE_EVENTS_IN_SOURCE;
	
	return NOERROR;
}

//----------------
// Destructor
//----------------
JEventSource_FileEVIO::~JEventSource_FileEVIO()
{
	// chan is closed and deleted by base class destructor
}
