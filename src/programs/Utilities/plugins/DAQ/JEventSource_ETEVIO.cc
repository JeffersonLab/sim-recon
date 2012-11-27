// $Id$
//
//    File: JEventSource_ETEVIO.h
// Created: Mon Nov 26 10:48:42 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2  x86_64)
//

#include "JEventSource_ETEVIO.h"
using namespace jana;

#include <evioFileChannel.hxx>
using namespace evio;


//----------------
// Constructor
//----------------
JEventSource_ETEVIO::JEventSource_ETEVIO(const char* source_name):JEventSource_EVIO(source_name)
{
	// open event source here
	chan = new evioFileChannel(source_name);
	if(chan)chan->open();
}

//----------------
// ReadEVIOEvent
//----------------
jerror_t JEventSource_ETEVIO::ReadEVIOEvent(void)
{
	if(!chan->read())return NO_MORE_EVENTS_IN_SOURCE;

	return NOERROR;
}

//----------------
// Destructor
//----------------
JEventSource_ETEVIO::~JEventSource_ETEVIO()
{
	// chan is closed and deleted by base class destructor
}
