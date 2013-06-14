// $Id$
//
//    File: JEventSourceGenerator_EVIO.cc
// Created: Tue May 21 14:05:48 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#include <string>
using std::string;

#include "JEventSourceGenerator_EVIO.h"
using namespace jana;

#include <evioFileChannel.hxx>


//---------------------------------
// Description
//---------------------------------
const char* JEventSourceGenerator_EVIO::Description(void)
{
	return "EVIO  - Reads EVIO formatted data from file or ET system";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double JEventSourceGenerator_EVIO::CheckOpenable(string source)
{
	// This should return a value between 0 and 1 inclusive
	// with 1 indicating it definitely can read events from
	// the specified source and 0 meaning it definitely can't.
	// Typically, this will just check the file suffix.

	// Try to open the file
	try {
		
		// create evio file channel object using first arg as file name
		evioFileChannel chan(source);
		
		// open the file. Throws exception if not successful
		chan.open();
		
		// close file
		chan.close();
		
		return 0.5;
		
	} catch (evioException &e) {

		// Could not open file. Check if name starts with "ET:"
		if(source.substr(0,3) == "ET:")return 0.1;
	}	

	// Doesn't seem to be a source we can open
	return 0.0;
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* JEventSourceGenerator_EVIO::MakeJEventSource(string source)
{
	return new JEventSource_EVIO(source.c_str());
}

