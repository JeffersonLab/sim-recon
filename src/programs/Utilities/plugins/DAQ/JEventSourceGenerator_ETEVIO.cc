// $Id$
//
//    File: JEventSourceGenerator_ETEVIO.cc
// Created: Mon Nov 26 11:01:01 EST 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#include <string>
using std::string;

#include "JEventSourceGenerator_ETEVIO.h"
using namespace jana;

#include <evioFileChannel.hxx>
#include <evioUtil.hxx>
using namespace evio;


//---------------------------------
// Description
//---------------------------------
const char* JEventSourceGenerator_ETEVIO::Description(void)
{
	return "ETEVIO - Reads EVIO from ET systems";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double JEventSourceGenerator_ETEVIO::CheckOpenable(string source)
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
		
		return 0.1;
		
	} catch (evioException *e) {
		
		// Could not open file
		return 0.0;
	}	
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* JEventSourceGenerator_ETEVIO::MakeJEventSource(string source)
{
	return new JEventSource_ETEVIO(source.c_str());
}

