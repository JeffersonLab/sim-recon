// $Id$
//
//    File: JEventSourceGenerator_FileEVIO.cc
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#include <string>
using std::string;

#include "JEventSourceGenerator_FileEVIO.h"
using namespace jana;

#include <evioFileChannel.hxx>
#include <evioUtil.hxx>
using namespace evio;



//---------------------------------
// Description
//---------------------------------
const char* JEventSourceGenerator_FileEVIO::Description(void)
{
	return "FileEVIO - Reads EVIO formatted files";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double JEventSourceGenerator_FileEVIO::CheckOpenable(string source)
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
		
		return 0.0;
	}	
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* JEventSourceGenerator_FileEVIO::MakeJEventSource(string source)
{
	return new JEventSource_FileEVIO(source.c_str());
}

