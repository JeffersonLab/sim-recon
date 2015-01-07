// $Id: JEventSourceGenerator_EVIO.cc 16643 2014-11-24 21:17:40Z davidl $
// $HeadURL: https://halldsvn.jlab.org/repos/branches/sim-recon-commissioning/src/programs/Utilities/plugins/DAQ/JEventSourceGenerator_EVIO.cc $
//
//    File: JEventSourceGenerator_EVIO.cc
// Created: Tue May 21 14:05:48 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//

#include <string>
using std::string;

#include "JEventSourceGenerator_EVIO.h"
using namespace jana;

#if HAVE_EVIO
#include <evioFileChannel.hxx>
#endif // HAVE_EVIO

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
#if HAVE_EVIO
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

#else  // HAVE_EVIO

	// EVIO support not enabled. We give a small probability
	// just so the JEventSource_EVIO constructor will get called
	// if not other JEventSource objects claim they can ready this
	// source. That way, the "you didn't compile in EVIO support"
	// message will get printed.
	return 1.0E-32;
#endif // HAVE_EVIO
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* JEventSourceGenerator_EVIO::MakeJEventSource(string source)
{
	return new JEventSource_EVIO(source.c_str());
}

