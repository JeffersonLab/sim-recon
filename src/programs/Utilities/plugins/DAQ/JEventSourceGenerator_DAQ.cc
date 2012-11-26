// $Id$
//
//    File: JEventSourceGenerator_DAQ.cc
// Created: Tue Aug  7 15:22:29 EDT 2012
// Creator: davidl (on Darwin harriet.jlab.org 11.4.0 i386)
//

#include <string>
using std::string;

#include "JEventSourceGenerator_DAQ.h"
#include "JEventSourceGenerator_ETEVIO.h"
#include "JFactoryGenerator_DAQ.h"
using namespace jana;

#include <evioFileChannel.hxx>
#include <evioUtil.hxx>
using namespace evio;

// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
extern "C"{
	void InitPlugin(JApplication *app){
		InitJANAPlugin(app);
		app->AddEventSourceGenerator(new JEventSourceGenerator_DAQ());
		app->AddEventSourceGenerator(new JEventSourceGenerator_ETEVIO());
		app->AddFactoryGenerator(new JFactoryGenerator_DAQ());
	}
} // "C"


//---------------------------------
// Description
//---------------------------------
const char* JEventSourceGenerator_DAQ::Description(void)
{
	return "DAQ";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double JEventSourceGenerator_DAQ::CheckOpenable(string source)
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
		
	} catch (evioException *e) {
		
		// Could not open file
		return 0.0;
	}	
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* JEventSourceGenerator_DAQ::MakeJEventSource(string source)
{
	return new JEventSource_DAQ(source.c_str());
}

