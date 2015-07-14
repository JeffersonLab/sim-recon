// $Id$
//
//    File: DEventSourceEVIOGenerator.cc
// Created: Sat May  8 13:54:46 EDT 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//

#include "DEventSourceEVIOGenerator.h"
#include "DEventSourceEVIO.h"

// Make this a plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddEventSourceGenerator(new DEventSourceEVIOGenerator());
  }
} // "extern C"


//---------------------------------
// DEventSourceEVIOGenerator    (Constructor)
//---------------------------------
DEventSourceEVIOGenerator::DEventSourceEVIOGenerator()
{

}

//---------------------------------
// ~DEventSourceEVIOGenerator    (Destructor)
//---------------------------------
DEventSourceEVIOGenerator::~DEventSourceEVIOGenerator()
{

}

//---------------------------------
// Description
//---------------------------------
const char* DEventSourceEVIOGenerator::Description(void)
{
	return "EVIO";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double DEventSourceEVIOGenerator::CheckOpenable(string source)
{
	return source.find(".evio",0)==string::npos ? 0.0:0.5;
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* DEventSourceEVIOGenerator::MakeJEventSource(string source)
{
	return new DEventSourceEVIO(source.c_str());
}
