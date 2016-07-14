// $Id$
//
//    File: DEventSourceEventStoreGenerator.cc
// Creator: sdobbs
//

#include "DEventSourceEventStoreGenerator.h"
#include "DEventSourceEventStore.h"
#include "DFactoryGenerator_DESSkimData.h"

/*
// Make this a plugin
extern "C" {
  void InitPlugin(JApplication *app) {
    InitJANAPlugin(app);
    app->AddEventSourceGenerator(new DEventSourceEventStoreGenerator());
    app->AddFactoryGenerator(new DFactoryGenerator_DESSkimData());
  }
} // "extern C"
*/

//---------------------------------
// DEventSourceEventStoreGenerator    (Constructor)
//---------------------------------
DEventSourceEventStoreGenerator::DEventSourceEventStoreGenerator()
{

}

//---------------------------------
// ~DEventSourceEventStoreGenerator    (Destructor)
//---------------------------------
DEventSourceEventStoreGenerator::~DEventSourceEventStoreGenerator()
{

}

//---------------------------------
// Description
//---------------------------------
const char* DEventSourceEventStoreGenerator::Description(void)
{
	return "EventStore";
}

//---------------------------------
// CheckOpenable
//---------------------------------
double DEventSourceEventStoreGenerator::CheckOpenable(string source)
{
	// See if this looks like the beginning of an eventstore query 
	return source.substr(0,10)=="eventstore" ? 0.5:0.0;
}

//---------------------------------
// MakeJEventSource
//---------------------------------
JEventSource* DEventSourceEventStoreGenerator::MakeJEventSource(string source)
{
	return new DEventSourceEventStore(source.c_str());
}
