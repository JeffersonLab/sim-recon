// Author: David Lawrence  June 24, 2004
//
//
// DEventLoop methods
//

#include <iostream>
using namespace std;

#include "DEventLoop.h"
#include "DEventSourceET.h"

//----------------
// Constructor
//----------------
DEventLoop::DEventLoop(int narg, char *argv[])
{
	/// Constructor for DEventLoop object 

	// Initialize variables
	proc_mask = EVENT_PROCESSOR_ALL;
	Nprocessors = 0;
	
	// Determine event source type
	DEventSource::EVENT_SOURCE_TYPE source_type;
	switch(DEventSource::GuessSourceType(narg, argv)){
		case DEventSource::EVENT_SOURCE_ET:
			source = new DEventSourceET(narg, argv);
			break;
		default:
			source = NULL;
	}
}

//----------------
// Destructor
//----------------
DEventLoop::~DEventLoop()
{
	if(source)delete source;
}

//----------------
// AddProcessor
//----------------
derror_t DEventLoop::AddProcessor(DEventProcessor *processor)
{
	if(Nprocessors >= MAX_EVENT_PROCESSORS){
		cerr<<"Maximum number of event processors defined in DEventLoop"<<endl;
		cerr<<"object ("<<MAX_EVENT_PROCESSORS<<")"<<endl;
		return MAX_EVENT_PROCESSORS_EXCEEDED;
	}
	
	processors[Nprocessors++] = processor;
}

//----------------
// SetStandardProcessors
//----------------
derror_t DEventLoop::SetStandardProcessors(EVENT_PROCESSORS processmask)
{
	proc_mask = processmask;
	
	return NOERROR;
}

//----------------
// Run
//----------------
derror_t DEventLoop::Run(void)
{
	// Add standard processors specified by proc_mask
}

//----------------
// Run
//----------------
derror_t DEventLoop::Run(DEventProcessor *proc)
{
	derror_t err;
	if(err=AddProcessor(proc))return err;
	
	return Run();
}

//----------------
// Run
//----------------
derror_t DEventLoop::Run(EVENT_PROCESSORS mask)
{
	derror_t err;
	if(err=SetStandardProcessors(mask))return err;
	
	return Run();
}

//----------------
// Run
//----------------
derror_t DEventLoop::Run(DEventProcessor *proc, EVENT_PROCESSORS mask)
{
	derror_t err;
	if(err=AddProcessor(proc))return err;
	if(err=SetStandardProcessors(mask))return err;
	
	return Run();
}

