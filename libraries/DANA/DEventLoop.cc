// Author: David Lawrence  June 24, 2004
//
//
// DEventLoop methods
//

#include <iostream>
using namespace std;

#include "DEventLoop.h"
#include "DEventSourceET.h"

#include "../BCAL/DBCAL.h"
#include "../CDC/DCDC.h"
#include "../FCAL/DFCAL.h"
#include "../TAGGER/DTagger.h"
#include "../TOF/DTOF.h"
#include "../UPV/DUPV.h"
#include "../TRIGGER/DTRIGGER.h"
#include "../VERTEX/DVERTEX.h"
#include "../FDC/DFDC.h"
#include "../CHERENKOV/DCHERENKOV.h"

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
		case DEventSource::EVENT_SOURCE_FILE:
			source = new DEventSourceFILE(narg, argv);
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
	// At this point, all of the "user" specified event processors
	// are in the processors list. However, we want the "standard"
	// processors to be at the top of the list so they are called
	// first allowing the "users" processors to see the reconstructed
	// data. Thus, we must make a copy of the current list, make
	// the list empty, and then append the user processors after
	// the standard ones have been added.
	DEventProcessor* processors_copy[MAX_EVENT_PROCESSORS];
	int Nprocessors_copy = 0;
	for(int i=0;i<Nprocessors;i++)processors_copy[Nprocessors_copy++] = processors[i];
	Nprocessors = 0;

	// Add standard processors specified by proc_mask
	if(proc_mask&EVENT_PROCESSOR_TRIGGER	)AddProcessor(new DTRIGGER);
	if(proc_mask&EVENT_PROCESSOR_TAGGER		)AddProcessor(new DTagger);
	if(proc_mask&EVENT_PROCESSOR_UPV			)AddProcessor(new DUPV);
	if(proc_mask&EVENT_PROCESSOR_VERTEX		)AddProcessor(new DVERTEX);
	if(proc_mask&EVENT_PROCESSOR_BCAL		)AddProcessor(new DBCAL);
	if(proc_mask&EVENT_PROCESSOR_CDC			)AddProcessor(new DCDC);
	if(proc_mask&EVENT_PROCESSOR_FDC			)AddProcessor(new DFDC);
	if(proc_mask&EVENT_PROCESSOR_CHERENKOV	)AddProcessor(new DCHERENKOV);
	if(proc_mask&EVENT_PROCESSOR_TOF			)AddProcessor(new DTOF);
	if(proc_mask&EVENT_PROCESSOR_FCAL		)AddProcessor(new DFCAL);
	
	// Re-append "user" processors
	for(int i=0;i<Nprocessors_copy;i++)AddProcessor(processors_copy[i++]);
	
	// Loop over events
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

