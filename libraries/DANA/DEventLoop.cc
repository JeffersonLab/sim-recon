// Author: David Lawrence  June 24, 2004
//
//
// DEventLoop methods
//
#include <stdio.h>
#include <iostream>
using namespace std;

#include "DEventLoop.h"
#include "DEventSourceET.h"
#include "DEventSourceFile.h"

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

int DONE = 0;

//----------------
// Constructor
//----------------
DEventLoop::DEventLoop(int narg, char *argv[])
{
	/// Constructor for DEventLoop object 

	// Initialize variables
	proc_mask = EVENT_PROCESSOR_ALL;
	Nprocessors = 0;
	last_print_rate_time = 0;
	
	// Determine event source type
	DEventSource::EVENT_SOURCE_TYPE source_type;
	switch(DEventSource::GuessSourceType(narg, argv)){
		case DEventSource::EVENT_SOURCE_ET:
			source = new DEventSourceET(narg, argv);
			break;
		case DEventSource::EVENT_SOURCE_FILE:
			source = new DEventSourceFile(narg, argv);
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
	
	return NOERROR;
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
	
	// Set the event_loop field for all processors
	for(int i=0;i<Nprocessors;i++)processors[i]->event_loop = this;
	
	// Loop over events
	do{
		derror_t err = source->NextEvent();
		switch(err){
			case NOERROR:
				break;
			case NO_MORE_EVENT_SOURCES:
				cout<<endl<<"No more event sources"<<endl;
				break;
			default:
				break;
		}
		if(err != NOERROR)break;
		
		// Call Event Processors
		DEventProcessor **p = processors;
		for(int i=0;i<Nprocessors;i++, p++){
			(*p)->hddm = source->hddm;
			(*p)->evnt(0);
		}
		
	}while(!DONE);
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

//----------------
// GetRate
//----------------
float DEventLoop::GetRate(void)
{
	/// This just returns GetRate from the private DEventSource member
	if(!source){
		cerr<<__FILE__<<":"<<__LINE__<<" No event source set (source=NULL)!"<<endl;
		return -1.0;
	}
	
	return source->prate_last_rate;
}

//----------------
// PrintRate
//----------------
derror_t DEventLoop::PrintRate(void)
{
	/// Update the screen with the current event processing rate.
	/// This governs itself to only update as often as the rate
	/// is recalculated by DEventSource so it is safe to call
	/// it every event.
	if(!source){
		cerr<<__FILE__<<":"<<__LINE__<<" No event source set (source=NULL)!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	//cout<<__FILE__<<":"<<__LINE__<<" last_print_rate_time="<<last_print_rate_time;
	//cout<<" source->prate_last_time="<<source->prate_last_time<<endl;

	if(last_print_rate_time != source->prate_last_time){
		float rate = source->prate_last_rate;
		char *rateunits = "Hz";
		if(rate>1500000.){
			rate/=1000000.;
			rateunits = "kHz";
		}else if(rate>1500.){
			rate/=1000;
			rateunits = "kHz";
		}
		float Nevents = (float)source->Nevents_read;
		char *eventunits = "";
		if(Nevents>10000.){
			Nevents/=1000.;
			eventunits = "k";
		}
		
		printf("  %3.1f%s events  (%3.1f%s)          \r", Nevents,eventunits, rate, rateunits);
		fflush(stdout);
		last_print_rate_time = source->prate_last_time;
	}
	
	return NOERROR;
}

