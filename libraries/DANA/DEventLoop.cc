	// Author: David Lawrence  June 24, 2004
//
//
// DEventLoop methods
//
#include <stdio.h>
#include <iostream>
using namespace std;

#include "DEventLoop.h"
#include "DEventProcessor.h"
#include "DEventSourceET.h"
#include "DEventSourceFile.h"


//----------------
// Constructor
//----------------
DEventLoop::DEventLoop(int narg, char *argv[])
{
	/// Constructor for DEventLoop object 

	// Initialize variables
	Nprocessors = 0;
	last_print_rate_time = 0;
	quit=0;
	goto_event = -1;
	
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
// Init
//----------------
derror_t DEventLoop::Init(void)
{

	// Let each detector system install factories
	extern derror_t BCAL_init(DEvent*);
	extern derror_t CDC_init(DEvent*);
	extern derror_t CHERENKOV_init(DEvent*);
	extern derror_t FCAL_init(DEvent*);
	extern derror_t FDC_init(DEvent*);
	extern derror_t TAGGER_init(DEvent*);
	extern derror_t TOF_init(DEvent*);
	extern derror_t TRIGGER_init(DEvent*);
	extern derror_t UPV_init(DEvent*);
	extern derror_t TRACKING_init(DEvent*);
	BCAL_init(this);
	CDC_init(this);
	CHERENKOV_init(this);
	FCAL_init(this);
	FDC_init(this);
	TAGGER_init(this);
	TOF_init(this);
	TRIGGER_init(this);
	UPV_init(this);
	TRACKING_init(this);
	
	// Set the event_loop field for all processors
	for(int i=0;i<Nprocessors;i++)processors[i]->eventLoop = this;
	
	// Call init Processors
	for(int i=0;i<Nprocessors;i++)processors[i]->init();

	return NOERROR;
}

//----------------
// OneEvent
//----------------
derror_t DEventLoop::OneEvent(void)
{
	/// Read in and process one event. If eventno is
	/// less than 0, then grab the next event from
	/// the source. Otherwise, jump to the specified event

	// Clear evnt_called flag in all factories
	ClearFactories();

	derror_t err;
	if(goto_event>=0){
		err = source->GotoEvent(goto_event);
		if(!err)_eventnumber = goto_event - 1; // eventnumber is incremented below
		goto_event = -1;
	}else{
		err = source->NextEvent();
	}
	
	switch(err){
		case NOERROR:
			break;
		case NO_MORE_EVENT_SOURCES:
			cout<<endl<<"No more event sources"<<endl;
			break;
		case EVENT_NOT_IN_MEMORY:
			cout<<endl<<"Event not in memory"<<endl;
			_eventnumber--;
			break;
		default:
			break;
	}
	if(err != NOERROR && err !=EVENT_NOT_IN_MEMORY)return err;
	
	// Copy pointer to hddm_s to our object (from DEvent inheritance)
	_hddm_s = source->hddm_s;
	
	// Need to extract the event number and run number here.
	//runnumber = 1;
	_eventnumber++;
	
	// Call Event Processors
	DEventProcessor **p = processors;
	for(int i=0;i<Nprocessors;i++, p++){
		(*p)->hddm_s = hddm_s();

		// Call brun routine if run number has changed or it's not been called
		if(_runnumber!=(*p)->GetBRUN_RunNumber()){
			if((*p)->brun_was_called() && !(*p)->erun_was_called()){
				(*p)->erun();
				(*p)->Set_erun_called();
			}
			(*p)->Clear_brun_called();
		}
		if(!(*p)->brun_was_called()){
			(*p)->brun(_runnumber);
			(*p)->Set_brun_called();
			(*p)->Clear_erun_called();
			(*p)->SetBRUN_RunNumber(_runnumber);
		}

		// Call the event routine
		derror_t err = (*p)->evnt(_eventnumber);

		// More will need to be done to truly filter out the event.
		// for now, this just provides a means to not completely process
		// events by adding a "filter" to the front of the processor list
		if(err == FILTER_EVENT_OUT)break;
	}
	
	// In order for the source to jump to an event via GotoEvent(),
	// it needs to know the event number which may only be obtained
	// after parsing the data.
	source->RecordEventNumber(_eventnumber);
		
	return NOERROR;
}

//----------------
// Fini
//----------------
derror_t DEventLoop::Fini(void)
{
	// Call fini Processors
	for(int i=0;i<Nprocessors;i++)processors[i]->fini();
	
	return NOERROR;
}

//----------------
// Run
//----------------
derror_t DEventLoop::Run(void)
{
	// Initialize
	derror_t err;
	err = Init();
	if(err!=NOERROR)return err;

	// Loop over events
	do{
		err=OneEvent();
		if(err!=NOERROR)break;
	}while(!quit);

	// Clean up
	err = Fini();
	
	cout<<endl<<source->Nevents_read<<" events processed total"<<endl;

	return err;
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
			rateunits = "MHz";
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

//----------------
// GetSourceName
//----------------
const char* DEventLoop::GetSourceName(void)
{
	return source->GetSourceName();
}

