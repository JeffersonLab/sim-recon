// $Id$
//
//    File: DApplication.cc
// Created: Wed Jun  8 12:00:20 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <iostream>
#include <iomanip>
using namespace std;

#include <signal.h>

#include "DApplication.h"
#include "DEventProcessor.h"
#include "DEventSource.h"
#include "DEventLoop.h"
#include "DEvent.h"

void* LaunchThread(void* arg);

int SIGINT_RECEIVED = 0;
int NTHREADS_COMMAND_LINE = 0;

//-----------------------------------------------------------------
// ctrlCHandle
//-----------------------------------------------------------------
void ctrlCHandle(int x)
{
	SIGINT_RECEIVED++;
	cerr<<endl<<"SIGINT received ("<<SIGINT_RECEIVED<<")....."<<endl;
	if(SIGINT_RECEIVED == 3){
		cerr<<endl<<"Three SIGINTS received! Attempting graceful exit ..."<<endl<<endl;
	}
}


//---------------------------------
// DApplication    (Constructor)
//---------------------------------
DApplication::DApplication(int narg, char* argv[])
{
	// Set up to catch SIGINTs for graceful exits
	signal(SIGINT,ctrlCHandle);

	// Initialize application level mutexes
	pthread_mutex_init(&app_mutex, NULL);
	pthread_mutex_init(&current_source_mutex, NULL);

	// Variables used for calculating the rate
	show_ticker = 1;
	NEvents = 0;
	last_NEvents = 0;
	avg_NEvents = 0;
	avg_time = 0.0;
	rate_instantaneous = 0.0;
	rate_average = 0.0;
	
	// Sources
	current_source = NULL;
	for(int i=1; i<narg; i++){
		const char *arg="--nthreads=";
		if(!strncmp(arg, argv[i],strlen(arg))){
			NTHREADS_COMMAND_LINE = atoi(&argv[i][strlen(arg)]);
			continue;
		}
		if(argv[i][0] == '-')continue;
		source_names.push_back(argv[i]);
	}
}

//---------------------------------
// ~DApplication    (Destructor)
//---------------------------------
DApplication::~DApplication()
{
	
}

//---------------------------------
// NextEvent
//---------------------------------
derror_t DApplication::NextEvent(DEvent &event)
{
	/// Call the current source's GetEvent() method. If no
	/// source is open, or the current source has no more events,
	/// then open the next source and recall ourself

	pthread_mutex_lock(&current_source_mutex);
	if(!current_source){
		derror_t err = OpenNext();
		pthread_mutex_unlock(&current_source_mutex);
		if(err != NOERROR)return err;
		return NextEvent(event);
	}
	
	// Read next event from source. If none, then simply set current_source
	// to NULL and recall ourself
	switch(current_source->GetEvent(event)){
		case NO_MORE_EVENTS_IN_SOURCE:
			current_source = NULL;
			pthread_mutex_unlock(&current_source_mutex);
			return NextEvent(event);
			break;
		default:
			break;
	}

	pthread_mutex_unlock(&current_source_mutex);

	// Event counter
	NEvents++;

	return NOERROR;
}

//---------------------------------
// AddProcessor
//---------------------------------
derror_t DApplication::AddProcessor(DEventProcessor *processor)
{
	processor->SetDApplication(this);
	processors.push_back(processor);

	return NOERROR;
}

//---------------------------------
// RemoveProcessor
//---------------------------------
derror_t DApplication::RemoveProcessor(DEventProcessor *processor)
{
	vector<DEventProcessor*>::iterator iter = processors.begin();
	for(; iter!=processors.end(); iter++){
		if((*iter) == processor){
			processors.erase(iter);
			break;
		}
	}

	return NOERROR;
}

//---------------------------------
// GetProcessors
//---------------------------------
derror_t DApplication::GetProcessors(vector<DEventProcessor*> &processors)
{
	processors = this->processors;

	return NOERROR;
}

//---------------------------------
// AddDEventLoop
//---------------------------------
derror_t DApplication::AddDEventLoop(DEventLoop *loop)
{
	loops.push_back(loop);

	return NOERROR;
}

//---------------------------------
// RemoveDEventLoop
//---------------------------------
derror_t DApplication::RemoveDEventLoop(DEventLoop *loop)
{
	vector<DEventLoop*>::iterator iter = loops.begin();
	for(; iter!=loops.end(); iter++){
		if((*iter) == loop){
			loops.erase(iter);
			break;
		}
	}

	return NOERROR;
}

//---------------------------------
// GetDEventLoops
//---------------------------------
derror_t DApplication::GetDEventLoops(vector<DEventLoop*> &loops)
{
	loops = this->loops;

	return NOERROR;
}

//----------------
// LaunchThread
//----------------
void* LaunchThread(void* arg)
{
	/// This is a global function that is used to create
	/// a new DEventLoop object which runs in its own thread.

	// Create DEventLoop object. He automatically registers himself
	// with the DApplication object. 
	DEventLoop *eventLoop = new DEventLoop((DApplication*)arg);

	// Loop over events until done. Catch any derror_t's thrown
	try{
		eventLoop->Loop();
	}catch(derror_t err){
		cerr<<__FILE__<<":"<<__LINE__<<" ERROR ("<<err<<") caught for thread "<<pthread_self()<<endl;
	}

	// Delete DEventLoop object. He automatically de-registers himself
	// with the DEventLoop Object
	delete eventLoop;

	return arg;
}

//---------------------------------
// Init
//---------------------------------
derror_t DApplication::Init(void)
{
	// Call init Processors (note: factories don't exist yet)
	try{
		for(unsigned int i=0;i<processors.size();i++)processors[i]->init();
	}catch(derror_t err){
		cerr<<endl;
		cerr<<__FILE__<<":"<<__LINE__<<" Error thrown ("<<err<<") from DEventProcessor::init()"<<endl;
		exit(-1);
	}

	return NOERROR;
}

//---------------------------------
// Run
//---------------------------------
derror_t DApplication::Run(DEventProcessor *proc, int Nthreads)
{
	// If a DEventProcessor was passed, then add it to our list first
	if(proc){
		derror_t err = AddProcessor(proc);
		if(err)return err;
	}

	// Call init() routines
	Init();

	// Launch all threads
	if(Nthreads<1){
		// If Nthreads is less than 1 then automatically set to 1
		Nthreads = 1;
	}
	if(NTHREADS_COMMAND_LINE>0){
		Nthreads = NTHREADS_COMMAND_LINE;
	}
	cout<<"Launching threads "; cout.flush();
	for(int i=0; i<Nthreads; i++){
		pthread_t thr;
		pthread_create(&thr, NULL, LaunchThread, this);
		cout<<".";cout.flush();
	}
	cout<<endl;
	
	// Do a sleepy loop so the threads can do their work
	struct timespec req, rem;
	req.tv_nsec = (int)0.5E9; // set to 1/2 second
	req.tv_sec = 0;
	double sleep_time = (double)req.tv_sec + (1.0E-9)*(double)req.tv_nsec;
	do{
		// Sleep for a specific amount of time and calculate the rate
		// on each iteration through the loop
		rem.tv_sec = rem.tv_nsec = 0;
		nanosleep(&req, &rem);
		if(rem.tv_sec == 0 && rem.tv_nsec == 0){
			// If there was no time remaining, then we must have slept
			// the whole amount
			int delta_NEvents = NEvents - last_NEvents;
			avg_NEvents += delta_NEvents;
			avg_time += sleep_time;
			rate_instantaneous = (double)delta_NEvents/sleep_time;
			rate_average = (double)avg_NEvents/avg_time;
		}else{
			cout<<__FILE__<<":"<<__LINE__<<" didn't sleep full "<<sleep_time<<" seconds!"<<endl;
		}
		last_NEvents = NEvents;
		
		// If show_ticker is set, then update the screen with the rate(s)
		if(show_ticker && loops.size()>0)PrintRate();
		
		if(SIGINT_RECEIVED)Quit();

		// When a DEventLoop runs out of events, it removes itself from
		// the list before returning from the thread.
	}while(loops.size() > 0);
	
	// Call erun() and fini() methods and delete event sources
	Fini();
	
	cout<<" "<<NEvents<<" events processed. Average rate: "
		<<Val2StringWithPrefix(rate_average)<<"Hz"<<endl;

	return NOERROR;
}

//---------------------------------
// Fini
//---------------------------------
derror_t DApplication::Fini(void)
{
	// Make sure erun is called
	for(unsigned int i=0;i<processors.size();i++){
		DEventProcessor *proc = processors[i];
		if(proc->brun_was_called() && !proc->erun_was_called()){
			try{
				proc->erun();
			}catch(derror_t err){
				cerr<<endl;
				cerr<<__FILE__<<":"<<__LINE__<<" Error thrown ("<<err<<") from DEventProcessor::erun()"<<endl;
			}
			proc->Set_erun_called();
		}
	}

	// Call fini Processors
	try{
		for(unsigned int i=0;i<processors.size();i++)processors[i]->fini();
	}catch(derror_t err){
		cerr<<endl;
		cerr<<__FILE__<<":"<<__LINE__<<" Error thrown ("<<err<<") from DEventProcessor::fini()"<<endl;
	}
	
	// Delete all sources allowing them to close cleanly
	for(unsigned int i=0;i<sources.size();i++)delete sources[i];
	sources.clear();

	return NOERROR;
}

//---------------------------------
// Pause
//---------------------------------
void DApplication::Pause(void)
{
	vector<DEventLoop*>::iterator iter = loops.begin();
	for(; iter!=loops.end(); iter++){
		(*iter)->Pause();
	}
}

//---------------------------------
// Resume
//---------------------------------
void DApplication::Resume(void)
{
	vector<DEventLoop*>::iterator iter = loops.begin();
	for(; iter!=loops.end(); iter++){
		(*iter)->Resume();
	}
}

//---------------------------------
// Quit
//---------------------------------
void DApplication::Quit(void)
{
	cout<<endl<<"Telling all threads to quit ..."<<endl;
	vector<DEventLoop*>::iterator iter = loops.begin();
	for(; iter!=loops.end(); iter++){
		(*iter)->Quit();
	}
}

//----------------
// Val2StringWithPrefix
//----------------
string DApplication::Val2StringWithPrefix(float val)
{
	char *units = "";
	if(val>1.5E9){
		val/=1.0E9;
		units = "G";
	}else 	if(val>1.5E6){
		val/=1.0E6;
		units = "M";
	}else if(val>1.5E3){
		val/=1.0E3;
		units = "k";
	}else if(val<1.0E-7){
		units = "";
	}else if(val<1.0E-4){
		val/=1.0E6;
		units = "u";
	}else if(val<1.0E-1){
		val/=1.0E3;
		units = "m";
	}
	
	char str[256];
	sprintf(str,"%3.1f%s", val, units);

	return string(str);
}
	
//----------------
// PrintRate
//----------------
void DApplication::PrintRate(void)
{
	string event_str = Val2StringWithPrefix(NEvents) + " events";
	string ir_str = Val2StringWithPrefix(rate_instantaneous) + "Hz";
	string ar_str = Val2StringWithPrefix(rate_average) + "Hz";
	cout<<"  "<<event_str<<"   "<<ir_str<<"  (average rate: "<<ar_str<<")   \r";
	cout.flush();
}

//---------------------------------
// OpenNext
//---------------------------------
derror_t DApplication::OpenNext(void)
{
	/// Open the next source in the list. If there are none,
	/// then return NO_MORE_EVENT_SOURCES
	
	if(sources.size() >= source_names.size())return NO_MORE_EVENT_SOURCES;
	
	// Get the source type for the next source name from the command line
	// Instantiate an object of the appropriate type.
	const char *sname = source_names[sources.size()];
	const string type = DEventSource::GuessSourceType(sname);
	current_source = DEventSource::OpenSource(sname);
	if(!current_source){
		cerr<<"Unknown source type \""<<sname<<"\""<<endl;
	}
	
	// Add source to list (even if it's NULL!)
	sources.push_back(current_source);
	
	return NOERROR;
}

