// $Id$
//
//    File: DEventLoop.cc
// Created: Wed Jun  8 12:30:51 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include <cstdio>
#include <iostream>
#include <iomanip>
using namespace std;

#include <signal.h>
#include <setjmp.h>
#include <unistd.h>

#include "DApplication.h"
#include "DEventLoop.h"
#include "DEvent.h"
#include "DFactory.h"

jmp_buf SETJMP_ENV;

// Thread commits suicide when it receives HUP signal
void thread_HUP_sighandler(int sig)
{
	cerr<<"Caught HUP signal for thread 0x"<<hex<<pthread_self()<<dec<<" thread exiting..."<<endl;
	pthread_exit(NULL);
}

//---------------------------------
// DEventLoop    (Constructor)
//---------------------------------
DEventLoop::DEventLoop(DApplication *app)
{
	// Last Resort exit strategy: If this thread stops responding, the
	// main thread will send it a HUP signal to tell it to exit immediately.
	signal(SIGHUP, thread_HUP_sighandler);

	this->app = app;
	app->AddDEventLoop(this, heartbeat);
	event.SetDEventLoop(this);
	pause = 0;
	quit = 0;
	auto_free = 1;
	pthread_id = pthread_self();
	
	// Let each detector system install factories
	extern derror_t BCAL_init(DEventLoop*);
	extern derror_t CDC_init(DEventLoop*);
	extern derror_t CHERENKOV_init(DEventLoop*);
	extern derror_t FCAL_init(DEventLoop*);
	extern derror_t FDC_init(DEventLoop*);
	extern derror_t TAGGER_init(DEventLoop*);
	extern derror_t TOF_init(DEventLoop*);
	extern derror_t TRIGGER_init(DEventLoop*);
	extern derror_t UPV_init(DEventLoop*);
	extern derror_t TRACKING_init(DEventLoop*);
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
	
	// Copy the event processor list to our local vector
	app->GetProcessors(processors);

	dparms.PrintParameters();
}

//---------------------------------
// ~DEventLoop    (Destructor)
//---------------------------------
DEventLoop::~DEventLoop()
{
	// Remove us from the DEventLoop's list. The application exits
	// when there are no more DEventLoops registered with it.
	app->RemoveDEventLoop(this);

	// Call all factories' fini methods
	for(unsigned int i=0; i<factories.size(); i++){
		try{
			factories[i]->fini();
		}catch(derror_t err){
			cerr<<endl;
			cerr<<__FILE__<<":"<<__LINE__<<" Error thrown ("<<err<<") from DFactory<";
			cerr<<factories[i]->dataClassName()<<">::fini()"<<endl;
		}
	}

	// Delete all of the factories
	for(unsigned int i=0; i<factories.size(); i++){
		try{
			delete factories[i];
		}catch(derror_t err){
			cerr<<endl;
			cerr<<__FILE__<<":"<<__LINE__<<" Error thrown ("<<err<<") while deleting DFactory<";
			cerr<<factories[i]->dataClassName()<<">"<<endl;
		}
	}

	factories.clear();
}

//-------------
// AddFactory
//-------------
derror_t DEventLoop::AddFactory(DFactory_base* factory)
{
	factory->SetDEventLoop(this);
	factory->SetDApplication(app);
	factories.push_back(factory);

	return NOERROR;
}

//-------------
// RemoveFactory
//-------------
derror_t DEventLoop::RemoveFactory(DFactory_base* factory)
{
	vector<DFactory_base*>::iterator iter = factories.begin();
	for(;iter!=factories.end(); iter++){
		if(*iter == factory){
			factories.erase(iter);
			break;
		}
	}

	return NOERROR;
}

//-------------
// GetFactory
//-------------
DFactory_base* DEventLoop::GetFactory(const string data_name, const char *tag)
{
	// Search for specified factory and return pointer to it
	vector<DFactory_base*>::iterator iter = factories.begin();
	for(; iter!=factories.end(); iter++){
		if(data_name == (*iter)->dataClassName()){
			if(!strcmp((*iter)->Tag(), tag)){
				return *iter;
			}
		}
	}

	// No factory found. Return NULL
	return NULL;
}

//-------------
// GetFactoryNames
//-------------
vector<string> DEventLoop::GetFactoryNames(void)
{
	/// Return a vector<string> whose members are 
	/// the names of the currently registered factories. 
	vector<string> names;
	vector<DFactory_base*>::iterator factory = factories.begin();
	for(; factory!=factories.end(); factory++){
		names.push_back((*factory)->dataClassName());
	}	
	
	return names;
}

//-------------
// ClearFactories
//-------------
derror_t DEventLoop::ClearFactories(void)
{
	/// Loop over all factories and call their Reset() methods.
	/// Amoung other things, this will clear their evnt_called flags
	/// This is called from DEventLoop at the
	/// begining of a new event.

	vector<DFactory_base*>::iterator iter = factories.begin();
	for(; iter!=factories.end(); iter++){
		(*iter)->Reset();
	}

	return NOERROR;
}

//-------------
// PrintFactories
//-------------
derror_t DEventLoop::PrintFactories(int sparsify)
{
	/// Print a list of all registered factories to the screen
	/// along with a little info about each.

	cout<<endl;
	cout<<"Registered factories: ("<<factories.size()<<" total)"<<endl;
	cout<<endl;
	cout<<"Name:             nrows:  tag:"<<endl;
	cout<<"---------------- ------- --------------"<<endl;

	for(unsigned int i=0; i<factories.size(); i++){
		DFactory_base *factory = factories[i];
		
		if(sparsify)
			if(factory->GetNrows()<1)continue;
		
		// To make things look pretty, copy all values into the buffer "str"
		string str(79,' ');
		string name = factory->dataClassName();
		str.replace(0, name.size(), name);

		char num[32];
		sprintf(num, "%d", factory->GetNrows());
		str.replace(22-strlen(num), strlen(num), num);

		const char *tag = factory->Tag();
		if(strlen(tag)){
			char tag_str[256];
			sprintf(tag_str, "\"%s\"", tag);
			str.replace(26, strlen(tag_str), tag_str);
		}
		
		cout<<str<<endl;
	}
	
	cout<<endl;

	return NOERROR;
}

//-------------
// Print
//-------------
derror_t DEventLoop::Print(const string data_name)
{
	/// Dump the data to stdout for the specified factory
	///
	/// Find the factory corresponding to data_name and send
	/// the return value of its toString() method to stdout.

	// Search for specified factory and return pointer to it's data container
	DFactory_base *factory = GetFactory(data_name);
	if(!factory){
		cerr<<" ERROR -- Factory not found for class \""<<data_name<<"\""<<endl;
		return NOERROR;
	}
	
	cout<<factory->toString();

	return NOERROR;
}

//-------------
// PrintCallStack
//-------------
void DEventLoop::PrintCallStack(void)
{
	// Create a list of the call strings while finding the longest one
	vector<string> routines;
	unsigned int max_length = 0;
	for(unsigned int i=0; i<call_stack.size(); i++){
		string routine = call_stack[i].factory_name;
		if(call_stack[i].tag){
			if(strlen(call_stack[i].tag)){
				routine = routine + ":" + call_stack[i].tag;
			}
		}
		if(routine.size()>max_length) max_length = routine.size();
		routines.push_back(routine);
	}

	stringstream sstr;
	sstr<<" Factory Call Stack"<<endl;
	sstr<<"============================"<<endl;
	for(unsigned int i=0; i<call_stack.size(); i++){
		string routine = routines[i];
		sstr<<" "<<routine<<string(max_length+2 - routine.size(),' ');
		if(call_stack[i].filename){
			sstr<<"--  "<<" line:"<<call_stack[i].line<<"  "<<call_stack[i].filename;
		}
		sstr<<endl;
	}
	sstr<<"----------------------------"<<endl;
	
	cout<<sstr.str();
}

//-------------
// Loop
//-------------
derror_t DEventLoop::Loop(void)
{
	/// Loop over events until Quit() method is called or we run
	/// out of events.
	
	do{
		// Let main thread know we're alive
		*heartbeat = 0.0;

		// Handle pauses and quits
		while(pause){
			*heartbeat = 0.0;	// Let main thread know we're alive
			usleep(500000);
			if(quit)break;
		}
		if(quit)break;
		
		// Read in a new event
		switch(OneEvent()){
			case NO_MORE_EVENT_SOURCES:
				// No more events. Time to quit
				quit = 1;
				break;
			case NOERROR:
				// Don't need to do anything here
				break;
			default:
				break;
		}
	
	}while(!quit);
	
	return NOERROR;
}

//-------------
// OneEvent
//-------------
derror_t DEventLoop::OneEvent(void)
{
	/// Read in and process one event. If eventno is
	/// less than 0, then grab the next event from
	/// the source. Otherwise, jump to the specified event

	// Clear evnt_called flag in all factories
	ClearFactories();

	// Try to read in an event
	derror_t err = app->NextEvent(event);
	
	// Here is a bit of voodoo. Calling setjmp() records the entire stack
	// so that a subsequent call to longjmp() will return us right here.
	// If we're returning here from a longjmp() call, then the return
	// value of setjmp() will be non-zero. The longjmp() call is made
	// from the SIGNINT interrupt handler so infinite loops can be
	// broken out of and the program will still gracefully exit.
	if(setjmp(SETJMP_ENV)){
		cerr<<endl<<"Uh-oh, seems the way-back machine was activated. Bailing this thread"<<endl<<endl;
		err = NO_MORE_EVENT_SOURCES;
	}
	
	switch(err){
		case NOERROR:
			break;
		case NO_MORE_EVENT_SOURCES:
			cout<<endl<<"No more event sources"<<endl;
			break;
		case EVENT_NOT_IN_MEMORY:
			cout<<endl<<"Event not in memory"<<endl;
			break;
		default:
			break;
	}
	if(err != NOERROR && err !=EVENT_NOT_IN_MEMORY)return err;
		
	// Call Event Processors
	int event_number = event.GetEventNumber();
	int run_number = event.GetRunNumber();
	vector<DEventProcessor*>::iterator p = processors.begin();
	for(; p!=processors.end(); p++){
		DEventProcessor *proc = *p;

		// Call brun routine if run number has changed or it's not been called
		proc->LockState();
		if(run_number!=proc->GetBRUN_RunNumber()){
			if(proc->brun_was_called() && !proc->erun_was_called()){
				proc->erun();
				proc->Set_erun_called();
			}
			proc->Clear_brun_called();
		}
		if(!proc->brun_was_called()){
			proc->brun(this, run_number);
			proc->Set_brun_called();
			proc->Clear_erun_called();
			proc->SetBRUN_RunNumber(run_number);
		}
		proc->UnlockState();

		// Initialize the factory call stack
		call_stack.clear();

		// Call the event routine
		try{
			proc->evnt(this, event_number);
		}catch(DException *exception){
			call_stack_t cs = {"DEventLoop", "OneEvent", __FILE__, __LINE__};
			call_stack.push_back(cs);
			PrintCallStack();
			throw exception;
		}
	}

	if(auto_free)event.FreeEvent();
			
	return NOERROR;
}

