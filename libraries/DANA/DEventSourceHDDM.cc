// $Id$
//
// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceHDDM methods
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DEventSourceHDDM.h"
#include "DFactory_base.h"
#include "DEventLoop.h"
#include "DEvent.h"

//----------------
// Constructor
//----------------
DEventSourceHDDM::DEventSourceHDDM(const char* source_name):DEventSource(source_name)
{
	/// Constructor for DEventSourceHDDM object
	fin = open_s_HDDM((char*)source_name);
	hddm_s = NULL;
	
	if(fin)source_is_open = 1;
	flush_on_free = true;
}

//----------------
// Destructor
//----------------
DEventSourceHDDM::~DEventSourceHDDM()
{
	//LockRead();

	if(event_buff.size() > 0){
		cout<<__FILE__<<":"<<__LINE__<<" Hmmm... looks like "<<event_buff.size();
		cout<<" s_HDDM_t structures weren't deleted. Is this a bug?"<<endl;
	}

	// Delete any hddm_s structures still in memory
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		s_HDDM_t *my_hddm_s = (*iter).hddm_s;
		if(my_hddm_s)flush_s_HDDM(my_hddm_s, fin);
	}

	// Set pointers to NULL and close file
	hddm_s = NULL;
	fin = NULL;
	if(fin)close_s_HDDM(fin);

	//UnlockRead();
}

//----------------
// GetEvent
//----------------
derror_t DEventSourceHDDM::GetEvent(DEvent &event)
{
	/// Implementation of DEventSource virtual function
		
	if(!fin){
		return EVENT_SOURCE_NOT_OPEN;
	}
	hddm_s = read_s_HDDM(fin);
	
	// --- HERE WE NEED TO EXTRACT THE EVENT AND RUN NUMBERS SOMEHOW! ---
	int event_number = ++Nevents_read;
	int run_number = 1;
	
	// Copy the reference info into the DEvent object
	event.SetDEventSource(this);
	event.SetEventNumber(event_number);
	event.SetRunNumber(run_number);
	event.SetRef(hddm_s);
	
	// Add this event to our private buffer
	if(hddm_s){
		event_buffer_t e;
		e.hddm_s = hddm_s;
		e.event_number = event_number;
		LockRead();
		event_buff.push_back(e);
		UnlockRead();
	}

	// each open HDDM file takes up about 1M of memeory so it's
	// worthwhile to close it as soon as we can
	if(!hddm_s){
		if(fin)close_s_HDDM(fin);
		fin = NULL;
	}

	return hddm_s==NULL ? NO_MORE_EVENTS_IN_SOURCE:NOERROR;
}

//----------------
// FreeEvent
//----------------
void DEventSourceHDDM::FreeEvent(void *ref)
{
	s_HDDM_t *my_hddm_s = (s_HDDM_t*)ref;
	LockRead();
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		if((*iter).hddm_s == my_hddm_s){
			if(flush_on_free)flush_s_HDDM((*iter).hddm_s, 0);
			event_buff.erase(iter);
			break;
		}
	}
	UnlockRead();
}

//----------------
// Goto
//----------------
derror_t DEventSourceHDDM::Goto(int eventno)
{
	/// (This is only half written so it won't work right)
	/// Look through the event buffer for the event with the
	/// specified number and set the hddm_s attribute to 
	/// point to it.

	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		if(iter->event_number == eventno){
			
			hddm_s = (*iter).hddm_s;
			if(!hddm_s){
				cerr<<__FILE__<<":"<<__LINE__<<" hddm_s is NULL!"<<endl;
			}
			
			return NOERROR;
		}
	}

	return EVENT_NOT_IN_MEMORY;
}

//----------------
// GetObjects
//----------------
derror_t DEventSourceHDDM::GetObjects(const char *name, vector<void*> &v, const char* tag, void* ref, DFactory_base *factory)
{
	// The extraction code actually resides in the factories.
	// Factories which don't define the virtual method Extract_HDDM()
	// will fall back to the definition in the DFactory_base class
	// which just returns OBJECT_NOT_AVAILABLE.

	if(!factory)throw RESOURCE_UNAVAILABLE;
	
	s_HDDM_t *my_hddm_s = (s_HDDM_t*)ref;
	if(!my_hddm_s)throw RESOURCE_UNAVAILABLE;

	// Verify that the hddm_s given to us is still valid
	LockRead();
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		if((*iter).hddm_s == my_hddm_s)break;
	}
	if(iter == event_buff.end()){
		cerr<<__FILE__<<":"<<__LINE__<<" hddm_s not in buffer!!"<<endl;
		UnlockRead();
		throw RESOURCE_UNAVAILABLE;
	}
	
	derror_t err = factory->Extract_HDDM(my_hddm_s, v);
	
	UnlockRead();
	
	return err;
}


