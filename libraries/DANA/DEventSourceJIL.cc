// $Id$
//
//    File: DEventSourceJIL.cc
// Created: Sun Dec  4 16:51:23 EST 2005
// Creator: davidl (on Linux localhost.localdomain 2.6.9-1.667 i686)
//

#ifdef JILIO

#include "hd_serializers.h"

#include "DEventSourceJIL.h"
#include "DFactory_base.h"
#include "DEventLoop.h"
#include "DEvent.h"

//---------------------------------
// DEventSourceJIL    (Constructor)
//---------------------------------
DEventSourceJIL::DEventSourceJIL(const char* source_name):DEventSource(source_name)
{
	jilstream = new JILStream(source_name, "r");
	event_id_counter=0;
}

//---------------------------------
// ~DEventSourceJIL    (Destructor)
//---------------------------------
DEventSourceJIL::~DEventSourceJIL()
{
	if(event_buff.size() > 0){
		cout<<__FILE__<<":"<<__LINE__<<" Hmmm... looks like "<<event_buff.size();
		cout<<" object buffers weren't deleted. Is this a bug?"<<endl;
	}

	// Delete any JILObjectRecords still in memory
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		jilstream->DeleteObjectRecords((*iter).objects);
	}

	delete jilstream;
}

//---------------------------------
// GetEvent
//---------------------------------
derror_t DEventSourceJIL::GetEvent(DEvent &event)
{
	// Thread is already locked before calling us and released
	// after we exit

	if(!jilstream)return EVENT_SOURCE_NOT_OPEN;
	
	// Get the next event from the JILStream and "adopt" the object records
	if(!jilstream->GetNamed("Event"))return NO_MORE_EVENTS_IN_SOURCE;
	event_buffer_t e;
	//jilstream->PrintObjectHisto();
	jilstream->AdoptObjectRecords(e.objects);
	e.event_id = event_id_counter++;
	
	// --- HERE WE NEED TO EXTRACT THE EVENT AND RUN NUMBERS SOMEHOW! ---
	e.event_number = ++Nevents_read;
	int run_number = 1;

	// Copy the reference info into the DEvent object
	event.SetDEventSource(this);
	event.SetEventNumber(e.event_number);
	event.SetRunNumber(run_number);
	event.SetRef((void*)e.event_id);
	
	event_buff.push_back(e);

	return NOERROR;
}

//---------------------------------
// GetObjects
//---------------------------------
derror_t DEventSourceJIL::GetObjects(const char *name, vector<void*> &v, const char* tag, void* ref, DFactory_base *factory)
{
	if(!factory)throw RESOURCE_UNAVAILABLE;

	unsigned long event_id = (unsigned long)ref;
	LockRead();
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		event_buffer_t &e = *iter;
		if(e.event_id == event_id){
			factory->StreamFromInput(jilstream, e.objects, v);
			UnlockRead();
			return v.size()>0 ? NOERROR:OBJECT_NOT_AVAILABLE;
		}
	}
	UnlockRead();

	// If we get to here, then the reference wasn't in our list
	throw RESOURCE_UNAVAILABLE;

	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void DEventSourceJIL::FreeEvent(void *ref)
{
	unsigned long event_id = (unsigned long)ref;
	LockRead();
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		if((*iter).event_id == event_id){
			jilstream->DeleteObjectRecords((*iter).objects);
			event_buff.erase(iter);
			break;
		}
	}
	UnlockRead();
}

#else	// JILIO
bool DEventSourceJIL_is_not_used = true;
#endif // JILIO
