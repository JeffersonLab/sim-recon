// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceFile methods
//

#include <iostream>
#include <iomanip>
using namespace std;

#include "DEventSourceFile.h"

//----------------
// Constructor
//----------------
DEventSourceFile::DEventSourceFile(int narg, char *argv[]):DEventSource(narg,argv)
{
	/// Constructor for DEventSourceFile object
	fin = NULL;
	hddm_s = NULL;
	max_event_buff = 10;
	current_event_index = -1;
}

//----------------
// Destructor
//----------------
DEventSourceFile::~DEventSourceFile()
{
	if(hddm_s)flush_s_HDDM(hddm_s, fin);
	if(fin)close_s_HDDM(fin);
}

//----------------
// Open
//----------------
derror_t DEventSourceFile::Open(char *source)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	Close();
	fin = open_s_HDDM(source);
	
	return fin==NULL ? ERROR_OPENING_EVENT_SOURCE:NOERROR;
}

//----------------
// Close
//----------------
derror_t DEventSourceFile::Close(void)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.
	if(fin)close_s_HDDM(fin);
	fin = NULL;
	
	return NOERROR;
}

//----------------
// GetEvent
//----------------
derror_t DEventSourceFile::GetEvent(void)
{
	/// Implementation of DEventSource virtual function
	/// Support only HDDM formatted files for now.

	// A previous call to GotoEvent() may mean we are looking
	// at an event in the middle of the buffer. The attribute
	// current_event_index holds the index to the event_buff
	// container of the current (or rather, previous after this
	// call) event. If it is not pointing to the last element,
	// Then call GotoEvent instead to use the values already in
	// memory. The value of current_event_index will be updated
	// in RecordEventNumber() below.
	if(current_event_index != (int)event_buff.size() - 1){
		const event_buffer_t &e = event_buff[current_event_index+1];
		return GotoEvent(e.event);
	}

	// Note: the data is freed ala flush_s_HDDM in RecordEventNumber() below
	hddm_s = read_s_HDDM(fin);
	if(!hddm_s)return NO_MORE_EVENTS_IN_SOURCE;
	
	return NOERROR;
}

//----------------
// GotoEvent
//----------------
derror_t DEventSourceFile::GotoEvent(int eventno)
{
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(; iter!=event_buff.end(); iter++){
		if((*iter).event == eventno){
			
			// Note: The value of current_event_index will be updated
			// in RecordEventNumber() below.
			hddm_s = (*iter).hddm_s;
			
			return NOERROR;
		}
	}

	return EVENT_NOT_IN_MEMORY;
}

//----------------
// SetEventBufferSize
//----------------
derror_t DEventSourceFile::SetEventBufferSize(int Nevents)
{
	max_event_buff = Nevents;

	// If we're reducing the number of events, delete oldest
	// records until we're within the new limit
	while(max_event_buff<(int)event_buff.size()){
		event_buffer_t &e = event_buff[0];
		flush_s_HDDM(e.hddm_s, 0);
		event_buff.erase(event_buff.begin(),event_buff.begin());
	}

	return NOERROR;
}

//----------------
// RecordEventNumber
//----------------
derror_t DEventSourceFile::RecordEventNumber(int eventno)
{
	/// This may seem a little convoluted, but here is the reasoning:
	/// The GotoEvent() method must be implemented specifically by
	/// each sub-class that chooses to implement it. This is because it will
	/// be specific to the source type. However, the source doesn't
	/// neccesarily know the event number. For example, if it is not
	/// kept in the file, but the DEventLoop object keeps track of it
	/// instead. Thus, DEventLoop must have some way of communicating
	/// the definitive event number back to the source.
	///
	/// One more complication is that by the time this method is
	/// called, the file pointer has already been advanced to the next
	/// event. Thus, the pointer has to be recorded in the GetEvent()
	/// method which corrsponds to the event number passed here.
	
	// Make sure we don't add an event twice
	vector<event_buffer_t>::iterator iter = event_buff.begin();
	for(int i=0; iter!=event_buff.end(); i++, iter++){
		if((*iter).event == eventno){
			current_event_index = i;
			return NOERROR;
		}
	}
	
	// Add this event to list
	event_buffer_t e;
	e.event = eventno;
	e.hddm_s = hddm_s;
	event_buff.push_back(e);

	// If we're over the limit, then delete an event
	if((int)event_buff.size() > max_event_buff){
		event_buffer_t &e = event_buff[0];
		flush_s_HDDM(e.hddm_s, 0);
		event_buff.erase(event_buff.begin());
	}
	
	current_event_index = (int)event_buff.size() - 1;

	return NOERROR;
}



