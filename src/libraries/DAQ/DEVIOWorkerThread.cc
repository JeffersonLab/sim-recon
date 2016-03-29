// $Id$
//
//    File: DEVIOWorkerThread.cc
// Created: Mon Mar 28 07:40:07 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#include "DEVIOWorkerThread.h"

#include <swap_bank.h>

using namespace std;
using namespace std::chrono;


//uint32_t MAX_PARSED_EVENTS = 128;
//mutex PARSED_EVENTS_MUTEX;
//condition_variable PARSED_EVENTS_CV;
//list<DParsedEvent*> parsed_events;


//---------------------------------
// DEVIOWorkerThread    (Constructor)
//---------------------------------
DEVIOWorkerThread::DEVIOWorkerThread(
	 list<DParsedEvent*>  &parsed_events
	 ,uint32_t            &MAX_PARSED_EVENTS
	 ,mutex               &PARSED_EVENTS_MUTEX
	 ,condition_variable  &PARSED_EVENTS_CV
	 ):
	 parsed_events(parsed_events)
	,MAX_PARSED_EVENTS(MAX_PARSED_EVENTS)
	,PARSED_EVENTS_MUTEX(PARSED_EVENTS_MUTEX)
	,PARSED_EVENTS_CV(PARSED_EVENTS_CV)
	,done(false)
	,thd(&DEVIOWorkerThread::Run,this)
{
	// n.b. in principal, the worker thread is started when the
	// above constructor is hit and so may already be in Run()
	// before executing anything below. The "done" variable is
	// therefore initialized first to guarantee that if that
	// happens, it gets to the cv.wait() call where it will wait
	// for someone to notify it. That won't happen before this
	// constructor completes so we do the remaining initializations
	// below.

	in_use  = false;
	jobtype = JOB_NONE;

	buff_len = 100; // this will grow as needed
	buff = new uint32_t[buff_len];
}

//---------------------------------
// ~DEVIOWorkerThread    (Destructor)
//---------------------------------
DEVIOWorkerThread::~DEVIOWorkerThread()
{
	if(buff) delete[] buff;
}

//---------------------------------
// Run
//---------------------------------
void DEVIOWorkerThread::Run(void)
{
	unique_lock<std::mutex> lck(mtx);

	// Loop waiting for jobs or until told to quit	
	while(!done){

		cv.wait_for(lck, std::chrono::milliseconds(1));

		if( jobtype & JOB_SWAP       ) swap_bank(buff, buff, swap32(buff[1]) );

		if( jobtype & JOB_FULL_PARSE ) ParseBank();
		
		if( jobtype & JOB_QUIT       ) break;
		
		// Reset and mark us as available for use
		jobtype = JOB_NONE;
		in_use  = false;
	}
}

//---------------------------------
// Finish
//---------------------------------
void DEVIOWorkerThread::Finish(bool wait_to_complete)
{
	/// Set the done flag so that the worker thread
	/// will exit once it is done processing its current
	/// job. The thread is notified to wake up in case
	/// it is currently idle. If the wait_to_complete
	/// flag is set (default), then the worker thread is
	/// joined to guarantee the current job's processing
	/// is completed before returning.
	done = true;
	cv.notify_all();
	if(wait_to_complete) {
		thd.join();
	} else {
		thd.detach();
	}
}

//---------------------------------
// ParseBank
//---------------------------------
void DEVIOWorkerThread::ParseBank(void)
{
	list<DParsedEvent*> my_parsed_events;
	
	//------- Parse event, making hit objects -------

	uint32_t *iptr = buff;
	uint32_t *iend = &buff[buff[0]+1];
	
	int M = 1;
	uint64_t event_num = 0;

	iptr++;
	uint32_t mask = 0xFF001000;
	if( ((*iptr)&mask) == mask ){
		// Physics event
		M = *(iptr)&0xFF;
		uint64_t eventnum_lo = iptr[4];
		uint64_t eventnum_hi = iptr[5];
		event_num = (eventnum_hi<<32) + (eventnum_lo);
	}

	for(int i=0; i<M; i++){
		DParsedEvent *pe = new DParsedEvent(istreamorder);
		
		pe->event_number = event_num++;
		
		my_parsed_events.push_back(pe);
		
this_thread::sleep_for(milliseconds(1));
	}
	//-----------------------------------------------
	
	// Lock mutex so other threads can't modify parsed_events
	unique_lock<mutex> lck(PARSED_EVENTS_MUTEX);
	
	// Make sure we don't exceed the maximum number of simultaneous
	// parsed events. If the done flag is set, go ahead and add
	// this regardless
	while( ((my_parsed_events.size()+parsed_events.size())>=MAX_PARSED_EVENTS) && !done ){
		PARSED_EVENTS_CV.wait_for(lck, std::chrono::milliseconds(1));
	}
	
	// Loop over all elements of parsed_events and insert
	// these based on istreamorder so that the front element
	// is the most recent.
	bool inserted = false;
	for(auto it = parsed_events.begin(); it!=parsed_events.end(); it++){
		if( istreamorder < (*it)->istreamorder ){
			parsed_events.insert(it, my_parsed_events.begin(), my_parsed_events.end());
			inserted = true;
			break;
		}
	}
	
	// In case this should go at end of list
	if(!inserted) parsed_events.insert(parsed_events.end(), my_parsed_events.begin(), my_parsed_events.end());
	
	PARSED_EVENTS_CV.notify_all();
}

