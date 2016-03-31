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

// This implements an alternative to thread_local. It is
// needed because Apple stupidly decided not to include
// thread_local support in their default compiler. This
// allocates a fixed size array of vectors to hold a pool
// of pointers to DParsedEvent objects. Each DEVIOWorkerThread
// object created (should be one per worker thread) will
// increment the counter and claim one of the elements.
// A check is made in the DEVIOWorkerThread constructor to
// make sure we don't override the array size.
#define TLS_PARSED_EVENT_MAX 512
atomic<int> TLS_PARSED_EVENT_MAX_IDX(0);
vector<DParsedEvent*> TLS_PARSED_EVENT[TLS_PARSED_EVENT_MAX];


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
	,parsed_event_pool(TLS_PARSED_EVENT[TLS_PARSED_EVENT_MAX_IDX.fetch_add(1)])
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

	// Check that we didn't overrun the allocated space	
	if(TLS_PARSED_EVENT_MAX_IDX >= TLS_PARSED_EVENT_MAX){
		_DBG_ << "-- ERROR: more than TLS_PARSED_EVENT_MAX (=" << TLS_PARSED_EVENT_MAX << ") DEVIOWorkerThread objects instantiated!" << endl;
		_DBG_ << "--        To fix this either run with fewer DEVIOWorkerThreads or modify the" << endl;
		_DBG_ << "          #define for TLS_PARSED_EVENT_MAX at the top of DEVIOWorkerThread.cc" << endl;
		exit(-1);
	}

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
	for(auto pe : parsed_event_pool) delete pe;
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

		if( jobtype & JOB_FULL_PARSE ) MakeEvents();
		
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
// MakeEvents
//---------------------------------
void DEVIOWorkerThread::MakeEvents(void)
{
	
	/// Make DParsedEvent objects from data currently in buff.
	
	uint32_t *iptr = buff;
	
	uint32_t M = 1;
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

	// Try and get M DParsedEvent objects from this thread's pool.
	current_parsed_events.clear();
	for(auto pe : parsed_event_pool){
		if(pe->in_use) continue;
		current_parsed_events.push_back(pe);
		if( current_parsed_events.size() >= M ) break;
	}
	
	// Create new DParsedEvent objects if needed
	while( current_parsed_events.size() < M ){
		DParsedEvent *pe = new DParsedEvent(0);
		current_parsed_events.push_back(pe);
		parsed_event_pool.push_back(pe);
	}
	
	// Set indexes for the parsed event objects
	// and flag them as being in use.
	for(auto pe : current_parsed_events){
	
		pe->istreamorder = istreamorder;
		pe->event_number = event_num++;
		pe->sync_flag = false;
		pe->Clear();
		pe->in_use = true;
	}
	
	// Parse data in buffer to create data objects
	ParseBank();
	
	
	//-----------------------------------------------
	
	// Lock mutex so other threads can't modify parsed_events
	unique_lock<mutex> lck(PARSED_EVENTS_MUTEX);
	
	// Make sure we don't exceed the maximum number of simultaneous
	// parsed events. If the done flag is set, go ahead and add
	// this regardless
	while( ((current_parsed_events.size()+parsed_events.size())>=MAX_PARSED_EVENTS) && !done ){
		PARSED_EVENTS_CV.wait_for(lck, std::chrono::milliseconds(1));
	}
	
	// Loop over all elements of parsed_events and insert
	// these based on istreamorder so that the front element
	// is the most recent.
	bool inserted = false;
	for(auto it = parsed_events.begin(); it!=parsed_events.end(); it++){
		if( istreamorder < (*it)->istreamorder ){
			parsed_events.insert(it, current_parsed_events.begin(), current_parsed_events.end());
			inserted = true;
			break;
		}
	}
	
	// In case this should go at end of list
	if(!inserted) parsed_events.insert(parsed_events.end(), current_parsed_events.begin(), current_parsed_events.end());

	lck.unlock();
	PARSED_EVENTS_CV.notify_all();
}

//---------------------------------
// ParseBank
//---------------------------------
void DEVIOWorkerThread::ParseBank(void)
{
//this_thread::sleep_for(milliseconds(1));

	uint32_t *iptr = buff;
	uint32_t *iend = &buff[buff[0]+1];

	while(iptr < iend){
		uint32_t event_len  = iptr[0];
		uint32_t event_head = iptr[1];
		uint32_t tag = (event_head >> 16) & 0xFFFF;

//_DBG_ << "tag=" << hex << tag << dec << endl;

		switch(tag){
			case 0x0056:    ParseEventTagBank(iptr, iend);    break;
			case 0x0060:       ParseEPICSbank(iptr, iend);    break;
			case 0x0070:         ParseBORbank(iptr, iend);    break;
			case 0xEE02:    ParseTSscalerBank(iptr, iend);    break;
			case 0xEE05:  Parsef250scalerBank(iptr, iend);    break;
			case 0xFF58:
			case 0xFF78:
							current_parsed_events.back()->sync_flag = true;
			case 0xFF50:     
			case 0xFF70:     
							 ParsePhysicsBank(iptr, iend);    break;
			default:
				_DBG_ << "Unknown outer EVIO bank tag: " << hex << tag << dec << endl;
				iptr = &iptr[event_len+1];
				if(event_len<1) iptr = iend;		
		}
	}
}

//---------------------------------
// ParseEventTagBank
//---------------------------------
void DEVIOWorkerThread::ParseEventTagBank(uint32_t* &iptr, uint32_t *iend)
{
	iptr = &iptr[(*iptr) + 1];
}

//---------------------------------
// ParseEPICSbank
//---------------------------------
void DEVIOWorkerThread::ParseEPICSbank(uint32_t* &iptr, uint32_t *iend)
{

	time_t timestamp=0;
	
	// Outer bank
	uint32_t *istart = iptr;
	uint32_t epics_bank_len = *iptr++;	
	if(epics_bank_len < 1){
		_DBG_ << "bank_len<1 in EPICS event!" << endl;
		iptr = iend;
		return;
	}
	
	uint32_t *iend_epics = &iptr[epics_bank_len];
	if( iend_epics < iend ) iend = iend_epics;
	
	// Advance to first daughter bank
	iptr++;
	
	// Get pointer to first DParsedEvent
	DParsedEvent *pe = current_parsed_events.front();
	
	// Loop over daughter banks
	while( iptr < iend_epics ){

		uint32_t bank_len =  (*iptr)&0xFFFF;
		uint32_t tag      = ((*iptr)>>24)&0xFF;
		iptr++;
	
		if(tag == 0x61){
			// timestamp bank
			timestamp = *iptr;
		}else if(tag == 0x62){
			// EPICS data value
			string nameval = (const char*)iptr;
//_DBG_<<"Making DEPICSvalue: " << nameval << endl;
			DEPICSvalue *epicsval = new DEPICSvalue(timestamp, nameval);
			pe->vDEPICSvalue.push_back(epicsval);
		}else{
			// Unknown tag. Bail
			_DBG_ << "Unknown tag 0x" << hex << tag << dec << " in EPICS event!" <<endl;
			DumpBinary(istart, iend_epics, 32, &iptr[-1]);
		}
		
		iptr = &iptr[bank_len];
	}	

	iptr = iend_epics;
}

//---------------------------------
// ParseBORbank
//---------------------------------
void DEVIOWorkerThread::ParseBORbank(uint32_t* &iptr, uint32_t *iend)
{
	iptr = &iptr[(*iptr) + 1];
}

//---------------------------------
// ParseTSscalerBank
//---------------------------------
void DEVIOWorkerThread::ParseTSscalerBank(uint32_t* &iptr, uint32_t *iend)
{
	iptr = &iptr[(*iptr) + 1];
}

//---------------------------------
// Parsef250scalerBank
//---------------------------------
void DEVIOWorkerThread::Parsef250scalerBank(uint32_t* &iptr, uint32_t *iend)
{
	iptr = &iptr[(*iptr) + 1];
}

//---------------------------------
// ParsePhysicsBank
//---------------------------------
void DEVIOWorkerThread::ParsePhysicsBank(uint32_t* &iptr, uint32_t *iend)
{
	iptr = &iptr[(*iptr) + 1];
}


//----------------
// DumpBinary
//----------------
void DEVIOWorkerThread::DumpBinary(const uint32_t *iptr, const uint32_t *iend, uint32_t MaxWords, const uint32_t *imark)
{
    /// This is used for debugging. It will print to the screen the words
    /// starting at the address given by iptr and ending just before iend
    /// or for MaxWords words, whichever comes first. If iend is NULL,
    /// then MaxWords will be printed. If MaxWords is zero then it is ignored
    /// and only iend is checked. If both iend==NULL and MaxWords==0, then
    /// only the word at iptr is printed.

    cout << "Dumping binary: istart=" << hex << iptr << " iend=" << iend << " MaxWords=" << dec << MaxWords << endl;

    if(iend==NULL && MaxWords==0) MaxWords=1;
    if(MaxWords==0) MaxWords = (uint32_t)0xffffffff;

    uint32_t Nwords=0;
    while(iptr!=iend && Nwords<MaxWords){

        // line1 is hex and line2 is decimal
        stringstream line1, line2;

        // print words in columns 8 words wide. First part is
        // reserved for word number
        uint32_t Ncols = 8;
        line1 << setw(5) << Nwords;
        line2 << string(5, ' ');

        // Loop over columns
        for(uint32_t i=0; i<Ncols; i++, iptr++, Nwords++){

            if(iptr == iend) break;
            if(Nwords>=MaxWords) break;

            stringstream iptr_hex;
            iptr_hex << hex << "0x" << *iptr;

            string mark = (iptr==imark ? "*":" ");

            line1 << setw(12) << iptr_hex.str() << mark;
            line2 << setw(12) << *iptr << mark;
        }

        cout << line1.str() << endl;
        cout << line2.str() << endl;
        cout << endl;
    }
}

