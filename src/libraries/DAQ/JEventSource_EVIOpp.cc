// $Id$
//
//    File: JEventSource_EVIOpp.cc
// Created: Tue Mar 29 08:14:42 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#include <unistd.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#include <forward_list>
#include <chrono>
#include <iostream>
#include <string>
#include <vector>
#include <stack>
#include <thread>
#include <mutex>
#include <condition_variable>

using namespace std;
using namespace std::chrono;



#include "JEventSource_EVIOpp.h"
using namespace jana;

//----------------
// Constructor
//----------------
JEventSource_EVIOpp::JEventSource_EVIOpp(const char* source_name):JEventSource(source_name)
{
	DONE = false;
	NEVENTS_PROCESSED = 0;
	NWAITS_FOR_THREAD = 0;
	NWAITS_FOR_PARSED_EVENT = 0;
	

	// Initialize dedicated JStreamLog used for debugging messages
	evioout.SetTag("--- EVIO ---: ");
	evioout.SetTimestampFlag();
	evioout.SetThreadstampFlag();
	
	// Define base set of status bits
	if(japp) DStatusBits::SetStatusBitDescriptions(japp);

	// Get configuration parameters
	VERBOSE = 0;
	NTHREADS = 2;
	MAX_PARSED_EVENTS = 128;
	LOOP_FOREVER = false;
	PRINT_STATS = true;

	gPARMS->SetDefaultParameter("EVIO:VERBOSE", VERBOSE, "Set verbosity level for processing and debugging statements while parsing. 0=no debugging messages. 10=all messages");
	gPARMS->SetDefaultParameter("EVIO:NTHREADS", NTHREADS, "Set the number of worker threads to use for parsing the EVIO data");
	gPARMS->SetDefaultParameter("EVIO:MAX_PARSED_EVENTS", MAX_PARSED_EVENTS, "Set maximum number of events to allow in EVIO parsed events queue");
	gPARMS->SetDefaultParameter("EVIO:LOOP_FOREVER", LOOP_FOREVER, "If reading from EVIO file, keep re-opening file and re-reading events forever (only useful for debugging) If reading from ET, this is ignored.");
	gPARMS->SetDefaultParameter("EVIO:PRINT_STATS", PRINT_STATS, "Print some additional stats from event source when it's finished processing events");

	// Try to open the file.
	if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as EVIO file..." <<endl;

	hdevio = new HDEVIO(this->source_name);
	
	if( hdevio->is_open ){
		if(VERBOSE>0) evioout << "Success opening event source \"" << this->source_name << "\"!" <<endl;
	}else{
		if(VERBOSE>0) evioout << "Unable to open event source \"" << this->source_name << "\"!" <<endl;
		cout << hdevio->err_mess.str() << endl;
		throw std::exception(); // throw exception indicating error
	}
	
	// Create dispatcher thread
	dispatcher_thread = new thread(&JEventSource_EVIOpp::Dispatcher, this);

	// Create worker threads
	for(int i=0; i<NTHREADS; i++){
		DEVIOWorkerThread *w = new DEVIOWorkerThread(parsed_events, MAX_PARSED_EVENTS, PARSED_EVENTS_MUTEX, PARSED_EVENTS_CV);
		worker_threads.push_back(w);
	}

	// Record start time
	tstart = high_resolution_clock::now();

}

//----------------
// Destructor
//----------------
JEventSource_EVIOpp::~JEventSource_EVIOpp()
{
	// Set DONE flag to tell dispatcher thread to quit
	// as well as anyone in a wait state
	DONE = true;
	PARSED_EVENTS_CV.notify_all();
	
	// Wait for dispatcher to complete
	if(dispatcher_thread){
		dispatcher_thread->join();
		delete dispatcher_thread;
	}

	// Wait for all worker threads to end and destroy them all
	for(uint32_t i=0; i<worker_threads.size(); i++){
		worker_threads[i]->Finish();
		delete worker_threads[i];
	}
	
	// Delete HDEVIO and print stats
	if(hdevio){
		if(VERBOSE>0) evioout << "Closing hdevio event source \"" << this->source_name << "\"" <<endl;
		
		if(PRINT_STATS){
			auto tdiff = duration_cast<duration<double>>(tend - tstart);
			double rate = (double)NEVENTS_PROCESSED/tdiff.count();

			cout << endl;
			cout << "   EVIO Processing rate = " << rate << " Hz" << endl;
			cout << "      NWAITS_FOR_THREAD = " << NWAITS_FOR_THREAD << endl;
			cout << "NWAITS_FOR_PARSED_EVENT = " << NWAITS_FOR_PARSED_EVENT << endl;
		}
		
		hdevio->PrintStats();
		delete hdevio;
	}

}

//----------------
// Dispatcher
//----------------
void JEventSource_EVIOpp::Dispatcher(void)
{
	/// This is run in a dedicated thread created by the constructor.
	/// It's job is to read in events and dispatch the processing of
	/// them to worker threads. When a worker thread is done, it adds
	/// the event to the parsed_events list and clears its own "in_use"
	/// flag thereby, marking itself as available for another job.
	/// The worker threads will stall if adding the event(s) it produced
	/// would make parsed_events contain more than MAX_PARSED_EVENTS.
	/// This creates backpressure here by having no worker threads
	/// available.
	
	bool allow_swap = false;
	uint32_t jobtype = DEVIOWorkerThread::JOB_SWAP | DEVIOWorkerThread::JOB_FULL_PARSE;
	uint64_t istreamorder = 0;
	while(true){
	
		if(japp->GetQuittingStatus()) break;

		// Get worker thread to handle this
		DEVIOWorkerThread *thr = NULL;
		while(!thr){
			for(uint32_t i=0; i<worker_threads.size(); i++){
				if(worker_threads[i]->in_use) continue;
				thr = worker_threads[i];
				break;
			}
			if(!thr) {
				NWAITS_FOR_THREAD++;
				this_thread::sleep_for(milliseconds(1));
			}
		}
		
		uint32_t* &buff     = thr->buff;
		uint32_t  &buff_len = thr->buff_len;

//		hdevio->readSparse(buff, buff_len, allow_swap);
		hdevio->readNoFileBuff(buff, buff_len, allow_swap);
//		hdevio->read(buff, buff_len, allow_swap);
		if(hdevio->err_code == HDEVIO::HDEVIO_USER_BUFFER_TOO_SMALL){
			delete[] buff;
			buff_len = hdevio->last_event_len;
			buff = new uint32_t[buff_len];
			continue;
		}else if(hdevio->err_code!=HDEVIO::HDEVIO_OK){
			if(LOOP_FOREVER && NEVENTS_PROCESSED>=1){
				if(hdevio){
					hdevio->rewind();
					continue;
				}
			}else{
				cout << hdevio->err_mess.str() << endl;
			}
			break;
		}
		
		// Wake up worker thread to handle event
		thr->in_use = true;
		thr->jobtype = (DEVIOWorkerThread::JOBTYPE)jobtype;
		thr->istreamorder = istreamorder++;

		thr->cv.notify_all();
	}
	
	// If we get here then there are no more events in the source.
	// Wait for all worker threads to become available so we know 
	// the system is drained of events. Then set the DONE flag so
	// GetEvent will properly return NO_MORE_EVENTS_IN_SOURCE.
	for(uint32_t i=0; i<worker_threads.size(); i++){
		if(japp->GetQuittingStatus()) worker_threads[i]->done = true;
		while(worker_threads[i]->in_use) this_thread::sleep_for(milliseconds(10));
	}
	
	tend = std::chrono::high_resolution_clock::now();
	DONE = true;
}

//----------------
// GetEvent
//----------------
jerror_t JEventSource_EVIOpp::GetEvent(JEvent &event)
{
	// Get next event from list, waiting if necessary
	unique_lock<std::mutex> lck(PARSED_EVENTS_MUTEX);
	while(parsed_events.empty()){
		if(DONE) return NO_MORE_EVENTS_IN_SOURCE;
		NWAITS_FOR_PARSED_EVENT++;
		PARSED_EVENTS_CV.wait_for(lck,std::chrono::milliseconds(1));
	}

	DParsedEvent *pe = parsed_events.front();
	parsed_events.pop_front();
	
	// Release mutex and notify workers they can use it again
	lck.unlock();
	PARSED_EVENTS_CV.notify_all();

	// Copy info for this parsed event into the JEvent
	event.SetJEventSource(this);
	event.SetEventNumber(pe->event_number);
	event.SetRunNumber(pe->run_number);
	event.SetRef(pe);
	
	return NOERROR;
}

//----------------
// FreeEvent
//----------------
void JEventSource_EVIOpp::FreeEvent(JEvent &event)
{
	DParsedEvent *pe = (DParsedEvent*)event.GetRef();
	pe->in_use = false;
//	delete pe;
	
	NEVENTS_PROCESSED++;
}

//----------------
// GetObjects
//----------------
jerror_t JEventSource_EVIOpp::GetObjects(JEvent &event, JFactory_base *factory)
{
	// This callback is called to extract objects of a specific type from
	// an event and store them in the factory pointed to by JFactory_base.
	// The data type desired can be obtained via factory->GetDataClassName()
	// and the tag via factory->Tag().
	//
	// If the object is not one of a type this source can provide, then
	// it should return OBJECT_NOT_AVAILABLE. Otherwise, it should return
	// NOERROR;
	
	// We must have a factory to hold the data
	if(!factory)throw RESOURCE_UNAVAILABLE;

	// Get name of data class we're trying to extract and the factory tag
	string dataClassName = factory->GetDataClassName();
	string tag = factory->Tag();
	
	DParsedEvent *pe = (DParsedEvent*)event.GetRef();
//_DBG_ << pe->vDEPICSvalue.size() << endl;
	
	// Example for providing objects of type XXX
	//
	// // Get pointer to object of type MyEvent (this is optional) 
	// MyEvent *myevent = (MyEvent*)event.GetRef();
	//
	// if(dataClassName == "XXX"){
	//
	//	 // Make objects of type XXX storing them in a vector
	//   vector<XXX*> xxx_objs;
	//	 ...
	//
	//	 JFactory<XXX> *fac = dynamic_cast<JFactory<XXX>*>(factory);
	//	 if(fac)fac->CopyTo(xxx_objs);
	//	
	//	 return NOERROR;
	// }
	
	// For all other object types, just return OBJECT_NOT_AVAILABLE to indicate
	// we can't provide this type of object
	return OBJECT_NOT_AVAILABLE;
}

