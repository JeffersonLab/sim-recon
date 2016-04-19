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


#include <TTAB/DTranslationTable_factory.h>

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
	SWAP = true;
	LINK = true;
	PARSE_F250 = true;
	PARSE_F125 = true;
	PARSE_F1TDC = true;
	PARSE_CAEN1290TDC = true;
	PARSE_CONFIG = true;
	PARSE_BOR = true;
	PARSE_EPICS = true;
	PARSE_EVENTTAG = true;
	PARSE_TRIGGER = true;
	

	gPARMS->SetDefaultParameter("EVIO:VERBOSE", VERBOSE, "Set verbosity level for processing and debugging statements while parsing. 0=no debugging messages. 10=all messages");
	gPARMS->SetDefaultParameter("EVIO:NTHREADS", NTHREADS, "Set the number of worker threads to use for parsing the EVIO data");
	gPARMS->SetDefaultParameter("EVIO:MAX_PARSED_EVENTS", MAX_PARSED_EVENTS, "Set maximum number of events to allow in EVIO parsed events queue");
	gPARMS->SetDefaultParameter("EVIO:LOOP_FOREVER", LOOP_FOREVER, "If reading from EVIO file, keep re-opening file and re-reading events forever (only useful for debugging) If reading from ET, this is ignored.");
	gPARMS->SetDefaultParameter("EVIO:PRINT_STATS", PRINT_STATS, "Print some additional stats from event source when it's finished processing events");

	gPARMS->SetDefaultParameter("EVIO:SWAP", SWAP, "Allow swapping automatic swapping. Turning this off should only be used for debugging.");
	gPARMS->SetDefaultParameter("EVIO:LINK", LINK, "Link associated objects. Turning this off should only be used for debugging.");

	gPARMS->SetDefaultParameter("EVIO:PARSE_F250", PARSE_F250, "Set this to 0 to disable parsing of data from F250 ADC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_F125", PARSE_F125, "Set this to 0 to disable parsing of data from F125 ADC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_F1TDC", PARSE_F1TDC, "Set this to 0 to disable parsing of data from F1TDC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_CAEN1290TDC", PARSE_CAEN1290TDC, "Set this to 0 to disable parsing of data from CAEN 1290 TDC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_CONFIG", PARSE_CONFIG, "Set this to 0 to disable parsing of ROC configuration data in the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_BOR", PARSE_BOR, "Set this to 0 to disable parsing of BOR events from the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_EPICS", PARSE_EPICS, "Set this to 0 to disable parsing of EPICS events from the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_EVENTTAG", PARSE_EVENTTAG, "Set this to 0 to disable parsing of event tag data in the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_TRIGGER", PARSE_TRIGGER, "Set this to 0 to disable parsing of the built trigger bank from CODA (for benchmarking/debugging)");

	jobtype = DEVIOWorkerThread::JOB_FULL_PARSE;
	if(SWAP) jobtype |= DEVIOWorkerThread::JOB_SWAP;
	if(LINK) jobtype |= DEVIOWorkerThread::JOB_ASSOCIATE;

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
		w->PARSE_F250        = PARSE_F250;
		w->PARSE_F125        = PARSE_F125;
		w->PARSE_F1TDC       = PARSE_F1TDC;
		w->PARSE_CAEN1290TDC = PARSE_CAEN1290TDC;
		w->PARSE_CONFIG      = PARSE_CONFIG;
		w->PARSE_BOR         = PARSE_BOR;
		w->PARSE_EPICS       = PARSE_EPICS;
		w->PARSE_EVENTTAG    = PARSE_EVENTTAG;
		w->PARSE_TRIGGER     = PARSE_TRIGGER;
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
	uint64_t istreamorder = 0;
	while(true){
	
		if(japp->GetQuittingStatus()) break;

		// Get worker thread to handle this
		DEVIOWorkerThread *thr = NULL;
		while( !thr){
			for(uint32_t i=0; i<worker_threads.size(); i++){
				if(worker_threads[i]->in_use) continue;
				thr = worker_threads[i];
				break;
			}
			if(!thr) {
				NWAITS_FOR_THREAD++;
				this_thread::sleep_for(milliseconds(1));
			}
			if(DONE) break;
		}
		if(DONE) break;
		
		uint32_t* &buff     = thr->buff;
		uint32_t  &buff_len = thr->buff_len;

//		hdevio->readSparse(buff, buff_len, allow_swap);
		hdevio->readNoFileBuff(buff, buff_len, allow_swap);
//		hdevio->read(buff, buff_len, allow_swap);
		thr->pos = hdevio->last_event_pos;
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
		worker_threads[i]->done = true;
//		if(japp->GetQuittingStatus()) worker_threads[i]->done = true;
		while(worker_threads[i]->in_use){
			this_thread::sleep_for(milliseconds(10));
		}
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
	// CAUTION: This is called by the EventBufferThread in JANA which
	// effectively serializes everything done here. (Don't delete all
	// objects in pe which can be slow.)
	DParsedEvent *pe = (DParsedEvent*)event.GetRef();
	pe->in_use = false; // return pe to pool
	
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
	
	// When first called, the data objects will still be in the DParsedEvent
	// object and not yet copied into their respective factories. This will
	// do that copy of the pointers and set the "copied_to_factories" flag
	// in DParsedEvent indicating this has been done. It will then check the
	// name of the class being requested and return the appropriate value
	// depending on whether we supply that class or not.
	
	// We must have a factory to hold the data
	if(!factory)throw RESOURCE_UNAVAILABLE;

	// Get name of data class we're trying to extract and the factory tag
	string dataClassName = factory->GetDataClassName();
	string tag = factory->Tag();
	if(tag.length()!=0) return OBJECT_NOT_AVAILABLE; // can't provide tagged factory objects
	
	// Copy source objects for all classes into their respective factories
	// if needed.
	DParsedEvent *pe = (DParsedEvent*)event.GetRef();
	JEventLoop *loop = event.GetJEventLoop();

	if(!pe->copied_to_factories) pe->CopyToFactories(loop);

	// Apply any translation tables that exist, writing the translated objects
	// into their respective factories
	bool isSuppliedType = pe->IsParsedDataType(dataClassName);
	vector<const DTranslationTable*> translationTables;
	DTranslationTable_factory *ttfac = static_cast<DTranslationTable_factory*>(loop->GetFactory("DTranslationTable"));
	if(ttfac) ttfac->Get(translationTables);
	bool isNotTaggedFactory = strlen(factory->Tag()) == 0;
	for(unsigned int i=0; i<translationTables.size(); i++){
		translationTables[i]->ApplyTranslationTable(loop);
		if(!isSuppliedType){
			if(translationTables[i]->IsSuppliedType(dataClassName)){
				if(isNotTaggedFactory){ // Don't allow tagged factories from Translation table
					isSuppliedType = true;
				}
			}
		}
	}

	// Check if this is a class we provide and return appropriate value
	if( isSuppliedType ){
		// We do produce this type
		return NOERROR;
	}else{
		// We do not produce this type
		return OBJECT_NOT_AVAILABLE;
	}	
}

