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
#include <cmath>

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
	USER_RUN_NUMBER = 0;
	ET_STATION_NEVENTS = 10;
	ET_STATION_CREATE_BLOCKING = false;
	PRINT_STATS = true;
	SWAP = true;
	LINK = true;
	PARSE = true;
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
	gPARMS->SetDefaultParameter("EVIO:RUN_NUMBER", USER_RUN_NUMBER, "User-supplied run number. Override run number from other sources with this.(will be ignored if set to zero)");
	gPARMS->SetDefaultParameter("EVIO:ET_STATION_NEVENTS", ET_STATION_NEVENTS, "Number of events to use if we have to create the ET station. Ignored if station already exists.");
	gPARMS->SetDefaultParameter("EVIO:ET_STATION_CREATE_BLOCKING", ET_STATION_CREATE_BLOCKING, "Set this to 0 to create station in non-blocking mode (default is to create it in blocking mode). Ignored if station already exists.");
	gPARMS->SetDefaultParameter("EVIO:PRINT_STATS", PRINT_STATS, "Print some additional stats from event source when it's finished processing events");

	gPARMS->SetDefaultParameter("EVIO:SWAP", SWAP, "Allow swapping automatic swapping. Turning this off should only be used for debugging.");
	gPARMS->SetDefaultParameter("EVIO:LINK", LINK, "Link associated objects. Turning this off should only be used for debugging.");

	gPARMS->SetDefaultParameter("EVIO:PARSE", PARSE, "Set this to 0 to disable parsing of event buffers and generation of any objects (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_F250", PARSE_F250, "Set this to 0 to disable parsing of data from F250 ADC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_F125", PARSE_F125, "Set this to 0 to disable parsing of data from F125 ADC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_F1TDC", PARSE_F1TDC, "Set this to 0 to disable parsing of data from F1TDC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_CAEN1290TDC", PARSE_CAEN1290TDC, "Set this to 0 to disable parsing of data from CAEN 1290 TDC modules (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_CONFIG", PARSE_CONFIG, "Set this to 0 to disable parsing of ROC configuration data in the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_BOR", PARSE_BOR, "Set this to 0 to disable parsing of BOR events from the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_EPICS", PARSE_EPICS, "Set this to 0 to disable parsing of EPICS events from the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_EVENTTAG", PARSE_EVENTTAG, "Set this to 0 to disable parsing of event tag data in the data stream (for benchmarking/debugging)");
	gPARMS->SetDefaultParameter("EVIO:PARSE_TRIGGER", PARSE_TRIGGER, "Set this to 0 to disable parsing of the built trigger bank from CODA (for benchmarking/debugging)");

	jobtype = DEVIOWorkerThread::JOB_NONE;
	if(PARSE) jobtype |= DEVIOWorkerThread::JOB_FULL_PARSE;
	if(LINK ) jobtype |= DEVIOWorkerThread::JOB_ASSOCIATE;
	
	source_type          = kNoSource;
	hdevio               = NULL;
	hdet                 = NULL;
	et_quit_next_timeout = false;

	// Open either ET system or file
	if(this->source_name.find("ET:") == 0){

		// Try to open ET system
		if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as ET (network) source..." <<endl;

		hdet = new HDET(this->source_name, ET_STATION_NEVENTS, ET_STATION_CREATE_BLOCKING);
		if( ! hdet->is_connected ){
			cerr << hdet->err_mess.str() << endl;
			throw JException("Failed to open ET system: " + this->source_name, __FILE__, __LINE__);
		}
		source_type = kETSource;


	}else{

		// Try to open the file.
		if(VERBOSE>0) evioout << "Attempting to open \""<<this->source_name<<"\" as EVIO file..." <<endl;

		hdevio = new HDEVIO(this->source_name);
		if( ! hdevio->is_open ){
			cerr << hdevio->err_mess.str() << endl;
			throw JException("Failed to open EVIO file: " + this->source_name, __FILE__, __LINE__); // throw exception indicating error
		}
		source_type = kFileSource;
	}

	if(VERBOSE>0) evioout << "Success opening event source \"" << this->source_name << "\"!" <<endl;
	
	// Create dispatcher thread
	dispatcher_thread = new thread(&JEventSource_EVIOpp::Dispatcher, this);

	// Create worker threads
	for(uint32_t i=0; i<NTHREADS; i++){
		DEVIOWorkerThread *w = new DEVIOWorkerThread(this, parsed_events, MAX_PARSED_EVENTS, PARSED_EVENTS_MUTEX, PARSED_EVENTS_CV);
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
	
	if(VERBOSE>0) evioout << "Closing hdevio event source \"" << this->source_name << "\"" <<endl;
	if(PRINT_STATS){
		auto tdiff = duration_cast<duration<double>>(tend - tstart);
		double rate = (double)NEVENTS_PROCESSED/tdiff.count();

		cout << endl;
		cout << "   EVIO Processing rate = " << rate << " Hz" << endl;
		cout << "      NWAITS_FOR_THREAD = " << NWAITS_FOR_THREAD << endl;
		cout << "NWAITS_FOR_PARSED_EVENT = " << NWAITS_FOR_PARSED_EVENT << endl;
	}
	
	// Delete all BOR objects
	for(auto p : borptrs_list) delete p;

	// Delete HDEVIO and print stats
	if(hdevio){
		hdevio->PrintStats();
		delete hdevio;
	}

	// Delete HDET and print stats
	if(hdet){
		hdet->PrintStats();
		delete hdet;
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
	
	bool allow_swap = false; // Defer swapping to DEVIOWorkerThread
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
		
		// Reduce average memory usage.
		if(++thr->Nrecycled%thr->MAX_RECYCLED == 0) thr->Prune();
		
		uint32_t* &buff     = thr->buff;
		uint32_t  &buff_len = thr->buff_len;
		
		bool swap_needed = false;

		if(source_type==kFileSource){
			// ---- Read From File ----
//			hdevio->read(buff, buff_len, allow_swap);
//			hdevio->readSparse(buff, buff_len, allow_swap);
			hdevio->readNoFileBuff(buff, buff_len, allow_swap);
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
			}else{
				// HDEVIO_OK
				swap_needed = hdevio->swap_needed;
			}
		}else{
			// ---- Read From ET ----
			hdet->read(buff, buff_len, allow_swap);
			thr->pos = 0;
			if(hdet->err_code == HDET::HDET_TIMEOUT){
				if(et_quit_next_timeout) break;
				this_thread::sleep_for(milliseconds(1));
				continue;
			}else if(hdet->err_code != HDET::HDET_OK){
				cout << hdet->err_mess.str() << endl;
				break;
			}else{
				// HDET_OK
				swap_needed = hdet->swap_needed;
			}
		}

		uint32_t myjobtype = jobtype;
		if(swap_needed && SWAP) myjobtype |= DEVIOWorkerThread::JOB_SWAP;
		
		// Wake up worker thread to handle event
		thr->in_use = true;
		thr->jobtype = (DEVIOWorkerThread::JOBTYPE)myjobtype;
		thr->istreamorder = istreamorder++;

		thr->cv.notify_all();
	}
	
	// If we get here then there are no more events in the source.
	// Wait for all worker threads to become available so we know 
	// the system is drained of events. Then set the DONE flag so
	// GetEvent will properly return NO_MORE_EVENTS_IN_SOURCE.
	for(uint32_t i=0; i<worker_threads.size(); i++){
		worker_threads[i]->done = true;
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
	
	// If this is a BOR event, then take ownership of
	// the DBORptrs object. If not, then copy a pointer
	// to the latest DBORptrs object into the event.
	if(pe->borptrs) borptrs_list.push_front(pe->borptrs);

	// Copy info for this parsed event into the JEvent
	event.SetJEventSource(this);
	event.SetEventNumber(pe->event_number);
	event.SetRunNumber(USER_RUN_NUMBER>0 ? USER_RUN_NUMBER:pe->run_number);
	event.SetRef(pe);

	// Set event status bits
	event.SetStatus(pe->event_status_bits);
	event.SetStatusBit(kSTATUS_EVIO);
	if( source_type == kFileSource ) event.SetStatusBit(kSTATUS_FROM_FILE);
	if( source_type == kETSource   ) event.SetStatusBit(kSTATUS_FROM_ET);

	// EPICS and BOR events are barrier events
	if(event.GetStatusBit(kSTATUS_EPICS_EVENT)) event.SetSequential();
	if(event.GetStatusBit(kSTATUS_BOR_EVENT  )) event.SetSequential();
	
	// Only add BOR events to physics events
	if(pe->borptrs==NULL)
		if(!borptrs_list.empty()) pe->borptrs = borptrs_list.front();
	
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
	
	// Get any translation tables we'll need to apply
	JEventLoop *loop = event.GetJEventLoop();
	vector<const DTranslationTable*> translationTables;
	DTranslationTable_factory *ttfac = static_cast<DTranslationTable_factory*>(loop->GetFactory("DTranslationTable"));
	if(ttfac) ttfac->Get(translationTables);
	
	// Copy pointers to all hits to appropriate factories.
	// Link BORconfig objects and apply translation tables if appropriate.
	DParsedEvent *pe = (DParsedEvent*)event.GetRef();
	if(!pe->copied_to_factories){

		// Copy all low-level hits to appropriate factories
		pe->CopyToFactories(loop);

		// Optionally link BOR object associations
		if(LINK && pe->borptrs) LinkBORassociations(pe);

		// Apply translation tables to create DigiHit objects
		for(auto tt : translationTables){
			tt->ApplyTranslationTable(loop);
		}
	}

	// Decide whether this is a data type the source supplies
	bool isSuppliedType = pe->IsParsedDataType(dataClassName);
	for(auto tt : translationTables){
		if(isSuppliedType) break;  // once this is set, it's set
		isSuppliedType = tt->IsSuppliedType(dataClassName);
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

//----------------------------
// LinkAssociationsBORConfigB
//----------------------------
template<class T, class U>
void LinkAssociationsBORConfigB(vector<T*> &a, vector<U*> &b)
{
	/// Template routine to loop over two vectors of pointers to
	/// objects derived from DDAQAddress. This will find any hits
	/// coming from the same DAQ module (channel number is not checked)
	/// When a match is found, the pointer from "a" will be added
	/// to "b"'s AssociatedObjects list. This will NOT do the inverse
	/// of adding "b" to "a"'s list. It is intended for adding a module
	/// level trigger time object to all hits from that module. Adding
	/// all of the hits to the trigger time object seems like it would
	/// be a little expensive with no real use case.
	///
	/// Note that this assumes the input vectors have been sorted by
	/// rocid, then slot in ascending order.

	// Bail early if nothing to link
	if(b.empty()) return;
	if(a.empty()) return;

	for(uint32_t i=0, j=0; i<b.size(); ){

		uint32_t rocid = b[i]->rocid;
		uint32_t slot  = b[i]->slot;

		// Find start and end of range in b
		uint32_t istart = i;
		uint32_t iend   = i+1; // index of first element outside of ROI
		for(; iend<b.size(); iend++){
			if(b[iend]->rocid != rocid) break;
			if(b[iend]->slot  != slot ) break;
		}
		i = iend; // setup for next iteration
		
		// Find start of range in a
		for(; j<a.size(); j++){
			if( a[j]->rocid > rocid ) break;
			if( a[j]->rocid == rocid ){
				if( a[j]->slot >= slot ) break;
			}
		}
		if(j>=a.size()) break; // exhausted all a's. we're done

		if( a[j]->rocid > rocid ) continue; // couldn't find rocid in a
		if( a[j]->slot  > slot  ) continue; // couldn't find slot in a
		
		// Find end of range in a
		uint32_t jend = j+1;
		for(; jend<a.size(); jend++){
			if(a[jend]->rocid != rocid) break;
			if(a[jend]->slot  != slot ) break;
		}

		// Loop over all combos of both ranges and make associations
		uint32_t jstart = j;
		for(uint32_t ii=istart; ii<iend; ii++){
			for(uint32_t jj=jstart; jj<jend; jj++){
				b[ii]->AddAssociatedObject(a[jj]);
			}
		}
		j = jend;

		if( i>=b.size() ) break;
		if( j>=a.size() ) break;
	}
}

//----------------
// LinkBORassociations
//----------------
void JEventSource_EVIOpp::LinkBORassociations(DParsedEvent *pe)
{
	/// Add BORConfig objects as associated objects
	/// to selected hit objects. Most other object associations
	/// are made in DEVIOWorkerThread::LinkAllAssociations. The
	/// BORConfig objects however, are not available when that
	/// is called.
	/// This is called from GetObjects() which is called from
	/// one of the processing threads.
	
	// n.b. all of the BORConfig object vectors will already
	// be sorted by rocid, then slot when the data was parsed.
	// The other hit objects are also already sorted (in
	// DEVIOWorkerThread::LinkAllAssociations).

	DBORptrs *borptrs = pe->borptrs;

	LinkAssociationsBORConfigB(borptrs->vDf250BORConfig, pe->vDf250WindowRawData);
	LinkAssociationsBORConfigB(borptrs->vDf250BORConfig, pe->vDf250PulseIntegral);

	LinkAssociationsBORConfigB(borptrs->vDf125BORConfig, pe->vDf125WindowRawData);
	LinkAssociationsBORConfigB(borptrs->vDf125BORConfig, pe->vDf125PulseIntegral);
	LinkAssociationsBORConfigB(borptrs->vDf125BORConfig, pe->vDf125CDCPulse);
	LinkAssociationsBORConfigB(borptrs->vDf125BORConfig, pe->vDf125FDCPulse);

	LinkAssociationsBORConfigB(borptrs->vDF1TDCBORConfig, pe->vDF1TDCHit);

	LinkAssociationsBORConfigB(borptrs->vDCAEN1290TDCBORConfig, pe->vDCAEN1290TDCHit);

}

////----------------
//// Cleanup
////----------------
//void JEventSource_EVIOpp::Cleanup(void)
//{
//	/// This is called internally by the JEventSource_EVIO class
//	/// once all events have been read in. Its purpose is to
//	/// free the hidden memory in all of the container class
//	/// members of the JEventSource_EVIO class. This is needed
//	/// for jobs that process a lot of input files and therefore
//	/// create a lot JEventSource_EVIO objects. JANA does not delete
//	/// these objects until the end of the job so this tends to
//	/// act like a memory leak. The data used can be substantial
//	/// (nearly 1GB per JEventSource_EVIO object).
//
//	if(hdevio) delete hdevio;
//	hdevio = NULL;
//
//	if(hdet) delete hdet;
//	hdet = NULL;
//
//}

