// $Id$
//
//    File: DEVIOWorkerThread.cc
// Created: Mon Mar 28 07:40:07 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#include <unistd.h>

#include "DEVIOWorkerThread.h"
#include "JEventSource_EVIOpp.h"
#include "LinkAssociations.h"

#include <swap_bank.h>

using namespace std;
using namespace std::chrono;



//---------------------------------
// DEVIOWorkerThread    (Constructor)
//---------------------------------
DEVIOWorkerThread::DEVIOWorkerThread(
	 JEventSource_EVIOpp  *event_source
	 ,list<DParsedEvent*> &parsed_events
	 ,uint32_t            &MAX_PARSED_EVENTS
	 ,mutex               &PARSED_EVENTS_MUTEX
	 ,condition_variable  &PARSED_EVENTS_CV
	 ):
	 event_source(event_source)
	,parsed_events(parsed_events)
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
	
	VERBOSE             = 1;
	Nrecycled           = 0;     // Incremented in JEventSource_EVIOpp::Dispatcher()
	MAX_EVENT_RECYCLES  = 1000;  // In EVIO events (not L1 trigger events!) overwritten in JEventSource_EVIOpp constructor
	MAX_OBJECT_RECYCLES = 1000;  // overwritten in JEventSource_EVIOpp constructor
	run_number_seed     = 0;     // Set in JEventSource_EVIOpp constructor

	in_use              = false;
	jobtype             = JOB_NONE;

	buff_len            = 100;   // this will grow as needed
	buff                = new uint32_t[buff_len];

	PARSE_F250          = true;
	PARSE_F125          = true;
	PARSE_F1TDC         = true;
	PARSE_CAEN1290TDC   = true;
	PARSE_CONFIG        = true;
	PARSE_BOR           = true;
	PARSE_EPICS         = true;
	PARSE_EVENTTAG      = true;
	PARSE_TRIGGER       = true;

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
		
		// In principle, in_use should never be false with a jobtype!=JOB_NONE
		// In practice, this has happened, possibly due to compiler optimization
		// reordering things in JEventSource_EVIOpp::Dispatcher. That led to 
		// attempting to process a buffer that was being written to. Avoid that
		// condition by checking the in_use flag is really set.
		if( !in_use ) continue;

		try {

			if( jobtype & JOB_SWAP       ) swap_bank(buff, buff, swap32(buff[0])+1 );

			if( jobtype & JOB_FULL_PARSE ) MakeEvents();
			
			if( jobtype & JOB_ASSOCIATE  ) LinkAllAssociations();
			
			if( !current_parsed_events.empty() ) PublishEvents();

		} catch (exception &e) {
			jerr << e.what() << endl;
			for(auto pe : parsed_event_pool) delete pe; // delete all parsed events any any objects they hold
			parsed_event_pool.clear();
			current_parsed_events.clear(); // (these are also in parsed_event_pool so were already deleted)
			//exit(-1);
		}
		
		// Reset and mark us as available for use
		jobtype = JOB_NONE;
		in_use  = false;

		if( jobtype & JOB_QUIT       ) break;		
	}

	in_use  = false;
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
// Prune
//---------------------------------
void DEVIOWorkerThread::Prune(void)
{
	/// Delete any DParsedEvent objects not currently in use. 
	/// If the DParsedEvent object pool and their internal
	/// hit object pools are allowed to continuously grow, it
	/// will appear as a though there is a memory leak. Occasional
	/// pruning will reduce the average memory footprint. 
	/// This is called from MakeEvents() every MAX_EVENT_RECYCLES
	/// EVIO events processed by this worker thread.
	/// Note that this is in EVIO events (i.e. possibly a block
	/// of events) not in L1 trigger events.
	///
	/// NOTE: We currently do NOT reduce the size of buff
	/// here if it is too big. We may wish to do that at some point!

	// Delete extra parsed events
	vector<DParsedEvent*> tmp_events = parsed_event_pool;
	parsed_event_pool.clear();
	for(auto pe : tmp_events) {
		if(pe->in_use)
			parsed_event_pool.push_back(pe);
		else
			delete pe;

	}
}

//---------------------------------
// MakeEvents
//---------------------------------
void DEVIOWorkerThread::MakeEvents(void)
{
	
	/// Make DParsedEvent objects from data currently in buff.
	/// This will look at the begining of the EVIO event to see
	/// how many L1 events are in it. It will then grab that many
	/// DParsedEvent objects from this threads pool , or create
	/// new ones and add them all to the current_parsed_events
	/// vector. These are then filled out later as the data is
	/// parsed.
	
	if(!current_parsed_events.empty()) throw JException("Attempting call to DEVIOWorkerThread::MakeEvents when current_parsed_events not empty!!", __FILE__, __LINE__);
	
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
	for(auto pe : parsed_event_pool){
		if(pe->in_use) continue;
		current_parsed_events.push_back(pe);
		if( current_parsed_events.size() >= M ) break;
	}
	
	// Create new DParsedEvent objects if needed
	while( current_parsed_events.size() < M ){
		DParsedEvent *pe = new DParsedEvent(MAX_OBJECT_RECYCLES);
		current_parsed_events.push_back(pe);
		parsed_event_pool.push_back(pe);
	}
	
	// Set indexes for the parsed event objects
	// and flag them as being in use.
	for(auto pe : current_parsed_events){
	
		pe->Clear(); // return previous event's objects to pools and clear vectors
		pe->istreamorder = istreamorder;
		pe->run_number   = run_number_seed;
		pe->event_number = event_num++;
		pe->sync_flag    = false;
		pe->in_use       = true;
		pe->copied_to_factories = false;
		pe->event_status_bits   = 0;
		pe->borptrs      = NULL; // may be set by either ParseBORbank or JEventSource_EVIOpp::GetEvent
	}

	// Parse data in buffer to create data objects
	ParseBank();
	
	// Occasionally prune extra DParsedEvent objects as well as objects
	// from the existing pools to reduce average memory usage. We do
	// this after parsing so that not everything is deleted (objects
	// being used this event will be returned to the pools later.)
	if(++Nrecycled%MAX_EVENT_RECYCLES == 0) Prune();
	for(auto pe : current_parsed_events){
		if( ++pe->Nrecycled%pe->MAX_RECYCLES == 0) pe->Prune();
	}
}	

//---------------------------------
// PublishEvents
//---------------------------------
void DEVIOWorkerThread::PublishEvents(void)
{	
	/// Copy our "current_parsed_events" pointers into the global "parsed_events"
	/// list making them available for consumption. 
	
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
	
	// Any events should now be published
	current_parsed_events.clear();
}

//---------------------------------
// ParseBank
//---------------------------------
void DEVIOWorkerThread::ParseBank(void)
{

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

			case 0xFFD0:
			case 0xFFD1:
			case 0xFFD2:
			case 0xFFD3:    ParseControlEvent(iptr, iend);    break;

			case 0xFF58:
			case 0xFF78: current_parsed_events.back()->sync_flag = true;
			case 0xFF50:     
			case 0xFF70:     ParsePhysicsBank(iptr, iend);    break;

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
	if(!PARSE_EPICS){ iptr = iend; return; }

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
	pe->event_status_bits |= (1<<kSTATUS_EPICS_EVENT);
	
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
			pe->NEW_DEPICSvalue(timestamp, nameval);
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
	/// Create BOR config objects from the EVIO bank and store them in
	/// the event (should only be one since BOR events are not entangled).
	/// These objects will eventually be inherited by the JEventSource_EVIOpp
	/// object and passed to all subsequent events.

	// Upon entry, iptr should point to length word of a bank of banks with tag=0x70
	// indicating BOR event. Each bank contained within will represent one crate and
	// will be a bank with tag=0x71 and num the rocid, containing tagsegments. Each tagsegment
	// represents a single module with the tag containing the module type (bits 0-4) and
	// slot (bits 5-10). The data in the tagsegments is uint32_t and maps to a data
	// structure in bor_roc.h depending on the module type. Below is a summary of
	// how this looks in memory:
	//
	//	BOR event length
	//	BOR header
	//	crate bank length
	//	crate header
	//	module bank len/header
	//	module data ...
	//	module bank len/header
	//	module data ...
	//  ...
	//	crate bank length
	//	crate header
	//	...

	if(!PARSE_BOR){ iptr = &iptr[(*iptr) + 1]; return; }

	// Make sure there is exactly 1 event in current_parsed_events
	if(current_parsed_events.size() != 1){
		stringstream ss;
		ss << "DEVIOWorkerThread::ParseBORbank called for EVIO event with " << current_parsed_events.size() << " events in it. (Should be exactly 1!)";
		throw JException(ss.str(), __FILE__, __LINE__);
	}
	
	// Create new DBORptrs object and set pointer to it in DParsedEvent
	// (see JEventSource_EVIOpp::GetEvent)
	DParsedEvent *pe = current_parsed_events.front();
	pe->event_status_bits |= (1<<kSTATUS_BOR_EVENT);
	pe->borptrs = new DBORptrs();
	DBORptrs* &borptrs = pe->borptrs;
	
	// Make sure we have full event
	uint32_t borevent_len = *iptr++;
	uint32_t bank_len = (uint32_t)((uint64_t)iend - (uint64_t)iptr)/sizeof(uint32_t);
	if(borevent_len > bank_len){
		stringstream ss;
		ss << "BOR: Size of bank doesn't match amount of data given (" << borevent_len << " > " << bank_len << ")";
		throw JException(ss.str(), __FILE__, __LINE__);
	}
	iend = &iptr[borevent_len]; // in case they give us too much data!
	
	// Make sure BOR header word is right
	uint32_t bor_header = *iptr++;
	if(bor_header != 0x700e01){
		stringstream ss;
		ss << "Bad BOR header: 0x" << hex << bor_header;
		throw JException(ss.str(), __FILE__, __LINE__);
	}

	// Loop over crates
	while(iptr<iend){
		uint32_t crate_len    = *iptr++;
		uint32_t *iend_crate  = &iptr[crate_len]; // points to first word after this crate
		uint32_t crate_header = *iptr++;
//		uint32_t rocid = crate_header&0xFF;
		
		// Make sure crate tag is right
		if( (crate_header>>16) != 0x71 ){
			stringstream ss;
			ss << "Bad BOR crate header: 0x" << hex << (crate_header>>16);
			throw JException(ss.str(), __FILE__, __LINE__);
		}

		// Loop over modules
		while(iptr<iend_crate){
			uint32_t module_header = *iptr++;
			uint32_t module_len    = module_header&0xFFFF;
			uint32_t modType       = (module_header>>20)&0x1f;
//			uint32_t slot          = (module_header>>25);
//			uint32_t *iend_module  = &iptr[module_len]; // points to first word after this module

			uint32_t *src          = iptr;
			uint32_t *dest         = NULL;
			uint32_t sizeof_dest   = 0;

			Df250BORConfig *f250conf = NULL;
			Df125BORConfig *f125conf = NULL;
			DF1TDCBORConfig *F1TDCconf = NULL;
			DCAEN1290TDCBORConfig *caen1190conf = NULL;

			switch(modType){
				case DModuleType::FADC250: // f250
					f250conf = new Df250BORConfig;
					dest = (uint32_t*)&f250conf->rocid;
					sizeof_dest = sizeof(f250config);
					break;
				case DModuleType::FADC125: // f125
					f125conf = new Df125BORConfig;
					dest = (uint32_t*)&f125conf->rocid;
					sizeof_dest = sizeof(f125config);
					break;

				case DModuleType::F1TDC32: // F1TDCv2
				case DModuleType::F1TDC48: // F1TDCv3
					F1TDCconf = new DF1TDCBORConfig;
					dest = (uint32_t*)&F1TDCconf->rocid;
					sizeof_dest = sizeof(F1TDCconfig);
					break;

				case DModuleType::CAEN1190: // CAEN 1190 TDC
				case DModuleType::CAEN1290: // CAEN 1290 TDC
					caen1190conf = new DCAEN1290TDCBORConfig;
					dest = (uint32_t*)&caen1190conf->rocid;
					sizeof_dest = sizeof(caen1190config);
					break;
				
				default:
					{
					stringstream ss;
					ss << "Unknown BOR module type: " << modType << "  (module_header=0x"<<hex<<module_header<<")";
					jerr << ss.str() << endl;
					throw JException(ss.str(), __FILE__, __LINE__);
					}
			}

			// Check that the bank size and data structure size match.
			if( module_len != (sizeof_dest/sizeof(uint32_t)) ){
				stringstream ss;
				ss << "BOR module bank size does not match structure! " << module_len << " != " << (sizeof_dest/sizeof(uint32_t)) << " for modType " << modType;
				throw JException(ss.str(), __FILE__, __LINE__);
			}

			// Copy bank data, assuming format is the same
			for(uint32_t i=0; i<module_len; i++) *dest++ = *src++;

			// Store object for use in this and subsequent events
			if(f250conf    ) borptrs->vDf250BORConfig.push_back(f250conf);
			if(f125conf    ) borptrs->vDf125BORConfig.push_back(f125conf);
			if(F1TDCconf   ) borptrs->vDF1TDCBORConfig.push_back(F1TDCconf);
			if(caen1190conf) borptrs->vDCAEN1290TDCBORConfig.push_back(caen1190conf);

			iptr = &iptr[module_len];
		}
		
		iptr = iend_crate; // ensure we're pointing past this crate
	}
	
	// Sort the BOR config events now so we don't have to do it for every event
	borptrs->Sort();

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
// ParseControlEvent
//---------------------------------
void DEVIOWorkerThread::ParseControlEvent(uint32_t* &iptr, uint32_t *iend)
{
	for(auto pe : current_parsed_events) pe->event_status_bits |= (1<<kSTATUS_CONTROL_EVENT);

	iptr = &iptr[(*iptr) + 1];
}

//---------------------------------
// ParsePhysicsBank
//---------------------------------
void DEVIOWorkerThread::ParsePhysicsBank(uint32_t* &iptr, uint32_t *iend)
{

	for(auto pe : current_parsed_events) pe->event_status_bits |= (1<<kSTATUS_PHYSICS_EVENT);

	uint32_t physics_event_len      = *iptr++;
	uint32_t *iend_physics_event    = &iptr[physics_event_len];
	iptr++;

	// Built Trigger Bank
	uint32_t built_trigger_bank_len  = *iptr;
	uint32_t *iend_built_trigger_bank = &iptr[built_trigger_bank_len+1];
	ParseBuiltTriggerBank(iptr, iend_built_trigger_bank);
	iptr = iend_built_trigger_bank;
	
	// Loop over Data banks
	while( iptr < iend_physics_event ) {

		uint32_t data_bank_len = *iptr;
		uint32_t *iend_data_bank = &iptr[data_bank_len+1];

		ParseDataBank(iptr, iend_data_bank);

		iptr = iend_data_bank;
	}

	iptr = iend_physics_event;
}

//---------------------------------
// ParseBuiltTriggerBank
//---------------------------------
void DEVIOWorkerThread::ParseBuiltTriggerBank(uint32_t* &iptr, uint32_t *iend)
{
	if(!PARSE_TRIGGER) return;

	iptr++; // advance past length word
	uint32_t mask = 0xFF202000;
	if( ((*iptr) & mask) != mask ){
		stringstream ss;
		ss << "Bad header word in Built Trigger Bank: " << hex << *iptr;
		throw JException(ss.str(), __FILE__, __LINE__);
	}
	
	uint32_t tag     = (*iptr)>>16; // 0xFF2X
	uint32_t Nrocs   = (*iptr++) & 0xFF;
	uint32_t Mevents = current_parsed_events.size();
	
	//-------- Common data (64bit)
	uint32_t common_header64 = *iptr++;
	uint32_t common_header64_len = common_header64 & 0xFFFF;
	uint64_t *iptr64 = (uint64_t*)iptr;
	iptr = &iptr[common_header64_len];

	// First event number
   uint64_t first_event_num = *iptr64++;

   // Hi and lo 32bit words in 64bit numbers seem to be
   // switched for events read from ET, but not read from
   // file. Not sure if this is in the swapping routine
//   if(source_type==kETSource) first_event_num = (first_event_num>>32) | (first_event_num<<32);

	// Average timestamps
   uint32_t Ntimestamps = (common_header64_len/2)-1;
   if(tag & 0x2) Ntimestamps--; // subtract 1 for run number/type word if present
	vector<uint64_t> avg_timestamps;
   for(uint32_t i=0; i<Ntimestamps; i++) avg_timestamps.push_back(*iptr64++);

   // run number and run type
	uint32_t run_number = 0;
	uint32_t run_type   = 0;
   if(tag & 0x02){
       run_number = (*iptr64) >> 32;
       run_type   = (*iptr64) & 0xFFFFFFFF;
		 iptr64++;
   }

	//-------- Common data (16bit)
	uint32_t common_header16 = *iptr++;
	uint32_t common_header16_len = common_header16 & 0xFFFF;
	uint16_t *iptr16 = (uint16_t*)iptr;
	iptr = &iptr[common_header16_len];

	vector<uint16_t> event_types;
   for(uint32_t i=0; i<Mevents; i++) event_types.push_back(*iptr16++);
	
	//-------- ROC data (32bit)
	for(uint32_t iroc=0; iroc<Nrocs; iroc++){
		uint32_t common_header32 = *iptr++;
		uint32_t common_header32_len = common_header32 & 0xFFFF;
		uint32_t rocid = common_header32 >> 24;

		uint32_t Nwords_per_event = common_header32_len/Mevents;
		for(auto pe : current_parsed_events){

			DCODAROCInfo *codarocinfo = pe->NEW_DCODAROCInfo();
			codarocinfo->rocid = rocid;

			uint64_t ts_low  = *iptr++;
			uint64_t ts_high = *iptr++;
			codarocinfo->timestamp = (ts_high<<32) + ts_low;
			codarocinfo->misc.clear(); // could be recycled from previous event
			for(uint32_t i=2; i<Nwords_per_event; i++) codarocinfo->misc.push_back(*iptr++);
			
			if(iptr > iend){
				throw JException("Bad data format in ParseBuiltTriggerBank!", __FILE__, __LINE__);
			}
		}
	}
	
	//-------- Make DCODAEventInfo objects
	uint64_t ievent = 0;
	for(auto pe : current_parsed_events){

		pe->run_number = run_number; // may be overwritten in JEventSource_EVIOpp::GetEvent()

		DCODAEventInfo *codaeventinfo = pe->NEW_DCODAEventInfo();
		codaeventinfo->run_number     = run_number;
		codaeventinfo->run_type       = run_type;
		codaeventinfo->event_number   = first_event_num + ievent;
		codaeventinfo->event_type     = event_types.empty() ? 0:event_types[ievent];
		codaeventinfo->avg_timestamp  = avg_timestamps.empty() ? 0:avg_timestamps[ievent];
		ievent++;
	}
}

//---------------------------------
// ParseDataBank
//---------------------------------
void DEVIOWorkerThread::ParseDataBank(uint32_t* &iptr, uint32_t *iend)
{
	// Physics Event's Data Bank header
	iptr++; // advance past data bank length word
	uint32_t rocid = ((*iptr)>>16) & 0xFFF;
	iptr++;
	
	// Loop over Data Block Banks
	while(iptr < iend){
		
		uint32_t data_block_bank_len     = *iptr++;
		uint32_t *iend_data_block_bank   = &iptr[data_block_bank_len];
		uint32_t data_block_bank_header  = *iptr++;
		
		// Not sure where this comes from, but it needs to be skipped if present
		while( (*iptr==0xF800FAFA) && (iptr<iend) ) iptr++;
		
		uint32_t det_id = (data_block_bank_header>>16) & 0xFFF;
        switch(det_id){

            case 20:
                ParseCAEN1190(rocid, iptr, iend_data_block_bank);
                break;

            case 0x55:
                ParseModuleConfiguration(rocid, iptr, iend_data_block_bank);
                break;

            case 0:
            case 1:
            case 3:
            case 6:  // flash 250 module, MMD 2014/2/4
            case 16: // flash 125 module (CDC), DL 2014/6/19
            case 26: // F1 TDC module (BCAL), MMD 2014-07-31
                ParseJLabModuleData(rocid, iptr, iend_data_block_bank);
                break;

			// These were implemented in the ROL for sync events
			// as 0xEE02 and 0xEE05. However, that violates the
			// spec. which reserves the top 4 bits as status bits
			// (the first "E" should really be a "1". We just check
			// other 12 bits here.
			case 0xE02:
					ParseTSscalerBank(iptr, iend);
					break;
			case 0xE05:
					Parsef250scalerBank(iptr, iend);
					break;
			case 0xE10:  // really wish Sascha would share when he does this stuff!
					Parsef250scalerBank(iptr, iend);
					break;

			case 5:
				// old ROL Beni used had this but I don't think its
				// been used for years. Run 10390 seems to have
				// this though (???)
				break;


			default:
				jerr<<"Unknown module type ("<<det_id<<" = " << hex << det_id << dec << " ) encountered" << endl;
//				if(VERBOSE>5){
					cout << "----- First few words to help with debugging -----" << endl;
					cout.flush(); cerr.flush();
					DumpBinary(&iptr[-2], iend, 32, &iptr[-1]);
//				}
		}

		iptr = iend_data_block_bank;
	}
	
}

//----------------
// ParseCAEN1190
//----------------
void DEVIOWorkerThread::ParseCAEN1190(uint32_t rocid, uint32_t* &iptr, uint32_t *iend)
{
	if(!PARSE_CAEN1290TDC){ iptr = &iptr[(*iptr) + 1]; return; }

    /// Parse data from a CAEN 1190 or 1290 module
    /// (See ppg. 72-74 of V1290_REV15.pdf manual)

    uint32_t slot = 0;
    uint32_t event_count = 0;
    uint32_t word_count = 0;
    uint32_t trigger_time_tag = 0;
    uint32_t tdc_num = 0;
    uint32_t event_id = 0;
    uint32_t bunch_id = 0;

    // We need to accomodate multi-event blocks where
    // events are entangled (i.e. hits from event 1
    // are mixed in between those of event 2,3,4,
    // etc... With CAEN modules, we only know which
    // event a hit came from by looking at the event_id
    // in the TDC header. This value is only 12 bits
    // and could roll over within an event block. This
    // means we need to keep track of the order we
    // encounter them in so it is maintained in the
    // "events" container. The event_id order is kept
    // in the "event_id_order" vector.
	map<uint32_t, DParsedEvent*> events_by_event_id;

	auto pe_iter = current_parsed_events.begin();
	DParsedEvent *pe = NULL;

    while(iptr<iend){

        // This word appears to be appended to the data.
        // Probably in the ROL. Ignore it if found.
        if(*iptr == 0xd00dd00d) {
            if(VERBOSE>7) cout << "         CAEN skipping 0xd00dd00d word" << endl;
            iptr++;
            continue;
        }

        uint32_t type = (*iptr) >> 27;
        uint32_t edge = 0; // 1=trailing, 0=leading
        uint32_t channel = 0;
        uint32_t tdc = 0;
        uint32_t error_flags = 0;
        switch(type){
            case 0b01000:  // Global Header
                slot = (*iptr) & 0x1f;
                event_count = ((*iptr)>>5) & 0xffffff;
                if(VERBOSE>7) cout << "         CAEN TDC Global Header (slot=" << slot << " , event count=" << event_count << ")" << endl;
                break;
            case 0b10000:  // Global Trailer
                slot = (*iptr) & 0x1f;
                word_count = ((*iptr)>>5) & 0x7ffff;
                if(VERBOSE>7) cout << "         CAEN TDC Global Trailer (slot=" << slot << " , word count=" << word_count << ")" << endl;
                slot = event_count = word_count = trigger_time_tag = tdc_num = event_id = bunch_id = 0;
                break;
            case 0b10001:  // Global Trigger Time Tag
                trigger_time_tag = ((*iptr)>>5) & 0x7ffffff;
                if(VERBOSE>7) cout << "         CAEN TDC Global Trigger Time Tag (tag=" << trigger_time_tag << ")" << endl;
                break;
            case 0b00001:  // TDC Header
                tdc_num = ((*iptr)>>24) & 0x03;
                event_id = ((*iptr)>>12) & 0x0fff;
                bunch_id = (*iptr) & 0x0fff;
				if(events_by_event_id.find(event_id) == events_by_event_id.end()){
					if(pe_iter == current_parsed_events.end()){
						_DBG_ << "CAEN1290TDC parser sees more events than CODA header! (>" << current_parsed_events.size() << ")" << endl;
						for( auto p : events_by_event_id) cout << "id=" << p.first << endl;
						iptr = iend;
						exit(-1); // should we exit, or try and continue??
						return;
					}
					pe = *pe_iter++;
					events_by_event_id[event_id] = pe;
				}else{
					pe = events_by_event_id[event_id];
				}				
                if(VERBOSE>7) cout << "         CAEN TDC TDC Header (tdc=" << tdc_num <<" , event id=" << event_id <<" , bunch id=" << bunch_id << ")" << endl;
                break;
            case 0b00000:  // TDC Measurement
                edge = ((*iptr)>>26) & 0x01;
                channel = ((*iptr)>>21) & 0x1f;
                tdc = ((*iptr)>>0) & 0x1fffff;
                if(VERBOSE>7) cout << "         CAEN TDC TDC Measurement (" << (edge ? "trailing":"leading") << " , channel=" << channel << " , tdc=" << tdc << ")" << endl;

                // Create DCAEN1290TDCHit object
                pe->NEW_DCAEN1290TDCHit(rocid, slot, channel, 0, edge, tdc_num, event_id, bunch_id, tdc);
                break;
            case 0b00100:  // TDC Error
                error_flags = (*iptr) & 0x7fff;
                if(VERBOSE>7) cout << "         CAEN TDC TDC Error (err flags=0x" << hex << error_flags << dec << ")" << endl;
                break;
            case 0b00011:  // TDC Trailer
                tdc_num = ((*iptr)>>24) & 0x03;
                event_id = ((*iptr)>>12) & 0x0fff;
                word_count = ((*iptr)>>0) & 0x0fff;
                if(VERBOSE>7) cout << "         CAEN TDC TDC Trailer (tdc=" << tdc_num <<" , event id=" << event_id <<" , word count=" << word_count << ")" << endl;
                tdc_num = event_id = bunch_id = 0;
                break;
            case 0b11000:  // Filler Word
                if(VERBOSE>7) cout << "         CAEN TDC Filler Word" << endl;
                break;
            default:
                cout << "Unknown datatype: 0x" << hex << type << " full word: "<< *iptr << dec << endl;
        }

        iptr++;
    }

}

//----------------
// ParseModuleConfiguration
//----------------
void DEVIOWorkerThread::ParseModuleConfiguration(uint32_t rocid, uint32_t* &iptr, uint32_t *iend)
{
	if(!PARSE_CONFIG){ iptr = &iptr[(*iptr) + 1]; return; }

    /// Parse a bank of module configuration data. These are configuration values
    /// programmed into the module at the beginning of the run that may be needed
    /// in the offline. For example, the number of samples to sum in a FADC pulse
    /// integral.
    ///
    /// The bank has one or more sections, each describing parameters applicable 
    /// to a number of modules as indicated by a 24bit slot mask.
    ///
    /// This bank should appear only once per DAQ event which, if in multi-event
    /// block mode, may have multiple L1 events. The parameters here will apply
    /// to all L1 events in the block. This method will put the config objects
	/// into each event in current_parsed_events. The config objects are duplicated
	/// as needed so each event has its own, indepenent set of config object.

    while(iptr < iend){
        uint32_t slot_mask = (*iptr) & 0xFFFFFF;
        uint32_t Nvals = ((*iptr) >> 24) & 0xFF;
        iptr++;

		// Events will be created in the first event (i.e. using its pool)
		// but pointers are saved so we can use them to construct identical
		// objects in all other event later
		DParsedEvent *pe = current_parsed_events.front();

        Df250Config *f250config = NULL;
        Df125Config *f125config = NULL;
        DF1TDCConfig *f1tdcconfig = NULL;
        DCAEN1290TDCConfig *caen1290tdcconfig = NULL;

        // Loop over all parameters in this section
        for(uint32_t i=0; i< Nvals; i++){
            if( iptr >= iend){
                _DBG_ << "DAQ Configuration bank corrupt! slot_mask=0x" << hex << slot_mask << dec << " Nvals="<< Nvals << endl;
                exit(-1);
            }

            daq_param_type ptype = (daq_param_type)((*iptr)>>16);
            uint16_t val = (*iptr) & 0xFFFF;

            if(VERBOSE>6) cout << "       DAQ parameter of type: 0x" << hex << ptype << dec << "  found with value: " << val << endl;

            // Create config object of correct type if needed and copy
            // parameter value into it.
            switch(ptype>>8){

                // f250
                case 0x05:
                    if( !f250config ) f250config = pe->NEW_Df250Config(rocid, slot_mask);
                    switch(ptype){
                        case kPARAM250_NSA            : f250config->NSA              = val; break;
                        case kPARAM250_NSB            : f250config->NSB              = val; break;
                        case kPARAM250_NSA_NSB        : f250config->NSA_NSB          = val; break;
                        case kPARAM250_NPED           : f250config->NPED             = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                    // f125
                case 0x0F:
                    if( !f125config ) f125config = pe->NEW_Df125Config(rocid, slot_mask);
                    switch(ptype){
                        case kPARAM125_NSA            : f125config->NSA              = val; break;
                        case kPARAM125_NSB            : f125config->NSB              = val; break;
                        case kPARAM125_NSA_NSB        : f125config->NSA_NSB          = val; break;
                        case kPARAM125_NPED           : f125config->NPED             = val; break;
                        case kPARAM125_WINWIDTH       : f125config->WINWIDTH         = val; break;
                        case kPARAM125_PL             : f125config->PL               = val; break;
                        case kPARAM125_NW             : f125config->NW               = val; break;
                        case kPARAM125_NPK            : f125config->NPK              = val; break;
                        case kPARAM125_P1             : f125config->P1               = val; break;
                        case kPARAM125_P2             : f125config->P2               = val; break;
                        case kPARAM125_PG             : f125config->PG               = val; break;
                        case kPARAM125_IE             : f125config->IE               = val; break;
                        case kPARAM125_H              : f125config->H                = val; break;
                        case kPARAM125_TH             : f125config->TH               = val; break;
                        case kPARAM125_TL             : f125config->TL               = val; break;
                        case kPARAM125_IBIT           : f125config->IBIT             = val; break;
                        case kPARAM125_ABIT           : f125config->ABIT             = val; break;
                        case kPARAM125_PBIT           : f125config->PBIT             = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                    // F1TDC
                case 0x06:
                    if( !f1tdcconfig ) f1tdcconfig = pe->NEW_DF1TDCConfig(rocid, slot_mask);
                    switch(ptype){
                        case kPARAMF1_REFCNT          : f1tdcconfig->REFCNT          = val; break;
                        case kPARAMF1_TRIGWIN         : f1tdcconfig->TRIGWIN         = val; break;
                        case kPARAMF1_TRIGLAT         : f1tdcconfig->TRIGLAT         = val; break;
                        case kPARAMF1_HSDIV           : f1tdcconfig->HSDIV           = val; break;
                        case kPARAMF1_BINSIZE         : f1tdcconfig->BINSIZE         = val; break;
                        case kPARAMF1_REFCLKDIV       : f1tdcconfig->REFCLKDIV       = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                    // caen1290
                case 0x10:
                    if( !caen1290tdcconfig ) caen1290tdcconfig = pe->NEW_DCAEN1290TDCConfig(rocid, slot_mask);
                    switch(ptype){
                        case kPARAMCAEN1290_WINWIDTH  : caen1290tdcconfig->WINWIDTH  = val; break;
                        case kPARAMCAEN1290_WINOFFSET : caen1290tdcconfig->WINOFFSET = val; break;
                        default: _DBG_ << "UNKNOWN DAQ Config Parameter type: 0x" << hex << ptype << dec << endl;
                    }
                    break;

                default:
                    _DBG_ << "Unknown module type: 0x" << hex << (ptype>>8) << endl;
                    exit(-1);
            }


            iptr++;
        }

		// Make copies of all config objects for all other events
		for(auto tpe : current_parsed_events){
		
			if(tpe == pe) continue; // first event already owns objects so skip it

        	if(f250config       ) tpe->NEW_Df250Config(f250config);
        	if(f125config       ) tpe->NEW_Df125Config(f125config);
        	if(f1tdcconfig      ) tpe->NEW_DF1TDCConfig(f1tdcconfig);
        	if(caen1290tdcconfig) tpe->NEW_DCAEN1290TDCConfig(caen1290tdcconfig);
		}
    }
}

//----------------
// ParseJLabModuleData
//----------------
void DEVIOWorkerThread::ParseJLabModuleData(uint32_t rocid, uint32_t* &iptr, uint32_t *iend)
{
	
	while(iptr<iend){
	
		// Get module type from next word (bits 18-21)
		uint32_t mod_id = ((*iptr) >> 18) & 0x000F;
		MODULE_TYPE type = (MODULE_TYPE)mod_id;
//		cout << "      rocid=" << rocid << "  Encountered module type: " << type << " (=" << DModuleType::GetModule(type).GetName() << ")  word=" << hex << (*iptr) << dec << endl;

        switch(type){
            case DModuleType::FADC250:
                Parsef250Bank(rocid, iptr, iend);
                break;

            case DModuleType::FADC125:
                Parsef125Bank(rocid, iptr, iend);
                break;

            case DModuleType::F1TDC32:
                ParseF1TDCBank(rocid, iptr, iend);
                break;

            case DModuleType::F1TDC48:
                ParseF1TDCBank(rocid, iptr, iend);
                break;

           case DModuleType::TID:
//                ParseTIBank(rocid, iptr, iend);
                break;

            case DModuleType::UNKNOWN:
            default:
                jerr<<"Unknown module type ("<<mod_id<<") iptr=0x" << hex << iptr << dec << endl;

                while(iptr<iend && ((*iptr) & 0xF8000000) != 0x88000000) iptr++; // Skip to JLab block trailer
                iptr++; // advance past JLab block trailer
                while(iptr<iend && *iptr == 0xF8000000) iptr++; // skip filler words after block trailer
                break;
        }
	}

}

//----------------
// Parsef250Bank
//----------------
void DEVIOWorkerThread::Parsef250Bank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend)
{
	if(!PARSE_F250){ iptr = &iptr[(*iptr) + 1]; return; }

	auto pe_iter = current_parsed_events.begin();
	DParsedEvent *pe = NULL;
	
	uint32_t slot = 0;
	uint32_t itrigger = -1;

    // Loop over data words
    for(; iptr<iend; iptr++){

        // Skip all non-data-type-defining words at this
        // level. When we do encounter one, the appropriate
        // case block below should handle parsing all of
        // the data continuation words and advance the iptr.
        if(((*iptr>>31) & 0x1) == 0)continue;

        uint32_t data_type = (*iptr>>27) & 0x0F;
        switch(data_type){
            case 0: // Block Header
                slot = (*iptr>>22) & 0x1F;
                if(VERBOSE>7) cout << "      FADC250 Block Header: slot="<<slot<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                break;
            case 1: // Block Trailer
                pe_iter = current_parsed_events.begin();
				pe = NULL;
                if(VERBOSE>7) cout << "      FADC250 Block Trailer"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                break;
            case 2: // Event Header
                itrigger = (*iptr>>0) & 0x3FFFFF;
				pe = *pe_iter++;
                if(VERBOSE>7) cout << "      FADC250 Event Header: itrigger="<<itrigger<<", rocid="<<rocid<<", slot="<<slot<<")" <<" ("<<hex<<*iptr<<dec<<")" <<endl;
                break;
            case 3: // Trigger Time
				{
					uint64_t t = ((*iptr)&0xFFFFFF)<<0;
					iptr++;
					if(((*iptr>>31) & 0x1) == 0){
						t += ((*iptr)&0xFFFFFF)<<24; // from word on the street: second trigger time word is optional!!??
						if(VERBOSE>7) cout << "       Trigger time high word="<<(((*iptr)&0xFFFFFF))<<" ("<<hex<<*iptr<<dec<<")"<<endl;
					}else{
						iptr--;
					}
					if(VERBOSE>7) cout << "      FADC250 Trigger Time: t="<<t<<" ("<<hex<<*iptr<<dec<<")"<<endl;
					if(pe) pe->NEW_Df250TriggerTime(rocid, slot, itrigger, t);
				}
                break;
            case 4: // Window Raw Data
                // iptr passed by reference and so will be updated automatically
                if(VERBOSE>7) cout << "      FADC250 Window Raw Data"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                MakeDf250WindowRawData(pe, rocid, slot, itrigger, iptr);
                break;
            case 5: // Window Sum
				{
					uint32_t channel = (*iptr>>23) & 0x0F;
					uint32_t sum = (*iptr>>0) & 0x3FFFFF;
					uint32_t overflow = (*iptr>>22) & 0x1;
					if(VERBOSE>7) cout << "      FADC250 Window Sum"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
					if(pe) pe->NEW_Df250WindowSum(rocid, slot, channel, itrigger, sum, overflow);
				}
                break;				
            case 6: // Pulse Raw Data
//                MakeDf250PulseRawData(objs, rocid, slot, itrigger, iptr);
                if(VERBOSE>7) cout << "      FADC250 Pulse Raw Data"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                break;
            case 7: // Pulse Integral
				{
					uint32_t channel = (*iptr>>23) & 0x0F;
					uint32_t pulse_number = (*iptr>>21) & 0x03;
					uint32_t quality_factor = (*iptr>>19) & 0x03;
					uint32_t sum = (*iptr>>0) & 0x7FFFF;
					uint32_t nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
					uint32_t nsamples_pedestal = 1;  // The firmware returns an already divided pedestal
					uint32_t pedestal = 0;  // This will be replaced by the one from Df250PulsePedestal in GetObjects
					if(VERBOSE>7) cout << "      FADC250 Pulse Integral: chan="<<channel<<" pulse_number="<<pulse_number<<" sum="<<sum<<" ("<<hex<<*iptr<<dec<<")"<<endl;
					if(pe) pe->NEW_Df250PulseIntegral(rocid, slot, channel, itrigger, pulse_number, quality_factor, sum, pedestal, nsamples_integral, nsamples_pedestal);
				}
                break;
            case 8: // Pulse Time
				{
					uint32_t channel = (*iptr>>23) & 0x0F;
					uint32_t pulse_number = (*iptr>>21) & 0x03;
					uint32_t quality_factor = (*iptr>>19) & 0x03;
					uint32_t pulse_time = (*iptr>>0) & 0x7FFFF;
					if(VERBOSE>7) cout << "      FADC250 Pulse Time: chan="<<channel<<" pulse_number="<<pulse_number<<" pulse_time="<<pulse_time<<" ("<<hex<<*iptr<<dec<<")"<<endl;
					if(pe) pe->NEW_Df250PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time);
				}
				break;
            case 9: // Streaming Raw Data
                // This is marked "reserved for future implementation" in the current manual (v2).
                // As such, we don't try handling it here just yet.
                if(VERBOSE>7) cout << "      FADC250 Streaming Raw Data (unsupported)"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                break;
            case 10: // Pulse Pedestal
				{
					uint32_t channel = (*iptr>>23) & 0x0F;
					uint32_t pulse_number = (*iptr>>21) & 0x03;
					uint32_t pedestal = (*iptr>>12) & 0x1FF;
					uint32_t pulse_peak = (*iptr>>0) & 0xFFF;
					if(VERBOSE>7) cout << "      FADC250 Pulse Pedestal chan="<<channel<<" pulse_number="<<pulse_number<<" pedestal="<<pedestal<<" pulse_peak="<<pulse_peak<<" ("<<hex<<*iptr<<dec<<")"<<endl;
					if(pe) pe->NEW_Df250PulsePedestal(rocid, slot, channel, itrigger, pulse_number, pedestal, pulse_peak);
				}
                break;
            case 13: // Event Trailer
                // This is marked "suppressed for normal readout – debug mode only" in the
                // current manual (v2). It does not contain any data so the most we could do here
                // is return early. I'm hesitant to do that though since it would mean
                // different behavior for debug mode data as regular data.
            case 14: // Data not valid (empty module)
            case 15: // Filler (non-data) word
                if(VERBOSE>7) cout << "      FADC250 Event Trailer, Data not Valid, or Filler word ("<<data_type<<")"<<" ("<<hex<<*iptr<<dec<<")"<<endl;
                break;
        }
    }

    // Chop off filler words
    for(; iptr<iend; iptr++){
        if(((*iptr)&0xf8000000) != 0xf8000000) break;
    }
}

//----------------
// MakeDf250WindowRawData
//----------------
void DEVIOWorkerThread::MakeDf250WindowRawData(DParsedEvent *pe, uint32_t rocid, uint32_t slot, uint32_t itrigger, uint32_t* &iptr)
{
    uint32_t channel = (*iptr>>23) & 0x0F;
    uint32_t window_width = (*iptr>>0) & 0x0FFF;

    Df250WindowRawData *wrd = pe->NEW_Df250WindowRawData(rocid, slot, channel, itrigger);

    for(uint32_t isample=0; isample<window_width; isample +=2){

        // Advance to next word
        iptr++;

        // Make sure this is a data continuation word, if not, stop here
        if(((*iptr>>31) & 0x1) != 0x0){
            iptr--; // calling method expects us to point to last word in block
            break;
        }

        bool invalid_1 = (*iptr>>29) & 0x1;
        bool invalid_2 = (*iptr>>13) & 0x1;
        uint16_t sample_1 = 0;
        uint16_t sample_2 = 0;
        if(!invalid_1)sample_1 = (*iptr>>16) & 0x1FFF;
        if(!invalid_2)sample_2 = (*iptr>>0) & 0x1FFF;

        // Sample 1
        wrd->samples.push_back(sample_1);
        wrd->invalid_samples |= invalid_1;
        wrd->overflow |= (sample_1>>12) & 0x1;

        if(((isample+2) == window_width) && invalid_2)break; // skip last sample if flagged as invalid

        // Sample 2
        wrd->samples.push_back(sample_2);
        wrd->invalid_samples |= invalid_2;
        wrd->overflow |= (sample_2>>12) & 0x1;
    }
}

//----------------
// Parsef125Bank
//----------------
void DEVIOWorkerThread::Parsef125Bank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend)
{
	if(!PARSE_F125){ iptr = &iptr[(*iptr) + 1]; return; }

	auto pe_iter = current_parsed_events.begin();
	DParsedEvent *pe = NULL;

    uint32_t slot=0;
    uint32_t itrigger = -1;
    uint32_t last_itrigger = -2;
    uint32_t last_pulse_time_channel=0;
    uint32_t last_slot = -1;
    uint32_t last_channel = -1;    

    // Loop over data words
    for(; iptr<iend; iptr++){

        // Skip all non-data-type-defining words at this
        // level. When we do encounter one, the appropriate
        // case block below should handle parsing all of
        // the data continuation words and advance the iptr.
        if(((*iptr>>31) & 0x1) == 0)continue;

        uint32_t data_type = (*iptr>>27) & 0x0F;
        switch(data_type){
            case 0: // Block Header
                slot = (*iptr>>22) & 0x1F;
                if(VERBOSE>7) cout << "      FADC125 Block Header: slot="<<slot<<endl;
                break;
            case 1: // Block Trailer
				pe_iter = current_parsed_events.begin();
				pe = NULL;
				break;
            case 2: // Event Header
                //slot_event_header = (*iptr>>22) & 0x1F;
                itrigger = (*iptr>>0) & 0x3FFFFFF;
				pe = *pe_iter++;
                if(VERBOSE>7) cout << "      FADC125 Event Header: itrigger="<<itrigger<<" last_itrigger="<<last_itrigger<<", rocid="<<rocid<<", slot="<<slot <<endl;
				break;
            case 3: // Trigger Time
				{
					uint64_t t = ((*iptr)&0xFFFFFF)<<0;
					iptr++;
					if(((*iptr>>31) & 0x1) == 0){
						t += ((*iptr)&0xFFFFFF)<<24; // from word on the street: second trigger time word is optional!!??
					}else{
						iptr--;
					}
					if(VERBOSE>7) cout << "      FADC125 Trigger Time (t="<<t<<")"<<endl;
					if(pe) pe->NEW_Df125TriggerTime(rocid, slot, itrigger, t);
				}
                break;
            case 4: // Window Raw Data
					// iptr passed by reference and so will be updated automatically
					if(VERBOSE>7) cout << "      FADC125 Window Raw Data"<<endl;
					MakeDf125WindowRawData(pe, rocid, slot, itrigger, iptr);
					break;

            case 5: // CDC pulse data (new)  (GlueX-doc-2274-v8)
				{
					// Word 1:
					uint32_t word1          = *iptr;
					uint32_t channel        = (*iptr>>20) & 0x7F;
					uint32_t pulse_number   = (*iptr>>15) & 0x1F;
					uint32_t pulse_time     = (*iptr>>4 ) & 0x7FF;
					uint32_t quality_factor = (*iptr>>3 ) & 0x1; //time QF bit
					uint32_t overflow_count = (*iptr>>0 ) & 0x7;
					if(VERBOSE>7){
						cout << "      FADC125 CDC Pulse Data word1: " << hex << (*iptr) << dec << endl;
						cout << "      FADC125 CDC Pulse Data (chan="<<channel<<" pulse="<<pulse_number<<" time="<<pulse_time<<" QF="<<quality_factor<<" OC="<<overflow_count<<")"<<endl;
					}

					// Word 2:
					++iptr;
					if(iptr>=iend){
						jerr << " Truncated f125 CDC hit (block ends before continuation word!)" << endl;
						continue;
					}
					if( ((*iptr>>31) & 0x1) != 0 ){
						jerr << " Truncated f125 CDC hit (missing continuation word!)" << endl;
						continue;
					}
					uint32_t word2      = *iptr;
					uint32_t pedestal   = (*iptr>>23) & 0xFF;
					uint32_t sum        = (*iptr>>9 ) & 0x3FFF;
					uint32_t pulse_peak = (*iptr>>0 ) & 0x1FF;
					if(VERBOSE>7){
						cout << "      FADC125 CDC Pulse Data word2: " << hex << (*iptr) << dec << endl;
						cout << "      FADC125 CDC Pulse Data (pedestal="<<pedestal<<" sum="<<sum<<" peak="<<pulse_peak<<")"<<endl;
					}

					// Create hit objects
					uint32_t nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
					uint32_t nsamples_pedestal = 1;  // The firmware pedestal divided by 2^PBIT where PBIT is a config. parameter

					if( pe ) {
						pe->NEW_Df125CDCPulse(rocid, slot, channel, itrigger
									, pulse_number        // NPK
									, pulse_time          // le_time
									, quality_factor      // time_quality_bit
									, overflow_count      // overflow_count
									, pedestal            // pedestal
									, sum                 // integral
									, pulse_peak          // first_max_amp
									, word1               // word1
									, word2               // word2
									, nsamples_pedestal   // nsamples_pedestal
									, nsamples_integral   // nsamples_integral
									, false);             // emulated
					}
				}
                break;

            case 6: // FDC pulse data-integral (new)  (GlueX-doc-2274-v8)
				{
					// Word 1:
					uint32_t word1          = *iptr;
					uint32_t channel        = (*iptr>>20) & 0x7F;
					uint32_t pulse_number   = (*iptr>>15) & 0x1F;
					uint32_t pulse_time     = (*iptr>>4 ) & 0x7FF;
					uint32_t quality_factor = (*iptr>>3 ) & 0x1; //time QF bit
					uint32_t overflow_count = (*iptr>>0 ) & 0x7;
					if(VERBOSE>7){
						cout << "      FADC125 FDC Pulse Data(integral) word1: " << hex << (*iptr) << dec << endl;
						cout << "      FADC125 FDC Pulse Data (chan="<<channel<<" pulse="<<pulse_number<<" time="<<pulse_time<<" QF="<<quality_factor<<" OC="<<overflow_count<<")"<<endl;
					}

					// Word 2:
					++iptr;
					if(iptr>=iend){
						jerr << " Truncated f125 FDC hit (block ends before continuation word!)" << endl;
						continue;
					}
					if( ((*iptr>>31) & 0x1) != 0 ){
						jerr << " Truncated f125 FDC hit (missing continuation word!)" << endl;
						continue;
					}
					uint32_t word2      = *iptr;
					uint32_t pulse_peak = 0;
					uint32_t sum        = (*iptr>>19) & 0xFFF;
					uint32_t peak_time  = (*iptr>>11) & 0xFF;
					uint32_t pedestal   = (*iptr>>0 ) & 0x7FF;
					if(VERBOSE>7){
						cout << "      FADC125 FDC Pulse Data(integral) word2: " << hex << (*iptr) << dec << endl;
						cout << "      FADC125 FDC Pulse Data (integral="<<sum<<" time="<<peak_time<<" pedestal="<<pedestal<<")"<<endl;
					}

					// Create hit objects
					uint32_t nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
					uint32_t nsamples_pedestal = 1;  // The firmware pedestal divided by 2^PBIT where PBIT is a config. parameter

					if( pe ) {
						pe->NEW_Df125FDCPulse(rocid, slot, channel, itrigger
									, pulse_number        // NPK
									, pulse_time          // le_time
									, quality_factor      // time_quality_bit
									, overflow_count      // overflow_count
									, pedestal            // pedestal
									, sum                 // integral
									, pulse_peak          // peak_amp
									, peak_time           // peak_time
									, word1               // word1
									, word2               // word2
									, nsamples_pedestal   // nsamples_pedestal
									, nsamples_integral   // nsamples_integral
									, false);             // emulated
					}
				}
                break;

            case 7: // Pulse Integral
				{
					if(VERBOSE>7) cout << "      FADC125 Pulse Integral"<<endl;
					uint32_t channel = (*iptr>>20) & 0x7F;
					uint32_t sum = (*iptr>>0) & 0xFFFFF;
					uint32_t quality_factor = 0;
					uint32_t nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
					uint32_t nsamples_pedestal = 1;  // The firmware returns an already divided pedestal
					uint32_t pedestal = 0;  // This will be replaced by the one from Df250PulsePedestal in GetObjects
					uint32_t pulse_number = 0;
					if (last_slot == slot && last_channel == channel) pulse_number = 1;
					last_slot = slot;
					last_channel = channel;
					if( pe ) pe->NEW_Df125PulseIntegral(rocid, slot, channel, itrigger, pulse_number, quality_factor, sum, pedestal, nsamples_integral, nsamples_pedestal);
				}
                break;
            case 8: // Pulse Time
				{
					if(VERBOSE>7) cout << "      FADC125 Pulse Time"<<endl;
					uint32_t channel = (*iptr>>20) & 0x7F;
					uint32_t pulse_number = (*iptr>>18) & 0x03;
					uint32_t pulse_time = (*iptr>>0) & 0xFFFF;
					uint32_t quality_factor = 0;
					if( pe ) pe->NEW_Df125PulseTime(rocid, slot, channel, itrigger, pulse_number, quality_factor, pulse_time);
					last_pulse_time_channel = channel;
				}
                break;

            case 9: // FDC pulse data-peak (new)  (GlueX-doc-2274-v8)
				{
					// Word 1:
					uint32_t word1          = *iptr;
					uint32_t channel        = (*iptr>>20) & 0x7F;
					uint32_t pulse_number   = (*iptr>>15) & 0x1F;
					uint32_t pulse_time     = (*iptr>>4 ) & 0x7FF;
					uint32_t quality_factor = (*iptr>>3 ) & 0x1; //time QF bit
					uint32_t overflow_count = (*iptr>>0 ) & 0x7;
					if(VERBOSE>7){
						cout << "      FADC125 FDC Pulse Data(peak) word1: " << hex << (*iptr) << dec << endl;
						cout << "      FADC125 FDC Pulse Data (chan="<<channel<<" pulse="<<pulse_number<<" time="<<pulse_time<<" QF="<<quality_factor<<" OC="<<overflow_count<<")"<<endl;
					}

					// Word 2:
					++iptr;
					if(iptr>=iend){
						jerr << " Truncated f125 FDC hit (block ends before continuation word!)" << endl;
						continue;
					}
					if( ((*iptr>>31) & 0x1) != 0 ){
						jerr << " Truncated f125 FDC hit (missing continuation word!)" << endl;
						continue;
					}
					uint32_t word2      = *iptr;
					uint32_t pulse_peak = (*iptr>>19) & 0xFFF;
					uint32_t sum        = 0;
					uint32_t peak_time  = (*iptr>>11) & 0xFF;
					uint32_t pedestal   = (*iptr>>0 ) & 0x7FF;
					if(VERBOSE>7){
						cout << "      FADC125 FDC Pulse Data(peak) word2: " << hex << (*iptr) << dec << endl;
						cout << "      FADC125 FDC Pulse Data (integral="<<sum<<" time="<<peak_time<<" pedestal="<<pedestal<<")"<<endl;
					}

					// Create hit objects
					uint32_t nsamples_integral = 0;  // must be overwritten later in GetObjects with value from Df125Config value
					uint32_t nsamples_pedestal = 1;  // The firmware pedestal divided by 2^PBIT where PBIT is a config. parameter

					if( pe ) {
						pe->NEW_Df125FDCPulse(rocid, slot, channel, itrigger
									, pulse_number        // NPK
									, pulse_time          // le_time
									, quality_factor      // time_quality_bit
									, overflow_count      // overflow_count
									, pedestal            // pedestal
									, sum                 // integral
									, pulse_peak          // peak_amp
									, peak_time           // peak_time
									, word1               // word1
									, word2               // word2
									, nsamples_pedestal   // nsamples_pedestal
									, nsamples_integral   // nsamples_integral
									, false);             // emulated
					}
				}
                break;

            case 10: // Pulse Pedestal (consistent with Beni's hand-edited version of Cody's document)
				{
					if(VERBOSE>7) cout << "      FADC125 Pulse Pedestal"<<endl;
					//channel = (*iptr>>20) & 0x7F;
					uint32_t channel = last_pulse_time_channel; // not enough bits to hold channel number so rely on proximity to Pulse Time in data stream (see "FADC125 dataformat 250 modes.docx")
					uint32_t pulse_number = (*iptr>>21) & 0x03;
					uint32_t pedestal = (*iptr>>12) & 0x1FF;
					uint32_t pulse_peak = (*iptr>>0) & 0xFFF;
					uint32_t nsamples_pedestal = 1;  // The firmware returns an already divided pedestal
					if( pe ) pe->NEW_Df125PulsePedestal(rocid, slot, channel, itrigger, pulse_number, pedestal, pulse_peak, nsamples_pedestal);
				}
                break;

            case 13: // Event Trailer
            case 14: // Data not valid (empty module)
            case 15: // Filler (non-data) word
                if(VERBOSE>7) cout << "      FADC125 ignored data type: " << data_type <<endl;
                break;
        }
    }

    // Chop off filler words
    for(; iptr<iend; iptr++){
        if(((*iptr)&0xf8000000) != 0xf8000000) break;
    }
}

//----------------
// MakeDf125WindowRawData
//----------------
void DEVIOWorkerThread::MakeDf125WindowRawData(DParsedEvent *pe, uint32_t rocid, uint32_t slot, uint32_t itrigger, uint32_t* &iptr)
{
    uint32_t channel = (*iptr>>20) & 0x7F;
    uint32_t window_width = (*iptr>>0) & 0x0FFF;

    Df125WindowRawData *wrd = pe->NEW_Df125WindowRawData(rocid, slot, channel, itrigger);

    for(uint32_t isample=0; isample<window_width; isample +=2){

        // Advance to next word
        iptr++;

        // Make sure this is a data continuation word, if not, stop here
        if(((*iptr>>31) & 0x1) != 0x0)break;

        bool invalid_1 = (*iptr>>29) & 0x1;
        bool invalid_2 = (*iptr>>13) & 0x1;
        uint16_t sample_1 = 0;
        uint16_t sample_2 = 0;
        if(!invalid_1)sample_1 = (*iptr>>16) & 0x1FFF;
        if(!invalid_2)sample_2 = (*iptr>>0) & 0x1FFF;

        // Sample 1
        wrd->samples.push_back(sample_1);
        wrd->invalid_samples |= invalid_1;
        wrd->overflow |= (sample_1>>12) & 0x1;

        if((isample+2) == window_width && invalid_2)break; // skip last sample if flagged as invalid

        // Sample 2
        wrd->samples.push_back(sample_2);
        wrd->invalid_samples |= invalid_2;
        wrd->overflow |= (sample_2>>12) & 0x1;
    }
}

//----------------
// ParseF1TDCBank
//----------------
void DEVIOWorkerThread::ParseF1TDCBank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend)
{
	if(!PARSE_F1TDC){ iptr = &iptr[(*iptr) + 1]; return; }

	uint32_t *istart = iptr;

	auto pe_iter = current_parsed_events.begin();
	DParsedEvent *pe = NULL;

    uint32_t slot = 0;
	uint32_t modtype = 0;
    uint32_t itrigger = -1;
	uint32_t trig_time_f1header = 0;

	// Some early data had a marker word at just before the actual F1 data
	if(*iptr == 0xf1daffff) iptr++;

    // Loop over data words
    for(; iptr<iend; iptr++){

        // Skip all non-data-type-defining words at this
        // level. When we do encounter one, the appropriate
        // case block below should handle parsing all of
        // the data continuation words and advance the iptr.
        if(((*iptr>>31) & 0x1) == 0)continue;

 		uint32_t data_type = (*iptr>>27) & 0x0F;
        switch(data_type){
			case 0: // Block Header
				slot = (*iptr)>>22 & 0x001F;
				modtype = (*iptr)>>18 & 0x000F;  // should match a DModuleType::type_id_t
				if(VERBOSE>7) cout << "      F1 Block Header: slot=" << slot << " modtype=" << modtype << endl;
				break;
		
            case 1: // Block Trailer
				pe_iter = current_parsed_events.begin();
				pe = NULL;
				if(VERBOSE>7) cout << "      F1 Block Trailer" << endl;
                break;

			case 2: // Event Header
				{
					pe = *pe_iter++;
					itrigger = (*iptr)>>0  & 0x0003FFFFF;
					if(VERBOSE>7) {
						uint32_t slot_event_header  = (*iptr)>>22 & 0x00000001F;
						cout << "      F1 Event Header: slot=" << slot_event_header << " itrigger=" << itrigger << endl;
					}
				}
				break;

			case 3: // Trigger time
				{
					uint64_t t = ((*iptr)&0xFFFFFF)<<0;
					iptr++;
					if(((*iptr>>31) & 0x1) == 0){
						t += ((*iptr)&0xFFFFFF)<<24; // from word on the street: second trigger time word is optional!!??
					}else{
						iptr--;
					}
					if(VERBOSE>7) cout << "      F1TDC Trigger Time (t="<<t<<")"<<endl;
					if(pe) pe->NEW_DF1TDCTriggerTime(rocid, slot, itrigger, t);
				}
				break;
			
			case 8: // F1 Chip Header
				trig_time_f1header    = ((*iptr)>> 7) & 0x1FF;
				if(VERBOSE>7) {
					uint32_t chip_f1header         = ((*iptr)>> 3) & 0x07;
					uint32_t chan_on_chip_f1header = ((*iptr)>> 0) & 0x07;  // this is always 7 in real data!
					uint32_t itrigger_f1header     = ((*iptr)>>16) & 0x3F;
					cout << "      Found F1 header: chip=" << chip_f1header << " chan=" << chan_on_chip_f1header << " itrig=" << itrigger_f1header << " trig_time=" << trig_time_f1header << endl;
				}
				break;

			case 7: // F1 Data
				{
					uint32_t chip         = (*iptr>>19) & 0x07;
					uint32_t chan_on_chip = (*iptr>>16) & 0x07;
					uint32_t time         = (*iptr>> 0) & 0xFFFF;
					uint32_t channel      = F1TDC_channel(chip, chan_on_chip, modtype);
					if(VERBOSE>7) cout << "      Found F1 data  : chip=" << chip << " chan=" << chan_on_chip  << " time=" << time << endl;
					if(pe) pe->NEW_DF1TDCHit(rocid, slot, channel, itrigger, trig_time_f1header, time, *iptr, MODULE_TYPE(modtype));
				}
				break;

			case 15: // Filler word
				if(VERBOSE>7) cout << "      F1 filler word" << endl;
			case 14: // Data not valid (how to handle this?)
				break;

			default:
				cerr<<endl;
				cout.flush(); cerr.flush();
				_DBG_<<"Unknown data word in F1TDC block. Dumping for debugging:" << endl;
				for(const uint32_t *iiptr = istart; iiptr<iend; iiptr++){
					_DBG_<<"0x"<<hex<<*iiptr<<dec;
					if(iiptr == iptr)cerr<<"  <----";
					switch( (*iiptr) & 0xF8000000 ){
						case 0x80000000: cerr << "   F1 Block Header"; break;
						case 0x90000000: cerr << "   F1 Event Header"; break;
						case 0x98000000: cerr << "   F1 Trigger time"; break;
						case 0xC0000000: cerr << "   F1 Header"; break;
						case 0xB8000000: cerr << "   F1 Data"; break;
						case 0x88000000: cerr << "   F1 Block Trailer"; break;
						case 0xF8000000: cerr << "   Filler word"; break;
						case 0xF0000000: cerr << "   <module has no valid data>"; break;
						default: break;
					}
					cerr<<endl;
					if(iiptr > (iptr+4)) break;
				}
				throw JException("Unexpected word type in F1TDC block!", __FILE__, __LINE__);
				break;
		}
	}	

	// Skip filler words
	while(iptr<iend && (*iptr&0xF8000000)==0xF8000000)iptr++;
}

//----------------
// LinkAllAssociations
//----------------
void DEVIOWorkerThread::LinkAllAssociations(void)
{
	/// Find objects that should be linked as "associated objects"
	/// of one another and add to each other's list.
	for( auto pe : current_parsed_events){

// auto svDf250PulseIntegral = pe->vDf250PulseIntegral;
// auto svDf125PulseIntegral = pe->vDf125PulseIntegral;
// auto svDf125CDCPulse      = pe->vDf125CDCPulse;
// auto svDf125FDCPulse      = pe->vDf125FDCPulse;
// auto svDF1TDCHit          = pe->vDF1TDCHit;
// auto svDCAEN1290TDCHit    = pe->vDCAEN1290TDCHit;

		//----------------- Sort all associations
		// fADC250
		if(pe->vDf250Config.size()>1       ) sort(pe->vDf250Config.begin(),        pe->vDf250Config.end(),        SortByROCID<Df250Config>              );
		if(pe->vDf250TriggerTime.size()>1  ) sort(pe->vDf250TriggerTime.begin(),   pe->vDf250TriggerTime.end(),   SortByModule<Df250TriggerTime>        );
		if(pe->vDf250PulseIntegral.size()>1) sort(pe->vDf250PulseIntegral.begin(), pe->vDf250PulseIntegral.end(), SortByPulseNumber<Df250PulseIntegral> );
		if(pe->vDf250PulseTime.size()>1    ) sort(pe->vDf250PulseTime.begin(),     pe->vDf250PulseTime.end(),     SortByPulseNumber<Df250PulseTime>     );
		if(pe->vDf250PulsePedestal.size()>1) sort(pe->vDf250PulsePedestal.begin(), pe->vDf250PulsePedestal.end(), SortByPulseNumber<Df250PulsePedestal> );

		// fADC125
		if(pe->vDf125Config.size()>1       ) sort(pe->vDf125Config.begin(),        pe->vDf125Config.end(),        SortByROCID<Df125Config>              );
		if(pe->vDf125TriggerTime.size()>1  ) sort(pe->vDf125TriggerTime.begin(),   pe->vDf125TriggerTime.end(),   SortByModule<Df125TriggerTime>        );
		if(pe->vDf125PulseIntegral.size()>1) sort(pe->vDf125PulseIntegral.begin(), pe->vDf125PulseIntegral.end(), SortByPulseNumber<Df125PulseIntegral> );
		if(pe->vDf125CDCPulse.size()>1     ) sort(pe->vDf125CDCPulse.begin(),      pe->vDf125CDCPulse.end(),      SortByChannel<Df125CDCPulse>          );
		if(pe->vDf125FDCPulse.size()>1     ) sort(pe->vDf125FDCPulse.begin(),      pe->vDf125FDCPulse.end(),      SortByChannel<Df125FDCPulse>          );
		if(pe->vDf125PulseTime.size()>1    ) sort(pe->vDf125PulseTime.begin(),     pe->vDf125PulseTime.end(),     SortByPulseNumber<Df125PulseTime>     );
		if(pe->vDf125PulsePedestal.size()>1) sort(pe->vDf125PulsePedestal.begin(), pe->vDf125PulsePedestal.end(), SortByPulseNumber<Df125PulsePedestal> );

		// F1TDC
		if(pe->vDF1TDCConfig.size()>1      ) sort(pe->vDF1TDCConfig.begin(),       pe->vDF1TDCConfig.end(),       SortByROCID<DF1TDCConfig>             );
		if(pe->vDF1TDCTriggerTime.size()>1 ) sort(pe->vDF1TDCTriggerTime.begin(),  pe->vDF1TDCTriggerTime.end(),  SortByModule<DF1TDCTriggerTime>       );
		if(pe->vDF1TDCHit.size()>1         ) sort(pe->vDF1TDCHit.begin(),          pe->vDF1TDCHit.end(),          SortByModule<DF1TDCHit>               );

		// CAEN1290TDC
		if(pe->vDCAEN1290TDCConfig.size()>1) sort(pe->vDCAEN1290TDCConfig.begin(), pe->vDCAEN1290TDCConfig.end(), SortByROCID<DCAEN1290TDCConfig>       );
		if(pe->vDCAEN1290TDCHit.size()>1   ) sort(pe->vDCAEN1290TDCHit.begin(),    pe->vDCAEN1290TDCHit.end(),    SortByModule<DCAEN1290TDCHit>         );


		//----------------- Link all associations
	
		// Connect Df250Config objects
		LinkConfigSamplesCopy(pe->vDf250Config, pe->vDf250PulseIntegral);

		// Connect Df125Config objects
		LinkConfigSamplesCopy(pe->vDf125Config, pe->vDf125PulseIntegral);
		LinkConfigSamplesCopy(pe->vDf125Config, pe->vDf125CDCPulse);
		LinkConfigSamplesCopy(pe->vDf125Config, pe->vDf125FDCPulse);

		// Connect DF1TDCConfig objects
		LinkConfig(pe->vDF1TDCConfig, pe->vDF1TDCHit);

		// Connect DCAEN1290TDCConfig objects
		LinkConfig(pe->vDCAEN1290TDCConfig, pe->vDCAEN1290TDCHit);

		// Connect Df250TriggerTime objects
		LinkModule(pe->vDf250TriggerTime, pe->vDf250PulseIntegral);

		// Connect Df125TriggerTime objects
		LinkModule(pe->vDf125TriggerTime, pe->vDf125PulseIntegral);
		LinkModule(pe->vDf125TriggerTime, pe->vDf125CDCPulse);
		LinkModule(pe->vDf125TriggerTime, pe->vDf125FDCPulse);

		// Connect DF1TDCTriggerTime objects
		LinkModule(pe->vDF1TDCTriggerTime, pe->vDF1TDCHit);

		// Connect Df250 pulse objects
		LinkPulse(pe->vDf250PulseTime,     pe->vDf250PulseIntegral);
		LinkPulsePedCopy(pe->vDf250PulsePedestal, pe->vDf250PulseIntegral);

		// Connect Df125 pulse objects
		LinkPulse(pe->vDf125PulseTime,     pe->vDf125PulseIntegral);
		LinkPulsePedCopy(pe->vDf125PulsePedestal, pe->vDf125PulseIntegral);

		// Connect Df250 window raw data objects
		if(!pe->vDf250WindowRawData.empty()){
			LinkConfig(pe->vDf250Config, pe->vDf250WindowRawData);
			LinkModule(pe->vDf250TriggerTime, pe->vDf250WindowRawData);
			LinkChannel(pe->vDf250WindowRawData, pe->vDf250PulseIntegral);
			LinkChannel(pe->vDf250WindowRawData, pe->vDf250PulseTime);
			LinkChannel(pe->vDf250WindowRawData, pe->vDf250PulsePedestal);
		}

		// Connect Df125 window raw data objects
		if(!pe->vDf125WindowRawData.empty()){
			LinkConfig(pe->vDf125Config, pe->vDf125WindowRawData);
			LinkModule(pe->vDf125TriggerTime, pe->vDf125WindowRawData);
			LinkChannel(pe->vDf125WindowRawData, pe->vDf125PulseIntegral);
			LinkChannel(pe->vDf125WindowRawData, pe->vDf125PulseTime);
			LinkChannel(pe->vDf125WindowRawData, pe->vDf125PulsePedestal);
			LinkChannel(pe->vDf125WindowRawData, pe->vDf125CDCPulse);
			LinkChannel(pe->vDf125WindowRawData, pe->vDf125FDCPulse);
		}

// pe->vDf250PulseIntegral = svDf250PulseIntegral;
// pe->vDf125PulseIntegral = svDf125PulseIntegral;
// pe->vDf125CDCPulse      = svDf125CDCPulse;
// pe->vDf125FDCPulse      = svDf125FDCPulse;
// pe->vDF1TDCHit          = svDF1TDCHit;
// pe->vDCAEN1290TDCHit    = svDCAEN1290TDCHit;

	}

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


