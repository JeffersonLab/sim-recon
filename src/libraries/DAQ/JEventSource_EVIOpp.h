// $Id$
//
//    File: JEventSource_EVIOpp.h
// Created: Tue Mar 29 08:14:42 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _JEventSource_EVIOpp_
#define _JEventSource_EVIOpp_

#include <atomic>
#include <chrono>
#include <cinttypes>

#include <JANA/jerror.h>
#include <JANA/JApplication.h>
#include <JANA/JEventSource.h>
#include <JANA/JEvent.h>
#include <JANA/JFactory.h>
#include <JANA/JStreamLog.h>

#include <DAQ/HDEVIO.h>
#include <DAQ/HDET.h>
#include <DAQ/DEVIOWorkerThread.h>
#include <DAQ/DParsedEvent.h>
#include <DAQ/DBORptrs.h>
#include <DAQ/Df250EmulatorAlgorithm.h>
#include <DAQ/Df125EmulatorAlgorithm.h>

#include <DANA/DStatusBits.h>

/// How this Event Source Works
/// ===================================================================
///
/// Overview
/// -------------------
/// This event source represents a complete rewrite of the original
/// JEventSource_EVIO class. It is highly optimized and while it
/// gets rid of some of the complexities that developed in the
/// original over time, in some ways it replaced them with complexities
/// of its own.
///
/// In a nut shell, this source launches a minimum of 2 dedicated 
/// threads for reading in and parsing the events. These are in
/// addition to any threads the JANA framework creates. There is
/// a single "dispatcher" thread that reads in the events and
/// assigns them to one of a group of "worker" threads.
///
/// JANA itself has a similar structure in that there is a single
/// "EventBufferThread" that grabs fully parsed events from this
/// class via its GetEvent method, and then hands them to one of 
/// a group of "processing" threads where the actual analysis begins.
///
/// One benefit of this design is that it allows certain objects,
/// methods, or member data to be modified by only a single thread
/// throughout its life. This avoids the use of locks which introduce
/// inefficiency in a multi-threaded application. 
///
///
/// BORconfig objects
/// --------------------
/// Typically, BOR events only occur at the begining of a file and
/// a given job will see only one. However, we must support the
/// possibility of multiple BOR events in case we are reading from 
/// an ET system or file of merged events.
///
/// BOR objects are applicable to all events in the stream until
/// another BOR event is encountered. We therefore need every event
/// to get a copy of the BOR config object pointers. Thus, they must be 
/// kept in a place that all worker threads see. This is the borptrs_list
/// object in the JEventSource_EVIOpp class.
///
/// Here's the sequence for how the BOR event is handled:
/// 1. When a BOR event is parsed it creates a DBORptrs object which
///    holds a complete set of BOR objects. The pointer to this 
///    DBORptrs object is left in the DParsedEvent object which
///    would otherwise be NULL. 
///
/// 2. When JEventSource_EVIOpp::GetEvent() is called, it will see the
///    non-NULL borptrs pointer in DParsedEvent. At this point it will
///    copy it to the front of borptrs_list, effectively giving ownership
///    to JEventSource_EVIOpp. Note that aside from the destructor,
///    borptrs_list is only accessed here which is only called from the
///    event reader thread. This means we don't have to use a lock.
///
/// 3. If the borptrs pointer is NULL, then JEventSource_EVIOpp::GetEvent
///    will copy the front pointer from borptrs_list into the
///    DParsedEvent. We rely on the events being in istream order here
///    to only apply the BORconfig objects to events after the BOR
///    event itself. 
///
/// 4. When GetObjects is called, it calls DParsedEvent::CopyToFactories
///    where the BOR objects are copied into the appropriate factories.
///    It then calls JEventSource_EVIOpp::LinkBORassociations to add
///    the BORConfig objects as associated objects to various hit objects.
///
/// 5. Only when the JEventSource_EVIOpp object is destroyed are any
///    BORConfig objects deleted. We don't implement a mechansim to
///    keep track of which DBORptrs objects are still in use so we can't
///    delete them sooner. This shouldn't be a problem though since BOR
///    events are rare.
///

class JEventSource_EVIOpp: public jana::JEventSource{
	public:

		enum EVIOSourceType{
			kNoSource,
			kFileSource,
			kETSource
		};

		enum EmulationModeType{
			kEmulationNone,
			kEmulationAlways,
			kEmulationAuto
		};


		                    JEventSource_EVIOpp(const char* source_name);
		           virtual ~JEventSource_EVIOpp();
		virtual const char* className(void){return static_className();}
		 static const char* static_className(void){return "JEventSource_EVIOpp";}
		
		               void Dispatcher(void);
		
		           jerror_t GetEvent(jana::JEvent &event);
		               void FreeEvent(jana::JEvent &event);
		           jerror_t GetObjects(jana::JEvent &event, jana::JFactory_base *factory);

		               void LinkBORassociations(DParsedEvent *pe);
		           uint64_t SearchFileForRunNumber(void);
		               void EmulateDf250Firmware(DParsedEvent *pe);
		               void EmulateDf125Firmware(DParsedEvent *pe);
		               void AddToCallStack(DParsedEvent *pe, JEventLoop *loop);
		               void AddSourceObjectsToCallStack(JEventLoop *loop, string className);
		               void AddEmulatedObjectsToCallStack(JEventLoop *loop, string caller, string callee);
		               void AddROCIDtoParseList(uint32_t rocid){ ROCIDS_TO_PARSE.insert(rocid); }
		      set<uint32_t> GetROCIDParseList(uint32_t rocid){ return ROCIDS_TO_PARSE; }

		
		bool DONE;
		std::chrono::high_resolution_clock::time_point tstart;
		std::chrono::high_resolution_clock::time_point tend;

		uint32_t MAX_PARSED_EVENTS;
		mutex PARSED_EVENTS_MUTEX;
		condition_variable PARSED_EVENTS_CV;
		list<DParsedEvent*> parsed_events;

		std::atomic<uint_fast64_t> NEVENTS_PROCESSED;
		std::atomic<uint_fast64_t> NDISPATCHER_STALLED;
		std::atomic<uint_fast64_t> NPARSER_STALLED;
		std::atomic<uint_fast64_t> NEVENTBUFF_STALLED;
		
		uint64_t MAX_EVENT_RECYCLES;
		uint64_t MAX_OBJECT_RECYCLES;

		EVIOSourceType source_type;
		HDEVIO *hdevio;
		HDET   *hdet;
		bool et_quit_next_timeout;

		vector<DEVIOWorkerThread*> worker_threads;
		thread *dispatcher_thread;

		JStreamLog evioout;
		
		uint32_t F250_EMULATION_MODE; // (EmulationModeType)
		uint32_t F125_EMULATION_MODE; // (EmulationModeType)
		uint32_t F250_EMULATION_VERSION;
		Df250EmulatorAlgorithm *f250Emulator;
		Df125EmulatorAlgorithm *f125Emulator;
		
		bool RECORD_CALL_STACK;
		set<uint32_t> ROCIDS_TO_PARSE;

		list<DBORptrs*> borptrs_list;

		bool     PARSE;
		bool     PARSE_F250;
		bool     PARSE_F125;
		bool     PARSE_F1TDC;
		bool     PARSE_CAEN1290TDC;
		bool     PARSE_CONFIG;
		bool     PARSE_BOR;
		bool     PARSE_EPICS;
		bool     PARSE_EVENTTAG;
		bool     PARSE_TRIGGER;
		bool     APPLY_TRANSLATION_TABLE;
		int      ET_STATION_NEVENTS;
		bool     ET_STATION_CREATE_BLOCKING;
		bool     LOOP_FOREVER;
		uint32_t USER_RUN_NUMBER;
		int      VERBOSE;
		int      VERBOSE_ET;
		float    TIMEOUT;
		uint32_t NTHREADS;
		bool     PRINT_STATS;
		bool     SWAP;
		bool     LINK;
		bool     LINK_TRIGGERTIME;
		bool     LINK_BORCONFIG;
		bool     LINK_CONFIG;
		bool     IGNORE_EMPTY_BOR;
		bool     TREAT_TRUNCATED_AS_ERROR;
		
		uint32_t jobtype;
};

#endif // _JEventSourceGenerator_EVIOpp_

