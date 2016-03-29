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
#include <DAQ/DEVIOWorkerThread.h>
#include <DAQ/DParsedEvent.h>

#include <DANA/DStatusBits.h>

class JEventSource_EVIOpp: public jana::JEventSource{
	public:
		JEventSource_EVIOpp(const char* source_name);
		virtual ~JEventSource_EVIOpp();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "JEventSource_EVIOpp";}
		
		void Dispatcher(void);
		
		jerror_t GetEvent(jana::JEvent &event);
		void FreeEvent(jana::JEvent &event);
		jerror_t GetObjects(jana::JEvent &event, jana::JFactory_base *factory);
		
		bool DONE;
		std::chrono::high_resolution_clock::time_point tstart;
		std::chrono::high_resolution_clock::time_point tend;

		uint32_t MAX_PARSED_EVENTS;
		mutex PARSED_EVENTS_MUTEX;
		condition_variable PARSED_EVENTS_CV;
		list<DParsedEvent*> parsed_events;

		std::atomic<uint_fast64_t> NEVENTS_PROCESSED;
		std::atomic<uint_fast64_t> NWAITS_FOR_THREAD;
		std::atomic<uint_fast64_t> NWAITS_FOR_PARSED_EVENT;
		
		HDEVIO *hdevio;
		vector<DEVIOWorkerThread*> worker_threads;
		thread *dispatcher_thread;


		JStreamLog evioout;

		bool  PARSE_EVIO_EVENTS;
		bool  PARSE_F250;
		bool  PARSE_F125;
		bool  PARSE_F1TDC;
		bool  PARSE_CAEN1290TDC;
		bool  PARSE_CONFIG;
		bool  PARSE_BOR;
		bool  PARSE_EPICS;
		bool  PARSE_EVENTTAG;
		bool  PARSE_TRIGGER;
		bool  MAKE_DOM_TREE;
		int   ET_STATION_NEVENTS;
		bool  ET_STATION_CREATE_BLOCKING;
		int   ET_DEBUG_WORDS_TO_DUMP;
		bool  LOOP_FOREVER;
		int   VERBOSE;
		float TIMEOUT;
		string MODTYPE_MAP_FILENAME;
		bool  ENABLE_DISENTANGLING;
		uint32_t N_WORKER_THREADS;
		bool  PRINT_STATS;
};

#endif // _JEventSourceGenerator_EVIOpp_

