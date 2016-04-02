// $Id$
//
//    File: DEVIOWorkerThread.h
// Created: Mon Mar 28 07:40:07 EDT 2016
// Creator: davidl (on Darwin harriet.jlab.org 13.4.0 i386)
//

#ifndef _DEVIOWorkerThread_
#define _DEVIOWorkerThread_

#include <stdint.h>

#include <atomic>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <list>
using namespace std;

#include <JANA/jerror.h>
#include <DAQ/HDEVIO.h>
#include <DAQ/DParsedEvent.h>
#include <DAQ/DModuleType.h>

// Parsed events are placed here and access
// controlled by the mutex. The CV is used 
// to send notification whenever a new event
// is placed in the list.
//extern list<DParsedEvent*> parsed_events;
//extern mutex PARSED_EVENTS_MUTEX;
//extern condition_variable PARSED_EVENTS_CV;


class DEVIOWorkerThread{
	public:
			
		enum JOBTYPE{
			JOB_NONE       = 0x0,
			JOB_QUIT       = 0x1,
			JOB_SWAP       = 0x2,
			JOB_FULL_PARSE = 0x4
		};

		DEVIOWorkerThread(
	 		list<DParsedEvent*>  &parsed_events
	 		,uint32_t            &MAX_PARSED_EVENTS
	 		,mutex               &PARSED_EVENTS_MUTEX
	 		,condition_variable  &PARSED_EVENTS_CV );
		virtual ~DEVIOWorkerThread();

		// These are owned by JEventSource and
		// are set in the constructor
		list<DParsedEvent*> &parsed_events;
		uint32_t            &MAX_PARSED_EVENTS;
		mutex               &PARSED_EVENTS_MUTEX;
		condition_variable  &PARSED_EVENTS_CV;
		
		// reference to element in TLS_PARSED_EVENT that
		// is unique to this thread
		vector<DParsedEvent*> &parsed_event_pool;
		
		// List of parsed events we are currently filling
		list<DParsedEvent*> current_parsed_events;
	
		int VERBOSE;
	
		atomic<bool> in_use;
		atomic<bool> done;
		JOBTYPE jobtype;
		uint64_t istreamorder;
		
		mutex mtx;
		condition_variable cv;
		thread thd;
		
		uint32_t buff_len;
		uint32_t *buff;
		streampos pos;
				
		
		void Run(void);
		void Finish(bool wait_to_complete=true);
		void MakeEvents(void);
		void ParseBank(void);
		
		void     ParseEventTagBank(uint32_t* &iptr, uint32_t *iend);
		void        ParseEPICSbank(uint32_t* &iptr, uint32_t *iend);
		void          ParseBORbank(uint32_t* &iptr, uint32_t *iend);
		void     ParseTSscalerBank(uint32_t* &iptr, uint32_t *iend);
		void   Parsef250scalerBank(uint32_t* &iptr, uint32_t *iend);
		void     ParseControlEvent(uint32_t* &iptr, uint32_t *iend);
		void      ParsePhysicsBank(uint32_t* &iptr, uint32_t *iend);
		void ParseBuiltTriggerBank(uint32_t* &iptr, uint32_t *iend);
		void         ParseDataBank(uint32_t* &iptr, uint32_t *iend);

		void        ParseJLabModuleData(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void              ParseCAEN1190(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void   ParseModuleConfiguration(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void              Parsef250Bank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void              Parsef125Bank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void             ParseF1TDCBank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);

		void DumpBinary(const uint32_t *iptr, const uint32_t *iend, uint32_t MaxWords=0, const uint32_t *imark=NULL);
		
	protected:
	
	
	private:

};

#endif // _DEVIOWorkerThread_

