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


#include <DAQ/HDEVIO.h>

#include <DParsedEvent.h>

// Parsed events are placed here and access
// controlled by the mutex. The CV is used 
// to send notification whenever a new event
// is placed in the list.
extern list<DParsedEvent*> parsed_events;
extern mutex PARSED_EVENTS_MUTEX;
extern condition_variable PARSED_EVENTS_CV;


class DEVIOWorkerThread{
	public:
			
		enum JOBTYPE{
			JOB_NONE       = 0x0,
			JOB_QUIT       = 0x1,
			JOB_SWAP       = 0x2,
			JOB_FULL_PARSE = 0x4
		};

		DEVIOWorkerThread();
		virtual ~DEVIOWorkerThread();
	
		atomic<bool> in_use;
		atomic<bool> done;
		JOBTYPE jobtype;
		uint64_t istreamorder;
		
		mutex mtx;
		condition_variable cv;
		thread thd;
		
		uint32_t buff_len;
		uint32_t *buff;
		
		void Run(void);
		void Finish(bool wait_to_complete=true);
		void ParseBank(void);

		
	protected:
	
	
	private:

};

#endif // _DEVIOWorkerThread_

