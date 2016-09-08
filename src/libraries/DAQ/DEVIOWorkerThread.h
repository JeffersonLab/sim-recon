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

class JEventSource_EVIOpp;


class DEVIOWorkerThread{
	public:
			
		enum JOBTYPE{
			JOB_NONE       = 0x0,
			JOB_QUIT       = 0x1,
			JOB_SWAP       = 0x2,
			JOB_FULL_PARSE = 0x4,
			JOB_ASSOCIATE  = 0x8
		};

		DEVIOWorkerThread(
			JEventSource_EVIOpp  *event_source
	 		,list<DParsedEvent*> &parsed_events
	 		,uint32_t            &MAX_PARSED_EVENTS
	 		,mutex               &PARSED_EVENTS_MUTEX
	 		,condition_variable  &PARSED_EVENTS_CV );
		virtual ~DEVIOWorkerThread();

		// These are owned by JEventSource and
		// are set in the constructor
		JEventSource_EVIOpp *event_source;
		list<DParsedEvent*> &parsed_events;
		uint32_t            &MAX_PARSED_EVENTS;
		mutex               &PARSED_EVENTS_MUTEX;
		condition_variable  &PARSED_EVENTS_CV;
		
		// Pool of parsed events
		vector<DParsedEvent*> parsed_event_pool;
		
		// List of parsed events we are currently filling
		list<DParsedEvent*> current_parsed_events;
	
		int VERBOSE;
		uint64_t Nrecycled;     // Incremented in JEventSource_EVIOpp::Dispatcher()
		uint64_t MAX_EVENT_RECYCLES;
		uint64_t MAX_OBJECT_RECYCLES;
	
		atomic<bool> in_use;
		atomic<bool> done;
		JOBTYPE jobtype;
		uint64_t istreamorder;
		uint64_t run_number_seed;
		
		mutex mtx;
		condition_variable cv;
		thread thd;
		
		uint32_t buff_len;
		uint32_t *buff;
		streampos pos;

		bool  PARSE_F250;
		bool  PARSE_F125;
		bool  PARSE_F1TDC;
		bool  PARSE_CAEN1290TDC;
		bool  PARSE_CONFIG;
		bool  PARSE_BOR;
		bool  PARSE_EPICS;
		bool  PARSE_EVENTTAG;
		bool  PARSE_TRIGGER;
		
		bool  LINK_TRIGGERTIME;
		bool  LINK_CONFIG;
		
		void Run(void);
		void Finish(bool wait_to_complete=true);
		void Prune(void);
		void MakeEvents(void);
		void PublishEvents(void);
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
        void      ParseDVertexBank(uint32_t* &iptr, uint32_t *iend);

		void        ParseJLabModuleData(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void                ParseTIBank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void              ParseCAEN1190(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void   ParseModuleConfiguration(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void              Parsef250Bank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void     MakeDf250WindowRawData(DParsedEvent *pe, uint32_t rocid, uint32_t slot, uint32_t itrigger, uint32_t* &iptr);
		void              Parsef125Bank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);
		void     MakeDf125WindowRawData(DParsedEvent *pe, uint32_t rocid, uint32_t slot, uint32_t itrigger, uint32_t* &iptr);
		void             ParseF1TDCBank(uint32_t rocid, uint32_t* &iptr, uint32_t *iend);

		void LinkAllAssociations(void);

		inline uint32_t F1TDC_channel(uint32_t chip, uint32_t chan_on_chip, int modtype);


		void DumpBinary(const uint32_t *iptr, const uint32_t *iend, uint32_t MaxWords=0, const uint32_t *imark=NULL);

	protected:
	
	
	private:

};

//----------------
// F1TDC_channel
//----------------
inline uint32_t DEVIOWorkerThread::F1TDC_channel(uint32_t chip, uint32_t chan_on_chip, int modtype)
{
    /// Convert a F1TDC chip number and channel on the chip to the
    /// front panel channel number. This is based on "Input Channel Mapping"
    /// section at the very bottom of the document F1TDC_V2_V3_4_29_14.pdf

    uint32_t channel_map[8] = {0, 0, 1, 1, 2, 2, 3, 3};
    switch(modtype){
        case DModuleType::F1TDC32:
            return (4 * chip) + channel_map[ chan_on_chip&0x7 ];
        case DModuleType::F1TDC48:
            return (chip <<3) | chan_on_chip;
        default:
            _DBG_ << "Calling F1TDC_channel for module type: " << DModuleType::GetName((DModuleType::type_id_t)modtype) << endl;
            throw JException("F1TDC_channel called for non-F1TDC module type");
    }
    return 1000000; // (should never get here)
}


#endif // _DEVIOWorkerThread_

