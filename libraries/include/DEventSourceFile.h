// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceFile
//
/// Implements DEventSource for HDDM files

#ifndef _DEVENT_SOURCEHDDM_H_
#define _DEVENT_SOURCEHDDM_H_

#include <vector>
#include <string>
using namespace std;

#include "DEventSource.h"
#include "derror.h"
#include "hddm_s.h"

class DEventSourceFile:public DEventSource
{
	public:
		DEventSourceFile(int narg, char *argv[]);
		~DEventSourceFile();
		
		derror_t Open(char *source);	///< Open the first (if more than one) event source 
		derror_t Close(void);			///< Close the current event source.
		derror_t GetEvent(void);		///< Copy event from ET system into buffer.
		derror_t GotoEvent(int eventno);				///< Jump to event given by eventno
		derror_t SetEventBufferSize(int Nevents);	///< Set number of event references to keep
		derror_t DEventSourceFile::RecordEventNumber(int eventno);	///< Current event's position cooresponds to this event number
		
		s_iostream_t *fin;
		
	private:
		typedef struct{
			int event;
			s_HDDM_t *hddm_s;
		}event_buffer_t;
		vector<event_buffer_t> event_buff;
		int max_event_buff;
		int current_event_index;
};

#endif //_DEVENT_SOURCEHDDM_H_
