// Author: David Lawrence  June 24, 2004
//
//
// DEventSourceHDDM
//
/// Implements DEventSource for HDDM files

#ifndef _DEVENT_SOURCEHDDM_H_
#define _DEVENT_SOURCEHDDM_H_

#include <vector>
#include <string>
using namespace std;

#include <pthread.h>

#include "DEventSource.h"
#include "derror.h"
#include "hddm_s.h"

class DEventSourceHDDM:public DEventSource
{
	public:
		DEventSourceHDDM(const char* source_name);
		virtual ~DEventSourceHDDM();		
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventSourceHDDM";}

		derror_t Goto(int runno, int eventno);
		derror_t GetObjects(const char *name, vector<void*> &v, const char* tag, void* ref, DFactory_base *factory);
		derror_t GetEvent(DEvent &event);
		void FreeEvent(void *ref);
		
		s_iostream_t *fin;
		s_HDDM_t *hddm_s;
		bool flush_on_free;
		
	private:
		typedef struct{
			int run_number;
			int event_number;
			s_HDDM_t *hddm_s;
		}event_buffer_t;
		vector<event_buffer_t> event_buff;
};

#endif //_DEVENT_SOURCEHDDM_H_
