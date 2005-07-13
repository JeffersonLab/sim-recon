// $Id$
//
//    File: DEventSource.h
// Created: Wed Jun  8 12:31:04 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DEventSource_
#define _DEventSource_

#include <vector>
#include <string>
using namespace std;

#include <pthread.h>

#include "derror.h"

class DFactory_base;
class DEvent;

class DEventSource{
	public:
		DEventSource(const char *source_name);
		virtual ~DEventSource();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventSource";}
		
		static const char* GuessSourceType(const char* source_name);
		static DEventSource* OpenSource(const char* name);
		virtual derror_t GetEvent(DEvent &event)=0;
		virtual void FreeEvent(void *ref){};
		virtual derror_t GetObjects(const char *name, vector<void*> &v, const char* tag, void* ref, DFactory_base *factory);

		inline void SetEventBufferSize(int size){event_buffer_size = size;}
		inline int GetEventBufferSize(void){return event_buffer_size;}
		inline const char* GetSourceName(void){return source_name;}

	protected:
		const char* source_name;
		int source_is_open;
		pthread_mutex_t read_mutex;
		int event_buffer_size;
		int Nevents_read;
	
		inline void LockRead(void){pthread_mutex_lock(&read_mutex);}
		inline void UnlockRead(void){pthread_mutex_unlock(&read_mutex);}
	
	private:

};

#endif // _DEventSource_

