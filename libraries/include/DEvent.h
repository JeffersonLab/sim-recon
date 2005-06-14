// $Id$
//
//    File: DEvent.h
// Created: Wed Jun  8 12:30:53 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DEvent_
#define _DEvent_

#include <cstdio>
#include <vector>
#include <string>
using namespace std;

#include "derror.h"
#include "DEventSource.h"

class DFactory_base;
class DEventLoop;

class DEvent{
	public:
		DEvent();
		virtual ~DEvent();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEvent";}
		template<class T> derror_t GetObjects(vector<const T*> &t, const char* tag, DFactory_base *factory=NULL);
		inline DEventSource* GetDEventSource(void){return source;}
		inline int GetEventNumber(void){return event_number;}
		inline int GetRunNumber(void){return run_number;}
		inline void* GetRef(void){return ref;}
		inline DEventLoop* GetDEventLoop(void){return loop;}
		inline void SetDEventSource(DEventSource *source){this->source=source;}
		inline void SetRunNumber(int run_number){this->run_number=run_number;}
		inline void SetEventNumber(int event_number){this->event_number=event_number;}
		inline void SetRef(void *ref){this->ref=ref;}
		inline void SetDEventLoop(DEventLoop *loop){this->loop=loop;}
		inline void FreeEvent(void){if(source)source->FreeEvent(ref);}
	
	private:
		DEventSource *source;
		int event_number;
		int run_number;
		void *ref;
		DEventLoop *loop;

};

//---------------------------------
// GetObjects
//---------------------------------
template<class T>
derror_t DEvent::GetObjects(vector<const T*> &t, const char* tag, DFactory_base *factory)
{
	/// Call the GetObjects() method of the associated source.
	
	// Make sure source is at least not NULL
	if(!source)throw EVENT_SOURCE_NOT_OPEN;
	
	// Get list of object pointers cast as void*. This is needed because
	// you can't have a virtual template and we need to call GetObjects
	// as a virtual method of DEventSource. 
	vector<void*> v;
	derror_t err = source->GetObjects(T::className(), v, tag, ref, factory);
	if(err)return err;
	
	// OK, must have found some objects (possibly even zero) in the source.
	// The pointers are cast back to type const T* here.
	for(unsigned int i=0; i<v.size(); i++)t.push_back((const T*)v[i]);

	return NOERROR;
}

#endif // _DEvent_

