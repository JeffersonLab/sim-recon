// $Id$
//
//    File: DEventSink.h
// Created: Mon Dec 19 16:15:49 EST 2005
// Creator: davidl (on Linux phecda 2.6.9-11.ELsmp athlon)
//

#pragma interface

#ifndef _DEventSink_
#define _DEventSink_

#include <pthread.h>
#include <string>
#include <typeinfo>
using namespace std;

#include "derror.h"
#include "DFactory_base.h"
#include "DEventProcessor.h"
#include "DEventLoop.h"

class DEventSink:public DEventProcessor{
	public:
		DEventSink();
		virtual ~DEventSink(){}
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventSink";}
		
	protected:
		virtual derror_t init(void){return NOERROR;}				///< Called once at program start.
		derror_t brun(DEventLoop *eventLoop, int runnumber);	///< Called everytime a new run number is detected.
		virtual derror_t brun_sink(DEventLoop *loop, int runnumber)=0;
		virtual derror_t evnt(DEventLoop *loop, int eventnumber){return NOERROR;}	///< Called every event.
		virtual derror_t erun(void){return NOERROR;}				///< Called everytime run number changes, provided brun has been called.
		virtual derror_t fini(void){return NOERROR;}				///< Called after last event of last event source has been processed.
		void AddToWriteList(string name, string tag);
		void AddAllToWriteList(DEventLoop *loop);
		void RemoveFromWriteList(string name, string tag);
		void ClearWriteList(void);

		inline void LockSink(void){pthread_mutex_lock(&sink_mutex);}
		inline void UnlockSink(void){pthread_mutex_unlock(&sink_mutex);}
		inline bool IsInWriteList(const string &name, const string &tag){
			for(unsigned int i=0; i<factories_to_write.size(); i++){
				if(factories_to_write[i].name == name){
					if(factories_to_write[i].tag == tag)return true;
				}
			}
			return false;
		}
	
	private:
		typedef struct{
			string name;
			string tag;
		}factory_name_spec_t;
		vector<factory_name_spec_t> factories_to_write;
		
		bool initialized;
		pthread_mutex_t sink_mutex;
};

#endif // _DEventSink_

