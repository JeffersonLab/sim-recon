// $Id$
//
//    File: DApplication.h
// Created: Wed Jun  8 12:00:20 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DApplication_
#define _DApplication_

#include <pthread.h>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

#include "derror.h"
#include "DParameter.h"


class DEventProcessor;
class DEventSource;
class DEventLoop;
class DEvent;
class DGeometry;

// These are for shared objects
typedef const char* GetDEventSourceType_t(void);
typedef DEventSource* MakeDEventSource_t(const char* name);
typedef void InitFactories_t(DEventLoop* eventLoop);

class DApplication{
	public:
		DApplication(int narg, char* argv[]);
		virtual ~DApplication();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DApplication";}
		
		derror_t NextEvent(DEvent &event);
		derror_t AddProcessor(DEventProcessor *processor);
		derror_t RemoveProcessor(DEventProcessor *processor);
		derror_t GetProcessors(vector<DEventProcessor*> &processors);
		derror_t AddDEventLoop(DEventLoop *loop);
		derror_t RemoveDEventLoop(DEventLoop *loop);
		derror_t GetDEventLoops(vector<DEventLoop*> &loops);
		DGeometry* GetGeometry(unsigned int run_number);
		derror_t Init(void);
		derror_t Run(DEventProcessor *proc=NULL, int Nthreads=0);
		derror_t Fini(void);
		void Pause(void);
		void Resume(void);
		void Quit(void);
		inline int GetNEvents(void){return NEvents;}
		inline float GetRate(void){return rate_instantaneous;}
		void PrintRate();
		void SetShowTicker(int what){show_ticker = what;}
		
	private:
	
		string Val2StringWithPrefix(float val);
		derror_t OpenNext(void);
		derror_t RegisterSharedObject(const char *soname);
		derror_t RegisterSharedObjectDirectory(const char *sodirname);

		vector<const char*> source_names;
		vector<DEventSource*> sources;
		DEventSource *current_source;
		pthread_mutex_t current_source_mutex;
	
		vector<DEventProcessor*> processors;
		vector<DEventLoop*> loops;
		vector<pthread_t> threads;
		pthread_mutex_t app_mutex;
		
		vector<DGeometry*> geometries;
		pthread_mutex_t geometry_mutex;
		
		typedef struct{
			const char* name;
			const char *soname;
			MakeDEventSource_t *MakeDEventSource;
		}EventSourceSharedObject_t;
		vector<EventSourceSharedObject_t> EventSourceSharedObjects;
		vector<InitFactories_t*> InitFactoriesProcs;

		int show_ticker;
		int NEvents;
		int last_NEvents;
		int avg_NEvents;
		double avg_time;
		double rate_instantaneous;
		double rate_average;
};


#endif // _DApplication_

