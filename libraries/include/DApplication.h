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
		
		template<typename K, typename V> DParameter* SetDefaultParameter(K key, V& val);
		template<typename K, typename V> DParameter* SetParameter(K key, V val);
		template<typename K> DParameter* GetParameterNoLock(K key);
		template<typename K> DParameter* DApplication::GetParameter(K key);
		template<typename K, typename V> DParameter* DApplication::GetParameter(K key, V &val);
		void PrintParameters(void);
		
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
		
		vector<DParameter*> parameters;
		pthread_mutex_t parameter_mutex;
		
		typedef struct{
			const char* name;
			const char *soname;
			MakeDEventSource_t *MakeDEventSource;
		}EventSourceSharedObject_t;
		vector<EventSourceSharedObject_t> EventSourceSharedObjects;
		vector<InitFactories_t*> InitFactoriesProcs;

		bool printDefaultParameters;
		int show_ticker;
		int NEvents;
		int last_NEvents;
		int avg_NEvents;
		double avg_time;
		double rate_instantaneous;
		double rate_average;
};


//---------------------------------
// SetDefaultParameter
//---------------------------------
template<typename K, typename V>
DParameter* DApplication::SetDefaultParameter(K key, V &val)
{
	V my_val = val;
	
	DParameter *p = GetParameter(key,val);
	if(!p){
		p = SetParameter(key, val);
		p->isdefault = true;
	}else{
		if(p->isdefault){
			cout<<" WARNING: Multiple calls to SetDefaultParameter with key=\""
				<<key<<"\" value= "<<val<<" and "<<my_val<<endl; 
		}
	}

	return p;
}

//---------------------------------
// SetParameter
//---------------------------------
template<typename K, typename V>
DParameter* DApplication::SetParameter(K key, V val)
{
	stringstream ss; // use a stringstream to convert type V into a string
	ss<<val;
	string skey(key); // key may be a const char* or a string
	string sval(ss.str());

	// block so one thread can't write while another reads
	pthread_mutex_lock(&parameter_mutex);

	DParameter *p = GetParameterNoLock(skey);
	if(!p){
		p = new DParameter(skey, sval);
		parameters.push_back(p);
	}else{
		p->SetValue(sval);
	}
	p->isdefault = false;

	// release the parameters mutex
	pthread_mutex_unlock(&parameter_mutex);

	return p;
}

//---------------------------------
// GetParameterNoLock
//---------------------------------
template<typename K>
DParameter* DApplication::GetParameterNoLock(K key)
{
	/// This is the thread un-safe routine for getting a DParameter*.
	/// Un-safe is mis-leading since this actually exists to provide
	/// thread safety. This is called by both SetParameter() and
	/// GetParameter(). Both of those routines lock the parameters
	/// mutex while calling this one. That way, we guarantee that
	/// the parameters list is not modified while it is being read.
	/// The only drawback is that reads are also serialized possibly
	/// losing a little in efficiency when running multi-threaded. 
	
	string skey(key);
	vector<DParameter*>::iterator iter = parameters.begin();
	for(; iter!= parameters.end(); iter++){
		if((*iter)->GetKey() == skey){
			return *iter;
		}
	}

	return NULL;
}

//---------------------------------
// GetParameter
//---------------------------------
template<typename K>
DParameter* DApplication::GetParameter(K key)
{
	/// Thread safe call to get a DParameter*
	
	// block so one thread can't write while another reads
	pthread_mutex_lock(&parameter_mutex);

	DParameter *p = GetParameterNoLock(key);

	// release the parameters mutex
	pthread_mutex_unlock(&parameter_mutex);
	
	return p;
}

//---------------------------------
// GetParameter
//---------------------------------
template<typename K, typename V>
DParameter* DApplication::GetParameter(K key, V &val)
{
	DParameter *p = GetParameter(key);
	if(p){
		// use stringstream to convert string into V
		stringstream ss(p->GetValue());
		ss>>val;
	}
	return p;
}

#endif // _DApplication_

