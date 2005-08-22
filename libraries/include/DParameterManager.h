// $Id$
//
//    File: DParameterManager.h
// Created: Tue Aug 16 14:30:24 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DParameterManager_
#define _DParameterManager_

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
using namespace std;

#include "DParameter.h"

#include "derror.h"

class DParameterManager{
	public:
		DParameterManager();
		virtual ~DParameterManager();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DParameterManager";}
	
		template<typename K, typename V> DParameter* SetDefaultParameter(K key, V& val);
		template<typename K, typename V> DParameter* SetParameter(K key, V val);
		template<typename K> DParameter* GetParameterNoLock(K key);
		template<typename K> DParameter* DParameterManager::GetParameter(K key);
		template<typename K, typename V> DParameter* DParameterManager::GetParameter(K key, V &val);
		void PrintParameters(void);
		void Dump(void);
		
	private:
		vector<DParameter*> parameters;
		pthread_mutex_t parameter_mutex;
		bool printParametersCalled;

};

// Global variable for accessing parameters (defined in DParameterManager.cc)
extern DParameterManager dparms;


//---------------------------------
// SetDefaultParameter
//---------------------------------
template<typename K, typename V>
DParameter* DParameterManager::SetDefaultParameter(K key, V &val)
{
	V my_val = val;
	
	DParameter *p = GetParameter(key,val);
	if(!p){
		p = SetParameter(key, val);
		p->isdefault = true;
	}else{
		// Warn user if two different default values are set
		if(p->isdefault){
			stringstream ss;
			ss<<val;
			if(ss.str() != p->GetValue()){
				cout<<" WARNING: Multiple calls to SetDefaultParameter with key=\""
				<<key<<"\" value= "<<val<<" and "<<my_val<<endl;
			}
		}
	}
	
	// Set the "hasdefault" flag so typos can be filtered later
	p->hasdefault = true;

	return p;
}

//---------------------------------
// SetParameter
//---------------------------------
template<typename K, typename V>
DParameter* DParameterManager::SetParameter(K key, V val)
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
	p->type = DParameter::DataType(val);
	
	// Tell the DParameterManager it needs to re-print if requested
	printParametersCalled = false;

	// release the parameters mutex
	pthread_mutex_unlock(&parameter_mutex);

	return p;
}

//---------------------------------
// GetParameterNoLock
//---------------------------------
template<typename K>
DParameter* DParameterManager::GetParameterNoLock(K key)
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
DParameter* DParameterManager::GetParameter(K key)
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
DParameter* DParameterManager::GetParameter(K key, V &val)
{
	DParameter *p = GetParameter(key);
	if(p){
		// use stringstream to convert string into V
		stringstream ss(p->GetValue());
		ss>>val;
	}
	return p;
}

#endif // _DParameterManager_

