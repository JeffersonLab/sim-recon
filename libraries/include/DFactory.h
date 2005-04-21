// $Id$

///
/// Hall-D Data Factory
///
/// All data (except that read from the data source) must
/// be derived from the data that was read from the source.
/// One DFactory object should exist for each type of data
/// that is "generated". For example: Clusters in the FCAL
/// are generated from Hits in the FCAL. Asking the factory
/// for a list of FCAL clusters will make the factory ask
/// for a list of FCAL Hits from the FCAL Hits factory which
/// may, in turn, ask for a list of raw data values from yet
/// another factory).
/// A list of DFactories is kept in the DEvent object so
/// all data may be accessed through it.

#ifndef _DFACTORY_H_
#define _DFACTORY_H_

#define HDCLASSDEF(T) \
	static const char* className(){return #T;}

#ifndef __CINT__  // disable this file from ROOT interpreter

#include <vector>
#include <string>
using namespace std;

#include "DFactory_base.h"
#include "DEvent.h"


//-----------------------
// class DFactory
//-----------------------
template<class T>
class DFactory:public DFactory_base{
	/// Templated base class which all factories are derived from.
	///
	/// At the user level, factory classes will be defined
	/// which inherit from this templated base class. By
	/// inheriting from a template, the derived factory class
	/// will automatically have a high degree of type safety
	/// since the "_data" vector will be specific to the type
	/// of objects the factory produces. This class (DFactory)
	/// inherits from DFactory_base so that all factories can
	/// be treated equally (polymorphism) by the DEvent object.
	/// Instantiating a DFactory<T> object itself would be pointless
	/// since they would use the default  brun(),evnt(),erun(),...
	/// methods from DEventProcessor which do nothing. Instead,
	/// A new class should be derived from this one which implements
	/// its own brun(),evnt(),erun(),... methods.
	
	public:
		DFactory();
		~DFactory();
		
		vector<void*>& Get();
		const int GetNrows(void);
		inline const char* className(void){return T::className();}
		inline const char* dataClassName(void){return className();}
		
	protected:
		vector<T*> _data;
		vector<void*>_vdata;
		
		derror_t Reset(void);
		derror_t HardReset(void);
		
};



//-------------
// DFactory
//-------------
template<class T>
DFactory<T>::DFactory()
{
	/// This is a base class that specific factories inherit from.
	/// my_devent will be kept and used to allow this factory to
	/// access other factories.

	// make sure vector is empty
	_data.clear(); // probably unnecessary

	// clear flags
	flags = DFACTORY_NULL;
	busy = 0;

	// Allow any factory to have its debug_level set via environment variable
	debug_level = 0;
	string envar = string() + "DEBUG_" + dataClassName();
	char *ptr = getenv(envar.c_str());
	if(ptr)debug_level = atoi(ptr);
}

//-------------
// ~DFactory
//-------------
template<class T>
DFactory<T>::~DFactory()
{
	/// Delete all objects in _data container.
	HardReset();
}

//-------------
// Get
//-------------
template<class T>
vector<void*>& DFactory<T>::Get()
{
	/// Return a reference to the vector object containing this
	/// factory's data, but type cast as a const vector<int*>&. The
	/// cast is done because the DEventLoop object has to
	/// keep a list of all factories and does so by treating them
	/// as DFactory_base objects. This Get() method is called
	/// by the DEventLoop calling the virtual DFactory_base method
	/// Get(). The DFactory_base class doesn't have type information
	/// in it. Thus, we must cast _data as a const vector<int*>&
	/// here and cast it back into the proper class type in
	/// DEventLoop::Get(). (make sense?).
	///
	/// This method will check first to make sure this factory hasn't
	/// already been called for this event. If so, it just returns the
	/// (type cast) _data pointer right away. Otherwise, it calls the
	/// factory's event() method to fill it first.
	///
	/// This also uses a busy flag to ensure were not called
	/// recursively. i.e. we call a factory who calls another
	/// factory who eventually calls us. An exception is thrown
	/// (type derror_t) with a value INFINITE_RECURSION if that
	/// situation is detected.
	
	// If evnt_called is set, then just return the _vdata pointer
	if(evnt_called)return _vdata;
	if(busy)throw(INFINITE_RECURSION);
	busy++;
	if(busy!=1)throw(INFINITE_RECURSION);  // Should we use a mutex here?
	
	// Copy DEvent's hddm_s pointer to our local one
	hddm_s = event->hddm_s();
	
	// Make sure we're initialized
	if(!init_called){
		init();
		init_called = 1;
	}
	
	// Call brun routine if run number has changed or it's not been called
	if(event->runnumber()!=brun_runnumber){
		if(brun_called && !erun_called){
			erun();
			erun_called = 1;
		}
		brun_called = 0;
	}
	if(!brun_called){
		brun(event->runnumber());
		brun_called = 1;
		erun_called = 0;
		brun_runnumber = event->runnumber();
	}
	
	// Call evnt routine to generate data
	evnt(event->eventnumber());
	evnt_called = 1;
	busy=0;
	
	// Since DFactory_base doesn't know anything about type T
	// we have to copy the pointers into a vector of void*s.
	// They will be cast back into the proper pointer type in
	// DEvent::Get().
	_vdata.clear();
	for(int i=0;i<_data.size();i++)_vdata.push_back((void*)_data[i]);
	
	return _vdata;
}

//-------------
// GetNrows
//-------------
template<class T>
const int DFactory<T>::GetNrows(void)
{
	// Call the Get() method to make sure the factory has produced
	// the data before returning the number of rows.
	Get();
	
	return _data.size();
}

//-------------
// Reset
//-------------
template<class T>
derror_t DFactory<T>::Reset(void)
{
	/// Clear out the factories current contents unless the
	/// PERSISTANT flag is set.
	if(flags & PERSISTANT)return NOERROR;
	
	return HardReset();
}

//-------------
// HardReset
//-------------
template<class T>
derror_t DFactory<T>::HardReset(void)
{
	/// Clear out the factories current contents.
	for(int i=0;i<_data.size();i++){
		delete _data[i];
	}
	_data.clear();
	_vdata.clear();

	evnt_called = 0;
	
	return NOERROR;
}

#endif // __CINT__

#endif // _DFACTORY_H_
