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

#include "DObject.h"

#ifndef __CINT__  // disable this file from ROOT interpreter

#include <vector>
#include <string>
using namespace std;

#include "DFactory_base.h"
#include "DEvent.h"
#include "DEventLoop.h"

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
		virtual const char* Tag(void){return "";}
		inline int dataClassSize(void){return sizeof(T);}
		inline int GetEventCalled(void){return evnt_called;}
		inline int CheckSourceFirst(void){return !use_factory;}
		derror_t CopyExternal(vector<const T*> data);
		
	protected:
		vector<T*> _data;
		vector<void*>_vdata;
		int use_factory;
		
		derror_t Reset(void);
		derror_t HardReset(void);
		
		enum data_origin_t{
			ORIGIN_NONE,
			ORIGIN_FACTORY,
			ORIGIN_EXTERNAL
		};
		data_origin_t data_origin;
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
	/// Return a STL vector of pointers to the objects produced by the
	/// factory. The pointers must be type cast as void* because
	/// the DEventLoop object deals with DFactory_base objects which
	/// don't know anything specific about
	/// the objects produced by the factory. (Note that this method
	/// is accessed via the DFactory_base class's virtual Get() method).
	/// The templatized DEvent::Get() method will typecast the
	/// pointers back to the appropriate type for the user.
	///
	/// This method will check first to make sure this factory hasn't
	/// already been called for this event. If so, it just returns the
	/// existing _vdata vector right away.
	///
	/// The normal way to get here is through a call to
	/// DEvent::Get(). That method handles checking for the objects
	/// in the data source so we don't have to worry about it here.
	///
	/// This also uses a busy flag to ensure we're not called
	/// recursively. i.e. we call a factory who calls another
	/// factory who eventually calls us. An exception is thrown
	/// (type derror_t) with a value INFINITE_RECURSION if that
	/// situation is detected.
	
	// If evnt_called is set, then just return the _vdata pointer
	if(evnt_called)return _vdata;
	if(busy)throw(INFINITE_RECURSION);
	busy++;
	if(busy!=1)throw(INFINITE_RECURSION);  // Should we use a mutex here?
	
	// Grab the current event and run numbers
	int event_number = eventLoop->GetDEvent().GetEventNumber();
	int run_number = eventLoop->GetDEvent().GetRunNumber();
	
	// Make sure we're initialized
	if(!init_called){
		init();
		init_called = 1;
	}
	
	// Call brun routine if run number has changed or it's not been called
	if(run_number!=brun_runnumber){
		if(brun_called && !erun_called){
			erun();
			erun_called = 1;
		}
		brun_called = 0;
	}
	if(!brun_called){
		brun(eventLoop, run_number);
		brun_called = 1;
		erun_called = 0;
		brun_runnumber = run_number;
	}
	
	// Call evnt routine to generate data
	evnt(eventLoop, event_number);
	evnt_called = 1;
	busy=0;
	
	// Since DFactory_base doesn't know anything about type T
	// we have to copy the pointers into a vector of void*s.
	// They will be cast back into the proper pointer type in
	// DEvent::GetFromFactory().
	_vdata.clear();
	for(unsigned int i=0;i<_data.size();i++)_vdata.push_back((void*)_data[i]);
	
	return _vdata;
}

//-------------
// GetNrows
//-------------
template<class T>
const int DFactory<T>::GetNrows(void)
{
	// Call the factory's Get() method to make sure the data has been
	// obtained (either via data souce or factory) before returning
	// the number of rows. Return the size of _vdata Since that
	// should be valid in all cases.

	if(!evnt_called){
		if(!eventLoop)throw RESOURCE_UNAVAILABLE;
		
		vector<const T*> d; // dummy placeholder
		eventLoop->Get(d, Tag());
	}
	
	return _vdata.size();
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

	// don't reset the evnt_called flag for persistent data because this
	// will force evnt to be called next event therby regenerating
	// the data
	evnt_called = 0;
	
	return HardReset();
}

//-------------
// HardReset
//-------------
template<class T>
derror_t DFactory<T>::HardReset(void)
{
	
	/// Clear out the factories current contents.
	if(!(flags & NOT_OBJECT_OWNER)){
		for(unsigned int i=0;i<_data.size();i++){
			delete _data[i];
		}
	}
	_data.clear();
	_vdata.clear();

	evnt_called = 0;
	
	return NOERROR;
}

//-------------
// CopyExternal
//-------------
template<class T>
derror_t DFactory<T>::CopyExternal(vector<const T*> data)
{
	// Set flag so subsequent calls for this event will return this
	// data.
	evnt_called = 1;
		
	// Just copy into the _vdata vector since _data is not used outside
	// of the factory.
	_data.clear();
	_vdata.clear();
	for(unsigned int i=0;i<data.size();i++){
		_data.push_back((T*)data[i]);
		_vdata.push_back((void*)data[i]);
	}

	return NOERROR;
}


#endif // __CINT__

#endif // _DFACTORY_H_
