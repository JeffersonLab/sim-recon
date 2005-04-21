// $Id$

/// Top-level Hall-D event object
///
/// This is actually a base class from which DEventLoop is derived.
/// Virtually all accesses to DEvent will be though its derived
/// form as a DEventLoop.


#ifndef _DEVENT_H_
#define _DEVENT_H_

#include <string>
#include <vector>
using namespace std;

#include "derror.h"
#include "hddm_s.h"

class DFactory_base;
template<class T> class DFactory;

class DEvent{
	public:
		DEvent();
		~DEvent();
		template<class T> DFactory<T>* Get(vector<T*> &t);
		DFactory_base* GetFactory(const string data_name);
		derror_t AddFactory(DFactory_base* factory);
		derror_t PrintFactories(void){return PrintFactories(0);};
		derror_t PrintFactories(int sparsify);
		derror_t ClearFactories(void);
		template<class T> const string toString(T*);
		const vector<string> GetFactoryNames(void);
		derror_t Print(const string data_name);
		
		inline const int runnumber(){return _runnumber;}
		inline const int eventnumber(){return _eventnumber;}
		inline s_HDDM_t* hddm_s(){return _hddm_s;}

	protected:
		s_HDDM_t *_hddm_s;
		int _runnumber;
		int _eventnumber;

	private:
		vector<DFactory_base*> factories;
};


//-------------
// Get
//-------------
template<class T> 
DFactory<T>* DEvent::Get(vector<T*> &t)
{
	/// Search through the list of factories and find the one
	/// who can supply the data type T. Call the factory's
	/// Get() method to get a reference to the data vector.
	/// We have to call the Get() method (which eventually
	/// calls the factory's evnt() method) through the
	/// DFactory_base virtual Get() method. The object pointers
	/// are  type cast as void* and placed in another vector before
	/// getting sent to us. We must cast it back to the proper
	/// type before passing it back to the user.

	// We need to find the factory providing data type T. Since
	// the factory will return the result of a call to the same
	// static className() method of class T, we only need to 
	// compare the pointers of the string and not the actual
	// string contents.
	const char* className = T::className();
	vector<DFactory_base*>::iterator iter=factories.begin();
	DFactory_base *factory = NULL;
	for(; iter!=factories.end(); iter++){
		if((*iter)->dataClassName() == className){
			factory = *iter;
			break;
		}
	}
	if(!factory)return (DFactory<T>*)NULL;
	
	// Get pointers to data from factory. The pointers must be
	// passed to us as void* since DFactory_base doesn't necessarily 
	// know about the class type they hold. Here, we re-cast them 
	// back into the proper type.
	vector<void*> vt = factory->Get();
	t.clear();
	for(int i=0;i<vt.size();i++)t.push_back((T*)vt[i]);

	return (DFactory<T>*)factory;
}

#endif // _DEVENT_H_
