// $Id$
//
//    File: DEventLoop.h
// Created: Wed Jun  8 12:30:51 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#ifndef _DEventLoop_
#define _DEventLoop_

#include <vector>
#include <string>
using namespace std;

#include "DEvent.h"
#include "DFactory_base.h"
#include "derror.h"

template<class T> class DFactory;
class DApplication;
class DEventProcessor;

class DEventLoop{
	public:
		DEventLoop(DApplication *app);
		virtual ~DEventLoop();
		virtual const char* className(void){return static_className();}
		static const char* static_className(void){return "DEventLoop";}
		
		virtual derror_t AddFactory(DFactory_base* factory);
		derror_t RemoveFactory(DFactory_base* factory);
		DFactory_base* GetFactory(const string data_name, const char *tag="");
		vector<DFactory_base*> GetFactories(void){return factories;}
		vector<string> GetFactoryNames(void);
		derror_t ClearFactories(void);
		derror_t PrintFactories(int sparsify=0);
		derror_t Print(const string data_name);

		derror_t Loop(void);
		derror_t OneEvent(void);
		inline void Pause(void){pause = 1;}
		inline void Resume(void){pause = 0;}
		inline void Quit(void){quit = 1;}
		
		template<class T> DFactory<T>* Get(vector<const T*> &t, const char *tag="");
		template<class T> DFactory<T>* GetFromFactory(vector<const T*> &t, const char *tag="");
		template<class T> derror_t GetFromSource(vector<const T*> &t, const char *tag="", DFactory_base *factory=NULL);
		inline DEvent& GetDEvent(void){return event;}
		inline void SetDEvent(DEvent *event){this->event = *event;}
		inline void SetAutoFree(int auto_free){this->auto_free = auto_free;}
		
	private:
		DEvent event;
		vector<DFactory_base*> factories;
		vector<DEventProcessor*> processors;
		DApplication *app;
		double *heartbeat;
		int pause;
		int quit;
		int auto_free;

};


//-------------
// Get
//-------------
template<class T> 
DFactory<T>* DEventLoop::Get(vector<const T*> &t, const char *tag)
{
	/// Retrieve or generate the array of objects of
	/// type T for the curent event being processed
	/// by this thread.
	///
	/// By default, preference is given to reading the
	/// objects from the data source(e.g. file) before generating
	/// them in the factory. A flag exists in the factory
	/// however to change this so that the factory is
	/// given preference.
	///
	/// Note that regardless of the setting of this flag,
	/// the data are only either read in or generated once.
	/// Ownership of the objects will always be with the
	/// factory so subsequent calls will always return pointers to
	/// the same data (unless the factory doesn't exist in which
	/// case the objects will be read from the file multiple
	/// times).
	///
	/// If the factory is called on to generate the data,
	/// it is done by calling the factory's Get() method
	/// which, in turn, calls the evnt() method. 
	/// 
	/// First, we just call the GetFromFactory() method.
	/// It will make the initial decision as to whether
	/// it should look in the source first or not. If
	/// it returns NULL, then the factory couldn't be
	/// found so we automatically try the file.
	
	DFactory<T>* factory = GetFromFactory(t, tag);
	if(!factory)GetFromSource(t, tag, NULL);
	
	return (DFactory<T>*)factory;
}

//-------------
// GetFromFactory
//-------------
template<class T> 
DFactory<T>* DEventLoop::GetFromFactory(vector<const T*> &t, const char *tag)
{
	// We need to find the factory providing data type T. Since
	// the factory will return the result of a call to the same
	// static_className() method of class T, we only need to 
	// compare the pointers of the string and not the actual
	// string contents.

	const char* className = T::className();
	vector<DFactory_base*>::iterator iter=factories.begin();
	DFactory<T> *factory = NULL;
	for(; iter!=factories.end(); iter++){
		const char *factory_name = (*iter)->dataClassName();
		if(factory_name == className){
			if(!strcmp((*iter)->Tag(), tag)){
				factory = (DFactory<T>*)*iter;
				break;
			}
		}
	}
	
	// If factory not found, just return now
	if(!factory)return NULL;
	
	// OK, we found the factory. If the evnt() routine has already
	// been called, then just call the factory's Get() routine
	// to return a copy of the existing data
	if(factory->evnt_was_called()){
		vector<void*> &vt = factory->Get();
		t.clear();
		for(unsigned int i=0;i<vt.size();i++)t.push_back((const T*)vt[i]);
		return factory;
	}
	
	// Next option is to get the objects from the data source
	if(factory->CheckSourceFirst()){
		// If the object type/tag is found in the source, it
		// will return NOERROR, even if there are zero instances
		// of it. If it is not available in the source then it
		// will return OBJECT_NOT_AVAILABLE.
			
		derror_t err = GetFromSource(t, tag, factory);
		if(err == NOERROR){
			// copy data into factory and have it mark itself as having
			// valid data.
			factory->CopyExternal(t);
			return factory;
		}
	}
		
	// OK. It looks like we have to have the factory make this.
	// Get pointers to data from the factory. The pointers must be
	// passed to us as void* since DFactory_base doesn't 
	// know about the class type they hold. Here, we re-cast them 
	// back into the proper (const) type.
	vector<void*> &vt = factory->Get();
	t.clear();
	for(unsigned int i=0;i<vt.size();i++)t.push_back((const T*)vt[i]);
	
	return factory;
}

//-------------
// GetFromSource
//-------------
template<class T> 
derror_t DEventLoop::GetFromSource(vector<const T*> &t, const char *tag, DFactory_base *factory)
{
	/// This tries to get objects from the event source.
	/// factory may be NULL if no factory is associated
	/// with this class. However, HDDM has the routines
	/// to extract the class info from the s_HDDM_t structure
	/// kept in the classes themselves (Extract_HDDM method).
	/// Hence, calling this for a HDDM file with factory=NULL
	/// will result in empty object lists. Use the Get() method
	/// instead. That will call this, but with the appropriate
	/// value for factory set.
	
	return event.GetObjects(t, tag, factory);
}


#endif // _DEventLoop_

