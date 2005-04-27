// $Id$

#include <iostream>
using namespace std;

#include "DEvent.h"
#include "DFactory_base.h"

//-------------
// DEvent
//-------------
DEvent::DEvent()
{
	factories.clear(); // probably don't need this
	
	_runnumber=0;
	_eventnumber=0;
	_hddm_s = NULL;
}

//-------------
// ~DEvent
//-------------
DEvent::~DEvent()
{
	// Delete all of the factories
	for(vector<DFactory_base*>::iterator i=factories.begin(); i!=factories.end(); i++){
		delete *i;
	}
	factories.clear();
}


//-------------
// GetFactory
//-------------
DFactory_base* DEvent::GetFactory(const string data_name)
{
	// Search for specified factory and return pointer to it
	vector<DFactory_base*>::iterator iter=factories.begin();
	for(; iter!=factories.end(); iter++){
		if(data_name == (*iter)->dataClassName()){
			return *iter;
			break;
		}
	}

	// No factory found. Return NULL
	return NULL;
}

//-------------
// AddFactory
//-------------
derror_t DEvent::AddFactory(DFactory_base* factory)
{
	factory->SetDEvent(this);
	factories.push_back(factory);

	return NOERROR;
}

//-------------
// PrintFactories
//-------------
derror_t DEvent::PrintFactories(int sparsify)
{
	/// Print a list of all registered factories to the screen
	/// along with a little info about each.

	cout<<endl;
	cout<<"Registered factories: ("<<factories.size()<<" total)"<<endl;
	cout<<endl;
	cout<<"Name:             nrows:"<<endl;
	cout<<"---------------- -------"<<endl;

	for(int i=0; i<factories.size(); i++){
		DFactory_base *factory = factories[i];
		
		if(_hddm_s)
			if(sparsify)
				if(factory->GetNrows()<1)continue;
		
		// To make things look pretty, copy all values into the buffer "str"
		string str(79,' ');
		string name = factory->dataClassName();
		str.replace(0, name.size(), name);

		if(_hddm_s){
			char num[32];
			sprintf(num, "%d", factory->GetNrows());
			str.replace(22-strlen(num), strlen(num), num);
		}

		cout<<str<<endl;
		//cout<<"0x"<<hex<<(unsigned long)*factory<<dec<<endl;
	}
	
	cout<<endl;

	return NOERROR;
}

//-------------
// ClearFactories
//-------------
derror_t DEvent::ClearFactories(void)
{
	/// Loop over all factories and clear their evnt_called flags
	/// (if apropriate). This is called from DEventLoop at the
	/// begining of a new event.

	for(int i=0; i<factories.size(); i++){
		factories[i]->Reset();
	}

	return NOERROR;
}

//-------------
// Print
//-------------
derror_t DEvent::Print(const string data_name)
{
	/// Dump the data to stdout for the specified factory
	///
	/// Find the factory corresponding to data_name and send
	/// the return value of its toString() method to stdout.

	// Search for specified factory and return pointer to it's data container
	DFactory_base *factory = GetFactory(data_name);
	if(factory){
		cout<<factory->toString();
		return NOERROR;
	}

	// No factory found.
	return NOERROR;
}

//-------------
// GetFactoryNames
//-------------
const vector<string> DEvent::GetFactoryNames(void)
{
	/// Return a vector<string> whose members are char*s pointing to
	/// the names of the currently registered factories. The 
	/// caller must delete the DContainer object (but not the
	/// char*s since they are still in use by the factories
	/// themselves.
	vector<string> names;
	vector<DFactory_base*>::iterator factory=factories.begin();
	for(; factory!=factories.end(); factory++){
		names.push_back((*factory)->dataClassName());
	}	
	
	return names;
}


//-------------
// Fini
//-------------
derror_t DEvent::Fini(void)
{
	/// Call fini() methods of all factories. Not all factories are
	/// accessed by every program so each one that implements a fini()
	/// method must not assume it's init(), brun(), evnt(), or erun()
	/// methods were called.
	vector<string> names;
	vector<DFactory_base*>::iterator factory=factories.begin();
	for(; factory!=factories.end(); factory++){
		(*factory)->fini();
	}

	return NOERROR;
}


