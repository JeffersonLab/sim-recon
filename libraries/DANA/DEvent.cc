// $Id$

#include "DEvent.h"

//-------------
// DEvent
//-------------
DEvent::DEvent()
{
	factories = new DContainer(NULL,sizeof(DFactory*), "factories");
	
	runnumber=0;
	eventnumber=0;
}

//-------------
// ~DEvent
//-------------
DEvent::~DEvent()
{

}

//-------------
// Get
//-------------
DContainer* DEvent::Get(char *data_name)
{
	/// Currently, this allows infinite recursion should two
	/// factories (eventually) reference one another. A mechanism
	/// needs to be written here to catch that error and print
	/// an error message.

	// Search for specified factory and return pointer to it's data container
	DFactory **factory = (DFactory**)*(factories->container_ptr);
	for(int i=0; i<factories->nrows; i++, factory++){
		if(!strcmp((*factory)->name, data_name)){
			(*factory)->hddm_s = hddm_s;
			return (*factory)->Get();
		}
	}

	// No factory found. Return NULL
	return NULL;
}

//-------------
// AddFactory
//-------------
derror_t DEvent::AddFactory(DFactory* factory)
{
	*(DFactory**)factories->Add() = factory;

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
	cout<<"Registered factories: ("<<factories->nrows<<" total)"<<endl;
	cout<<"----------------------"<<endl;
	cout<<"Name:            nrows: rowsize: maxrows:"<<endl;
	cout<<endl;

	DFactory **factory = (DFactory**)*(factories->container_ptr);
	for(int i=0; i<factories->nrows; i++, factory++){
		
		if(sparsify)
			if((*factory)->GetNrows()<1)continue;
		
		// To make things look pretty, copy all values into the buffer "str"
		char str[80];
		memset(str,' ',80);
		str[79] = 0;
		strncpy(str,(*factory)->name, strlen((*factory)->name));

		char num[32];
		sprintf(num, "%d", (*factory)->GetNrows());
		strncpy(&str[22-strlen(num)], num, strlen(num));

		sprintf(num, "%d", (*factory)->GetRowsize());
		strncpy(&str[31-strlen(num)], num, strlen(num));

		sprintf(num, "%d", (*factory)->GetMaxrows());
		strncpy(&str[40-strlen(num)], num, strlen(num));

		cout<<str<<endl;
		//cout<<"0x"<<hex<<(unsigned long)*factory<<dec<<endl;
	}
	
	cout<<"----------------------"<<endl;
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

	DFactory **factory = (DFactory**)*(factories->container_ptr);
	for(int i=0; i<factories->nrows; i++, factory++){
		int flags = (*factory)->GetContainerFlags();
		if(flags & DContainer::PERSISTANT)continue;
		(*factory)->Clear_evnt_called();
		(*factory)->ResetNrows();
	}

	return NOERROR;
}

//-------------
// Print
//-------------
derror_t DEvent::Print(char *data_name)
{
	/// Print the current data from this factory.
	///
	/// This method
	/// first invokes the Get() method for the factory to make
	/// sure the data has been generated for this event. Then,
	/// it just calls the factory's Print() method to do the
	/// actual printing to the screen.

	// Search for specified factory and return pointer to it's data container
	DFactory **factory = (DFactory**)*(factories->container_ptr);
	for(int i=0; i<factories->nrows; i++, factory++){
		if(!strcmp((*factory)->name, data_name)){
			(*factory)->hddm_s = hddm_s;
			(*factory)->Print();
			return NOERROR;
		}
	}

	// No factory found.
	return NOERROR;
}

//-------------
// GetFactoryNames
//-------------
DContainer* DEvent::GetFactoryNames(void)
{
	/// Return a DContainer whose members are char*s pointing to
	/// the names of the currently registered factories. The 
	/// caller must delete the DContainer object (but not the
	/// char*s since they are still in use by the factories
	/// themselves.
	DContainer *names = new DContainer(NULL, sizeof(char*), "factoryNames");
	DFactory **factory = (DFactory**)*(factories->container_ptr);
	for(int i=0; i<factories->nrows; i++, factory++){
		*(char**)names->Add() = (*factory)->name;
	}	
	
	return names;
}



