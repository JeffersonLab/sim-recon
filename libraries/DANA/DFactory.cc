// $Id$

#include "DFactory.h"
#include "DEvent.h"

//-------------
// DFactory
//-------------
DFactory::DFactory(DEvent *my_devent, char *my_name, int rowsize)
{
	/// This is a base class that specific factories inherit from.
	/// my_devent will be kept and used to allow this factory to
	/// access other factories.

	event = my_devent;
	name = my_name;
	_data = new DContainer(NULL, rowsize, my_name);
}

//-------------
// ~DFactory
//-------------
DFactory::~DFactory()
{

}

//-------------
// Get
//-------------
DContainer* DFactory::Get()
{
	/// Return the list of values for the type of data this factory
	/// generates. This is function checks first if the data already
	/// exists and then calls the GetData method if it doesn't.
	
	// If evnt_called is set, then just return the _data pointer
	if(evnt_called)return _data;
	
	// Make sure we're initialized
	if(!init_called){
		init();
		init_called = 1;
	}
	
	// Call brun routine if run number has changed or it's not been called
	if(event->runnumber!=brun_runnumber){
		if(brun_called && !erun_called){
			erun();
			erun_called = 1;
		}
		brun_called = 0;
	}
	if(!brun_called){
		brun(event->runnumber);
		brun_called = 1;
		erun_called = 0;
		brun_runnumber = event->runnumber;
	}
	
	// Call evnt routine to generate data
	evnt(event->eventnumber);
	evnt_called = 1;
	
	return _data;
}

