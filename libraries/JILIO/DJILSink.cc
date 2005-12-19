// Author: David Lawrence  Dec. 9, 2005
//
//
// DJILSink.cc
//

#include <iostream>
using namespace std;

#ifdef JILIO

#include "JILStreamPBF.h"
#include "DJILSink.h"
#include "DFactory_base.h"
#include "DEvent.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t DJILSink::init(void)
{
	s= new JILStreamPBF(filename, "w");
	s->SetPointerTracking(JILStream::PTR_NONE);

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
derror_t DJILSink::evnt(DEventLoop *loop, int eventnumber)
{
	// Get list of all factories
	vector<DFactory_base*> factories = loop->GetFactories();

	// Lock out other threads and write out the event
	LockState();
	s->StartNamedWrite("Event");
	for(unsigned int i=0; i<factories.size(); i++){
		if(!factories[i]->TestFactoryFlag(DFactory_base::WRITE_TO_OUTPUT))continue;
		factories[i]->StreamToOutput(s);
	}
	(*s)<<JILStream::END_NAMED;
	UnlockState();
	
	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
derror_t DJILSink::fini(void)
{
	delete s;

	return NOERROR;
}

#else // JILIO
bool DJILSink_is_defined = false;   // avoids warnings about objects with no symbols!
#endif // JILIO

