// Author: David Lawrence  Dec. 9, 2005
//
//
// DEventSinkJIL.cc
//

#include <iostream>
using namespace std;

#ifdef JILIO

#include "JILStreamPBF.h"
#include "DEventSinkJIL.h"
#include "DFactory_base.h"
#include "DEvent.h"

//------------------------------------------------------------------
// init   -Open output file here (e.g. a ROOT file)
//------------------------------------------------------------------
derror_t DEventSinkJIL::init(void)
{
	s= new JILStreamPBF(filename, "w");
	s->SetPointerTracking(JILStream::PTR_NONE);

	return NOERROR;
}

//------------------------------------------------------------------
// brun_sink
//------------------------------------------------------------------
derror_t DEventSinkJIL::brun_sink(DEventLoop *loop, int runnumber)
{
	AddAllToWriteList(loop);
	//RemoveFromWriteList("DFDCHit", "");
	//RemoveFromWriteList("DCDCHit", "");

	return NOERROR;
}

//------------------------------------------------------------------
// evnt   -Fill histograms here
//------------------------------------------------------------------
derror_t DEventSinkJIL::evnt(DEventLoop *loop, int eventnumber)
{
	// Get list of all factories
	vector<DFactory_base*> factories = loop->GetFactories();

	// Lock out other threads and write out the event
	LockSink();
	s->StartNamedWrite("Event");
	for(unsigned int i=0; i<factories.size(); i++){
		if(!IsInWriteList(factories[i]->dataClassName(), factories[i]->Tag()))continue;
		if(factories[i]->GetNrows()==0)continue; // This actually invokes the factory
		factories[i]->StreamToOutput(s);
	}
	(*s)<<JILStream::END_NAMED;
	UnlockSink();
	
	return NOERROR;
}

//------------------------------------------------------------------
// fini   -Close output file here
//------------------------------------------------------------------
derror_t DEventSinkJIL::fini(void)
{
	delete s;

	return NOERROR;
}

#else // JILIO
bool DEventSinkJIL_is_defined = false;   // avoids warnings about objects with no symbols!
#endif // JILIO

