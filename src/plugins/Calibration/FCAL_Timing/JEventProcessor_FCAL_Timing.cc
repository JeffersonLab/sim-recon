// $Id$
//
//    File: JEventProcessor_FCAL_Timing.cc
// Created: Mon Jan  4 13:03:20 EST 2016
// Creator: mstaib (on Linux ifarm1401 2.6.32-431.el6.x86_64 x86_64)
//

#include "JEventProcessor_FCAL_Timing.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_FCAL_Timing());
}
} // "C"


//------------------
// JEventProcessor_FCAL_Timing (Constructor)
//------------------
JEventProcessor_FCAL_Timing::JEventProcessor_FCAL_Timing()
{

}

//------------------
// ~JEventProcessor_FCAL_Timing (Destructor)
//------------------
JEventProcessor_FCAL_Timing::~JEventProcessor_FCAL_Timing()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_FCAL_Timing::init(void)
{
	// This is called once at program startup. If you are creating
	// and filling historgrams in this plugin, you should lock the
	// ROOT mutex like this:
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();
	//

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_FCAL_Timing::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_FCAL_Timing::evnt(JEventLoop *loop, uint64_t eventnumber)
{
	// This is called for every event. Use of common resources like writing
	// to a file or filling a histogram should be mutex protected. Using
	// loop->Get(...) to get reconstructed objects (and thereby activating the
	// reconstruction algorithm) should be done outside of any mutex lock
	// since multiple threads may call this method at the same time.
	// Here's an example:
	//
	// vector<const MyDataClass*> mydataclasses;
	// loop->Get(mydataclasses);
	//
	// japp->RootWriteLock();
	//  ... fill historgrams or trees ...
	// japp->RootUnLock();


	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_FCAL_Timing::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_FCAL_Timing::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

