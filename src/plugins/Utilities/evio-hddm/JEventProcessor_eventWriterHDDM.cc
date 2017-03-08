// $Id$
//
//    File: JEventProcessor_eventWriterHDDM.cc
// Created: Mon Mar  6 14:11:50 EST 2017
// Creator: tbritton (on Linux halld03.jlab.org 3.10.0-514.6.1.el7.x86_64 x86_64)
//

#include "JEventProcessor_eventWriterHDDM.h"
using namespace jana;


// Routine used to create our JEventProcessor
#include <JANA/JApplication.h>
#include <JANA/JFactory.h>
extern "C"{
void InitPlugin(JApplication *app){
	InitJANAPlugin(app);
	app->AddProcessor(new JEventProcessor_eventWriterHDDM());
}
} // "C"


//------------------
// JEventProcessor_eventWriterHDDM (Constructor)
//------------------
JEventProcessor_eventWriterHDDM::JEventProcessor_eventWriterHDDM()
{

}

//------------------
// ~JEventProcessor_eventWriterHDDM (Destructor)
//------------------
JEventProcessor_eventWriterHDDM::~JEventProcessor_eventWriterHDDM()
{

}

//------------------
// init
//------------------
jerror_t JEventProcessor_eventWriterHDDM::init(void)
{
	// This is called once at program startup. 

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t JEventProcessor_eventWriterHDDM::brun(JEventLoop *eventLoop, int32_t runnumber)
{
	// This is called whenever the run number changes
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t JEventProcessor_eventWriterHDDM::evnt(JEventLoop *loop, uint64_t eventnumber)
{
  // Write this event to the rest output stream.                                                                                                                                                             
    vector<const DEventWriterHDDM*> locEventWriterHDDMVector;
  loop->Get(locEventWriterHDDMVector);
  locEventWriterHDDMVector[0]->Write_HDDMEvent(loop, ""); 
  //  std::cout<<"done"<<std::endl;
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t JEventProcessor_eventWriterHDDM::erun(void)
{
	// This is called whenever the run number changes, before it is
	// changed to give you a chance to clean up before processing
	// events from the next run number.
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t JEventProcessor_eventWriterHDDM::fini(void)
{
	// Called before program exit after event processing is finished.
	return NOERROR;
}

