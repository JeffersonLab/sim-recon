// Author: David Lawrence  June 25, 2004
//
//
// DEventProcessor methods
//

#include <iostream>
using namespace std;

#include "DEventProcessor.h"

//----------------
// Constructor
//----------------
DEventProcessor::DEventProcessor(void)
{
	/// Constructor for DEventProcessor object 
	
	event_loop = NULL; // This will be set by the event loop object (if used)
	init_called = 0;
	brun_called = 0;
	evnt_called = 0;
	erun_called = 0;
	fini_called = 0;
	brun_runnumber = -1; // ensure brun is called
	_data = NULL;
}

//----------------
// Destructor
//----------------
DEventProcessor::~DEventProcessor()
{
	/// Destructor for DEventProcessor object 
}

//----------------
// GetStatus
//----------------
int DEventProcessor::GetStatus(void)
{
	int status = 0;
	if(init_called)status |= 0x1<<0;
	if(brun_called)status |= 0x1<<1;
	if(evnt_called)status |= 0x1<<2;
	if(erun_called)status |= 0x1<<3;
	if(fini_called)status |= 0x1<<4;
	
	return status;
}

//----------------
// init
//----------------
derror_t DEventProcessor::init(void)
{
	return NOERROR;
}

//----------------
// brun
//----------------
derror_t DEventProcessor::brun(int runnumber)
{
	return NOERROR;
}

//----------------
// evnt
//----------------
derror_t DEventProcessor::evnt(int eventnumber)
{
	return NOERROR;
}

//----------------
// erun
//----------------
derror_t DEventProcessor::erun(void)
{
	return NOERROR;
}

//----------------
// fini
//----------------
derror_t DEventProcessor::fini(void)
{
	return NOERROR;
}


