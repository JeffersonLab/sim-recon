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
}

//----------------
// Destructor
//----------------
DEventProcessor::~DEventProcessor()
{
	/// Destructor for DEventProcessor object 
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


