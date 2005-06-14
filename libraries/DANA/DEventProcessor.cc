// $Id$
//
//    File: DEventProcessor.cc
// Created: Wed Jun  8 12:31:12 EDT 2005
// Creator: davidl (on Darwin wire129.jlab.org 7.8.0 powerpc)
//

#include "DEventProcessor.h"

//---------------------------------
// DEventProcessor    (Constructor)
//---------------------------------
DEventProcessor::DEventProcessor(void)
{

	init_called = 0;
	brun_called = 0;
	evnt_called = 0;
	erun_called = 0;
	fini_called = 0;
	brun_runnumber = -1; // ensure brun is called
	pthread_mutex_init(&state_mutex, NULL);
	app = NULL;
}

//---------------------------------
// ~DEventProcessor    (Destructor)
//---------------------------------
DEventProcessor::~DEventProcessor()
{

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
derror_t DEventProcessor::brun(DEventLoop *loop, int runnumber)
{
	return NOERROR;
}

//----------------
// evnt
//----------------
derror_t DEventProcessor::evnt(DEventLoop *eventLoop, int eventnumber)
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
