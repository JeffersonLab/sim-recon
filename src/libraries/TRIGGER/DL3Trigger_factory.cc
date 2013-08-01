// $Id$
//
//    File: DL3Trigger_factory.cc
// Created: Wed Jul 31 14:34:24 EDT 2013
// Creator: davidl (on Darwin harriet.jlab.org 11.4.2 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DL3Trigger_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DL3Trigger_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DL3Trigger_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DL3Trigger_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DL3Trigger *myDL3Trigger = new DL3Trigger;
	// myDL3Trigger->x = x;
	// myDL3Trigger->y = y;
	// ...
	// _data.push_back(myDL3Trigger);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DL3Trigger_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DL3Trigger_factory::fini(void)
{
	return NOERROR;
}

