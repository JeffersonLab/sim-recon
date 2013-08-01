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
	// Simple pass-through L3 trigger
	// algorithm = 0x1

	DL3Trigger *l3trig = new DL3Trigger(DL3Trigger::kKEEP_EVENT, 0x0L, 0x1);
	_data.push_back(l3trig);

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

