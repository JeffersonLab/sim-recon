// $Id$
//
//    File: DPSCTDCHit_factory.cc
// Created: Wed Oct 15 16:48:32 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DPSCTDCHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPSCTDCHit_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPSCTDCHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPSCTDCHit_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DPSCTDCHit *myDPSCTDCHit = new DPSCTDCHit;
	// myDPSCTDCHit->x = x;
	// myDPSCTDCHit->y = y;
	// ...
	// _data.push_back(myDPSCTDCHit);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPSCTDCHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPSCTDCHit_factory::fini(void)
{
	return NOERROR;
}

