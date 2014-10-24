// $Id$
//
//    File: DPSHit_factory.cc
// Created: Wed Oct 15 16:45:01 EDT 2014
// Creator: staylor (on Linux gluon05.jlab.org 2.6.32-358.18.1.el6.x86_64 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DPSHit_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPSHit_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPSHit_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPSHit_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DPSHit *myDPSHit = new DPSHit;
	// myDPSHit->x = x;
	// myDPSHit->y = y;
	// ...
	// _data.push_back(myDPSHit);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPSHit_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPSHit_factory::fini(void)
{
	return NOERROR;
}

