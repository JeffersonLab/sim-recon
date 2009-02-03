// $Id$
//
//    File: DPhoton_factory_HDParSim.cc
// Created: Tue Feb  3 11:29:30 EST 2009
// Creator: davidl (on Darwin harriet.jlab.org 9.6.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DPhoton_factory_HDParSim.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DPhoton_factory_HDParSim::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DPhoton_factory_HDParSim::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DPhoton_factory_HDParSim::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DPhoton *myDPhoton = new DPhoton;
	// myDPhoton->x = x;
	// myDPhoton->y = y;
	// ...
	// _data.push_back(myDPhoton);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DPhoton_factory_HDParSim::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DPhoton_factory_HDParSim::fini(void)
{
	return NOERROR;
}

