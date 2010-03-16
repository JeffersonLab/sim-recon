// $Id$
//
//    File: DChargedTruthMatch_factory.cc
// Created: Sun Jan 31 08:45:38 EST 2010
// Creator: davidl (on Darwin Amelia.local 9.8.0 i386)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DChargedTruthMatch_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DChargedTruthMatch_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DChargedTruthMatch_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DChargedTruthMatch_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DChargedTruthMatch *myDChargedTruthMatch = new DChargedTruthMatch;
	// myDChargedTruthMatch->x = x;
	// myDChargedTruthMatch->y = y;
	// ...
	// _data.push_back(myDChargedTruthMatch);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DChargedTruthMatch_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DChargedTruthMatch_factory::fini(void)
{
	return NOERROR;
}

