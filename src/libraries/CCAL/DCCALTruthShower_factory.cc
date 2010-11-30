// $Id$
//
//    File: DCCALTruthShower_factory.cc
// Created: Tue Nov 30 15:02:26 EST 2010
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DCCALTruthShower_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DCCALTruthShower_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCCALTruthShower_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCCALTruthShower_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DCCALTruthShower *myDCCALTruthShower = new DCCALTruthShower;
	// myDCCALTruthShower->x = x;
	// myDCCALTruthShower->y = y;
	// ...
	// _data.push_back(myDCCALTruthShower);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DCCALTruthShower_factory::erun(void)
{
	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DCCALTruthShower_factory::fini(void)
{
	return NOERROR;
}

