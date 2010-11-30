// $Id$
//
//    File: DCCALGeometry_factory.cc
// Created: Tue Nov 30 15:42:41 EST 2010
// Creator: davidl (on Linux ifarml6 2.6.18-128.el5 x86_64)
//


#include <iostream>
#include <iomanip>
using namespace std;

#include "DCCALGeometry_factory.h"
using namespace jana;

//------------------
// init
//------------------
jerror_t DCCALGeometry_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DCCALGeometry_factory::brun(jana::JEventLoop *eventLoop, int runnumber)
{
	flags = PERSISTANT;
	_data.push_back( new DCCALGeometry() );

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DCCALGeometry_factory::evnt(JEventLoop *loop, int eventnumber)
{

	// Code to generate factory data goes here. Add it like:
	//
	// DCCALGeometry *myDCCALGeometry = new DCCALGeometry;
	// myDCCALGeometry->x = x;
	// myDCCALGeometry->y = y;
	// ...
	// _data.push_back(myDCCALGeometry);
	//
	// Note that the objects you create here will be deleted later
	// by the system and the _data vector will be cleared automatically.

	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DCCALGeometry_factory::erun(void)
{
	for(unsigned int i=0; i<_data.size(); i++)delete _data[i];
	_data.clear();

	return NOERROR;
}

//------------------
// fini
//------------------
jerror_t DCCALGeometry_factory::fini(void)
{
	return NOERROR;
}

