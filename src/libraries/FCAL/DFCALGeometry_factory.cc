// $Id$
//
//    File: DFCALGeometry_factory.cc
// Created: Wed Aug 24 10:09:27 EST 2005
// Creator: shepherd (on Darwin 129-79-159-16.dhcp-bl.indiana.edu 8.2.0 powerpc)
//

#include <cassert>

#include "DFCALGeometry_factory.h"
#include "DFCALGeometry.h"

//------------------
// brun
//------------------
jerror_t DFCALGeometry_factory::brun(JEventLoop *loop, int runnumber)
{
	assert( _data.size() == 0 );

	flags = PERSISTANT;
	_data.push_back( new DFCALGeometry() );
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DFCALGeometry_factory::erun(void)
{
	for(unsigned int i=0; i<_data.size(); i++)delete _data[i];
	_data.clear();
	
	return NOERROR;
}

