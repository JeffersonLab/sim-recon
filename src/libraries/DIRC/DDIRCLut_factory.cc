// $Id$
//
//    File: DDIRCLut_factory.cc
//

#include <cassert>

#include "DDIRCLut_factory.h"
#include "DDIRCLut.h"

//------------------
// brun
//------------------
jerror_t DDIRCLut_factory::brun(JEventLoop *loop, int32_t runnumber)
{
	assert( _data.size() == 0 );

	flags = PERSISTANT;
	_data.push_back( new DDIRCLut() );
	
	return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DDIRCLut_factory::erun(void)
{
	for(unsigned int i=0; i<_data.size(); i++)delete _data[i];
	_data.clear();
	
	return NOERROR;
}

