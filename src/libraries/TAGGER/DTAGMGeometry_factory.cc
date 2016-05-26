//
// File: DTAGMGeometry_factory.cc
// Created: Sat Jul 5 10:09:27 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#include "DTAGMGeometry_factory.h"
#include "DTAGMGeometry.h"

using namespace std;

//------------------
// brun
//------------------
jerror_t DTAGMGeometry_factory::brun(JEventLoop *loop, int32_t runnumber)
{
	if(!_data.empty())
	{
		//for change in run #
		delete _data[0];
		_data.clear();
	}

   flags = PERSISTANT;
   _data.push_back( new DTAGMGeometry(loop) );
   
   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTAGMGeometry_factory::erun(void)
{
   for (unsigned int i=0; i < _data.size(); i++)
      delete _data[i];
   _data.clear();
   
   return NOERROR;
}
