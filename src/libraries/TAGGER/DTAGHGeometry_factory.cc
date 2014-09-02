//
// File: DTAGHGeometry_factory.cc
// Created: Sat Jul 5 10:09:27 EST 2014
// Creator: jonesrt on gluey.phys.uconn.edu
//

#include "DTAGHGeometry_factory.h"
#include "DTAGHGeometry.h"

//------------------
// brun
//------------------
jerror_t DTAGHGeometry_factory::brun(JEventLoop *loop, int runnumber)
{
   assert( _data.size() == 0 );

   flags = PERSISTANT;
   _data.push_back( new DTAGHGeometry(loop, factory_tag, runnumber) );
   
   return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DTAGHGeometry_factory::erun(void)
{
   for (unsigned int i=0; i < _data.size(); i++)
      delete _data[i];
   _data.clear();
   
   return NOERROR;
}
