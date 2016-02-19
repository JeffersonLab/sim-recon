// $Id$
//
//    File: DFactory_DBCALGeometry.cc
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>	

#include "DBCALGeometry_factory.h"


//------------------
// brun
//------------------
jerror_t DBCALGeometry_factory::brun(JEventLoop *loop, int32_t runnumber)
{
    assert( _data.size() == 0 );

    flags = PERSISTANT;
    _data.push_back( new DBCALGeometry(runnumber) );
        
    return NOERROR;
}

//------------------
// erun
//------------------
jerror_t DBCALGeometry_factory::erun(void)
{
    for(unsigned int i=0; i<_data.size(); i++)delete _data[i];
    _data.clear();
        
    return NOERROR;
}
