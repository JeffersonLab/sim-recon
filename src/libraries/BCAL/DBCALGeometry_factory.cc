// $Id$
//
//    File: DFactory_DBCALGeometry.cc
// Created: Thu Nov 17 15:10:51 CST 2005
// Creator: gluexuser (on Linux hydra.phys.uregina.ca 2.4.20-8smp i686)
//

#include <cassert>	

#include "DBCALGeometry_factory.h"

//------------------
// evnt
//------------------
jerror_t DBCALGeometry_factory::evnt(JEventLoop *loop, int eventnumber)
{

  DBCALGeometry *bcalGeom = new DBCALGeometry;
     
  _data.push_back(bcalGeom);

	return NOERROR;
}
