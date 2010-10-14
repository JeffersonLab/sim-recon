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

   bcalGeom->NBCALMODS  = 48;
   bcalGeom->NBCALLAYS1 =  2;
   bcalGeom->NBCALLAYS2 =  2; 
   bcalGeom->NBCALSECS1 =  4; 
   bcalGeom->NBCALSECS2 =  2;
   bcalGeom->BCALINNERRAD = 64.3;   
   bcalGeom->BCALMIDRAD   = 76.3;
   bcalGeom->BCALOUTERRAD = 86.17;
   bcalGeom->BCALFIBERLENGTH = 390.0;
   bcalGeom->GLOBAL_CENTER = 212;
   
   bcalGeom->ATTEN_LENGTH = 300.;
   
   bcalGeom->C_EFFECTIVE  = 16.75;
   
  _data.push_back(bcalGeom);


	return NOERROR;
}
