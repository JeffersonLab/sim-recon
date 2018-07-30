// $Id$
//
//    File: DTOFGeometry_factory.cc
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#include "DTOFGeometry_factory.h"
#include <HDGEOMETRY/DGeometry.h>
#include <DANA/DApplication.h>

jerror_t DTOFGeometry_factory::brun(jana::JEventLoop *loop, int32_t runnumber)
{

  flags = PERSISTANT;

  // Get the geometry
  DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
  const DGeometry* locGeometry = dapp->GetDGeometry(runnumber);

  DTOFGeometry *myDTOFGeometry = new DTOFGeometry(locGeometry);
   
  _data.push_back(myDTOFGeometry);
    
  return NOERROR;
}  


