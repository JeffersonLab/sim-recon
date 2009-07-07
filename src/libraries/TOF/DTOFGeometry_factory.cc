// $Id$
//
//    File: DTOFGeometry_factory.cc
// Created: Mon Jul 18 11:43:31 EST 2005
// Creator: remitche (on Linux mantrid00 2.4.20-18.8 i686)
//

#include "DTOFGeometry_factory.h"

jerror_t DTOFGeometry_factory::init(void)
{

  flags = PERSISTANT;

  DTOFGeometry *myDTOFGeometry = new DTOFGeometry;

  myDTOFGeometry->NLONGBARS        = 40;
  myDTOFGeometry->NSHORTBARS       = 4;
  myDTOFGeometry->LONGBARLENGTH    = 252.0;
  myDTOFGeometry->SHORTBARLENGTH   = 120.0;
  myDTOFGeometry->BARWIDTH         = 6.0;
  
  _data.push_back(myDTOFGeometry);
    
  return NOERROR;
}  


