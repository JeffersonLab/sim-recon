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

  myDTOFGeometry->NLONGBARS        = 44;
  myDTOFGeometry->NWIDEBARS        = 36;
  myDTOFGeometry->NSHORTBARS       = 4;
  myDTOFGeometry->LONGBARLENGTH    = 252.0;
  myDTOFGeometry->SHORTBARLENGTH   = 120.0;
  myDTOFGeometry->BARWIDTH         = 6.0;
  myDTOFGeometry->NBARS            = myDTOFGeometry->NLONGBARS + myDTOFGeometry->NSHORTBARS;

  myDTOFGeometry->CenterVPlane =  617.52;
  myDTOFGeometry->CenterHPlane =  620.10;
  myDTOFGeometry->CenterMPlane = (myDTOFGeometry->CenterVPlane +  myDTOFGeometry->CenterHPlane)/2.;


  for (int k=0;k<myDTOFGeometry->NBARS;k++){

    for (int k=0;k<18;k++){ // first 18 long wide bars
      myDTOFGeometry->YPOS[k] = -123.0 + k*myDTOFGeometry->BARWIDTH;
    }
    for (int k=18;k<22;k++){ // then 4 long narrow bars
      myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[17] + (k-17+0.5)*myDTOFGeometry->BARWIDTH/2.;
    }
    for (int k=22;k<24;k++){ // then two short bars up/down left/right
      myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[21] + 0.5*myDTOFGeometry->BARWIDTH/2. +
	(k-21-0.5)*myDTOFGeometry->BARWIDTH;
    }
    for (int k=24;k<28;k++){ // then again 4 long narrow bars
      myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[23] + 0.5*myDTOFGeometry->BARWIDTH +
	(k-23-0.5)*myDTOFGeometry->BARWIDTH/2.;
    }
    for (int k=28;k<46;k++){ // then the last 18 long wide bars
      myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[27] + 0.5*myDTOFGeometry->BARWIDTH/2. +
	(k-27-0.5)*myDTOFGeometry->BARWIDTH;
    }

  }

  _data.push_back(myDTOFGeometry);
    
  return NOERROR;
}  


