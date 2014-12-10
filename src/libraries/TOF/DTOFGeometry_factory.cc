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

  myDTOFGeometry->NLAYERS        = 2;
  myDTOFGeometry->NENDS          = 2;

  myDTOFGeometry->NLONGBARS        = 42;
  myDTOFGeometry->NWIDEBARS        = 38;
  myDTOFGeometry->NSHORTBARS       = 4;
  myDTOFGeometry->LONGBARLENGTH    = 252.0;
  myDTOFGeometry->SHORTBARLENGTH   = 120.0;
  myDTOFGeometry->BARWIDTH         = 6.0;
  myDTOFGeometry->NBARS            = myDTOFGeometry->NLONGBARS + myDTOFGeometry->NSHORTBARS;

  myDTOFGeometry->FirstShortBar = 22;
  myDTOFGeometry->LastShortBar = myDTOFGeometry->NSHORTBARS/4 + myDTOFGeometry->FirstShortBar;


  myDTOFGeometry->CenterVPlane =  617.52;
  myDTOFGeometry->CenterHPlane =  620.10;
  myDTOFGeometry->CenterMPlane = (myDTOFGeometry->CenterVPlane +  myDTOFGeometry->CenterHPlane)/2.;

  // YPos[barnumber] gives y position centeral location of the bar. barnumber = 1 - 46
  for (int k=1;k<20;k++){ // first 19 long wide bars
    myDTOFGeometry->YPOS[k] = -129.0 + k*myDTOFGeometry->BARWIDTH;
  }
  for (int k=20;k<22;k++){ // then 2 long narrow bars
    myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[19] + (k-19+0.5)*myDTOFGeometry->BARWIDTH/2.;
  }
  for (int k=22;k<24;k++){ // then two short bars up/down left/right
    myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[21] + 0.5*myDTOFGeometry->BARWIDTH/2. +
      (k-21-0.5)*myDTOFGeometry->BARWIDTH;
  }
  for (int k=24;k<26;k++){ // then again 2 long narrow bars
    myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[23] + 0.5*myDTOFGeometry->BARWIDTH +
      (k-23-0.5)*myDTOFGeometry->BARWIDTH/2.;
  }
  for (int k=26;k<45;k++){ // then the last 19 long wide bars
    myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[25] + 0.5*myDTOFGeometry->BARWIDTH/2. +
      (k-25-0.5)*myDTOFGeometry->BARWIDTH;
  }
  for (int k=45;k<47;k++){ // then two short bars up/down left/right
    myDTOFGeometry->YPOS[k] = myDTOFGeometry->YPOS[21] + 0.5*myDTOFGeometry->BARWIDTH/2. +
      (k-21-0.5)*myDTOFGeometry->BARWIDTH;
  }
  
  _data.push_back(myDTOFGeometry);
    
  return NOERROR;
}  


