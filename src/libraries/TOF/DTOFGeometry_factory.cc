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
   
  // Get the geometry
  DApplication* dapp = dynamic_cast<DApplication*>(loop->GetJApplication());
  DGeometry* locGeometry = dapp->GetDGeometry(runnumber);

  // Store the z position for both planes
  vector<double>tof_face;
  locGeometry->Get("//section/composition/posXYZ[@volume='ForwardTOF']/@X_Y_Z",tof_face);
  vector<double>tof_plane0;
  locGeometry->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='0']", tof_plane0);
  vector<double>tof_plane1;
  locGeometry->Get("//composition[@name='ForwardTOF']/posXYZ[@volume='forwardTOF']/@X_Y_Z/plane[@value='1']", tof_plane1);
  myDTOFGeometry->CenterVPlane=tof_face[2]+tof_plane1[2];
  myDTOFGeometry->CenterHPlane=tof_face[2]+tof_plane0[2];
  myDTOFGeometry->CenterMPlane=0.5*(myDTOFGeometry->CenterHPlane
				   +myDTOFGeometry->CenterVPlane);

  // Next fill array of bar positions within a plane
  // YPos[barnumber] gives y position centeral location of the bar. barnumber = 1 - 46
  double y0,dy;

  // First 19 long bars
  locGeometry->Get("//composition[@name='forwardTOF_bottom1']/mposY/@Y0",y0);
  locGeometry->Get("//composition[@name='forwardTOF_bottom1']/mposY/@dY",dy);
  vector<double>tof_bottom1;
  locGeometry->Get("//composition[@name='forwardTOF']/posXYZ[@volume='forwardTOF_bottom1']/@X_Y_Z",tof_bottom1);  
  for (unsigned int k=1;k<20;k++){
    myDTOFGeometry->YPOS[k]=y0+tof_bottom1[1]+dy*double(k-1);
  }
  
  // two narrow long bars
  locGeometry->Get("//composition[@name='forwardTOF_bottom2']/mposY/@Y0",y0);
  locGeometry->Get("//composition[@name='forwardTOF_bottom2']/mposY/@dY",dy);
  vector<double>tof_bottom2;
  locGeometry->Get("//composition[@name='forwardTOF']/posXYZ[@volume='forwardTOF_bottom2']/@X_Y_Z",tof_bottom2);  
  for (unsigned int k=20;k<22;k++){
    myDTOFGeometry->YPOS[k]=y0+tof_bottom2[1]+dy*double(k-20);
  }

  // two short wide bars
  locGeometry->Get("//composition[@name='forwardTOF_north']/mposY/@Y0",y0);
  locGeometry->Get("//composition[@name='forwardTOF_north']/mposY/@dY",dy);
  vector<double>tof_north;
  locGeometry->Get("//composition[@name='forwardTOF']/posXYZ[@volume='forwardTOF_north']/@X_Y_Z",tof_north);  
  for (unsigned int k=22;k<24;k++){
    myDTOFGeometry->YPOS[k]=y0+tof_north[1]+dy*double(k-22);
  }

  // two narrow long bars
  locGeometry->Get("//composition[@name='forwardTOF_top2']/mposY/@Y0",y0);
  locGeometry->Get("//composition[@name='forwardTOF_top2']/mposY/@dY",dy);
  vector<double>tof_top2;
  locGeometry->Get("//composition[@name='forwardTOF']/posXYZ[@volume='forwardTOF_top2']/@X_Y_Z",tof_top2);  
  for (unsigned int k=24;k<26;k++){
    myDTOFGeometry->YPOS[k]=y0+tof_top2[1]+dy*double(k-24);
  }

  // Last 19 long bars
  locGeometry->Get("//composition[@name='forwardTOF_top1']/mposY/@Y0",y0);
  locGeometry->Get("//composition[@name='forwardTOF_top1']/mposY/@dY",dy);
  vector<double>tof_top1;
  locGeometry->Get("//composition[@name='forwardTOF']/posXYZ[@volume='forwardTOF_top1']/@X_Y_Z",tof_top1);  
  for (unsigned int k=26;k<45;k++){
    myDTOFGeometry->YPOS[k]=y0+tof_top1[1]+dy*double(k-26);
  }

  // two more short wide bars
  locGeometry->Get("//composition[@name='forwardTOF_south']/mposY/@Y0",y0);
  locGeometry->Get("//composition[@name='forwardTOF_south']/mposY/@dY",dy);
  vector<double>tof_south;
  locGeometry->Get("//composition[@name='forwardTOF']/posXYZ[@volume='forwardTOF_south']/@X_Y_Z",tof_south);  
  for (unsigned int k=45;k<47;k++){
    myDTOFGeometry->YPOS[k]=y0+tof_south[1]+dy*double(k-45);
  }

  _data.push_back(myDTOFGeometry);
    
  return NOERROR;
}  


