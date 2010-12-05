//    File: DBCALGeometry.h
// Created: Fri Nov 26 15:10:51 CST 2010
// Creator: dwbennet

#include "DBCALGeometry.h"

int DBCALGeometry::NBCALMODS = 48;
int DBCALGeometry::NBCALLAYS1 =  6;
int DBCALGeometry::NBCALLAYS2 =  4; 
int DBCALGeometry::NBCALSECS1 =  4; 
int DBCALGeometry::NBCALSECS2 =  4;
int DBCALGeometry::BCALMID = 7; //Enter the first layer in "outer cells" (default 7)
float DBCALGeometry::BCALINNERRAD = 64.3;   
float DBCALGeometry::BCALOUTERRAD = 86.17;
float DBCALGeometry::BCALFIBERLENGTH = 390.0;
float DBCALGeometry::GLOBAL_CENTER = 212;
  
float DBCALGeometry::ATTEN_LENGTH = 300.;
   
float DBCALGeometry::C_EFFECTIVE  = 16.75;

float DBCALGeometry::m_radius[] = { 64.3, 
				  66.3,
				  68.3,
				  70.3,
				  72.3,
				  74.3,
				  76.3,
				  78.77,
				  81.24,
				  83.70,
				  86.17};

DBCALGeometry::DBCALGeometry()  
{
  BCALMIDRAD = m_radius[BCALMID-1];
  /// End if groupings do not evenly divide SiPM cells
  if(  (BCALMID-1) % NBCALLAYS1 != 0 
       || (11-BCALMID) % NBCALLAYS2 != 0
       ||  4 % NBCALSECS1 != 0
       ||  4 % NBCALSECS2 != 0)
    {
      std::cout<<"ERROR: Bad BCAL fADC groupings, do not evenly divide cells";
      assert (false);
    }
}
