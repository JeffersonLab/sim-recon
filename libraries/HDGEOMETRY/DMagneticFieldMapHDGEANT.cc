// $Id$
//
//    File: DMagneticFieldMapHDGEANT.cc
// Created: Thu Dec 21 14:03:31 EST 2006
// Creator: davidl (on Linux alkaid 2.6.9-42.0.2.ELsmp x86_64)
//

#include <iostream>
using namespace std;

#include "DMagneticFieldMapHDGEANT.h"

extern "C"{
 void gufld_(float*, float*); // FORTRAN
};

//---------------------------------
// DMagneticFieldMapHDGEANT    (Constructor)
//---------------------------------
DMagneticFieldMapHDGEANT::DMagneticFieldMapHDGEANT()
{

}

//---------------------------------
// ~DMagneticFieldMapHDGEANT    (Destructor)
//---------------------------------
DMagneticFieldMapHDGEANT::~DMagneticFieldMapHDGEANT()
{

}

//---------------------------------
// GetField
//---------------------------------
void DMagneticFieldMapHDGEANT::GetField(double x, double y, double z, double &Bx, double &By, double &Bz, int method=0) const
{
	float u[3], B[3];
	u[0] = x;
	u[1] = y;
	u[2] = z;
	gufld_(u,B);
	Bx = B[0]/10.0;
	By = B[1]/10.0;
	Bz = B[2]/10.0;
}
