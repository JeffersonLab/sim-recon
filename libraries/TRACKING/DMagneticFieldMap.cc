// $Id$
//
//    File: DMagneticFieldMap.cc
// Created: Thu Dec 21 12:50:04 EST 2006
// Creator: davidl (on Linux alkaid 2.6.9-42.0.2.ELsmp x86_64)
//

#include "DMagneticFieldMap.h"

void DMagneticFieldMap::GetGradient(double x, double y, double z, double &dBdx, double &dBdy, double &dBdz, int method) const
{
	// Calculate gradient by finding difference in magnetic field vectors
	// at 6 points around the given point and taking average.
	double ds = 0.5;
	double B1[3], B2[3];
	
	double dBx[3];
	GetField(x+ds, y, z, B1[0], B1[1], B1[2]);
	GetField(x-ds, y, z, B2[0], B2[1], B2[2]);
	for(int i=0; i<3; i++)dBx[i] = B1[i]-B2[i];

	double dBy[3];
	GetField(x, y+ds, z, B1[0], B1[1], B1[2]);
	GetField(x, y-ds, z, B2[0], B2[1], B2[2]);
	for(int i=0; i<3; i++)dBy[i] = B1[i]-B2[i];

	double dBz[3];
	GetField(x, y, z+ds, B1[0], B1[1], B1[2]);
	GetField(x, y, z-ds, B2[0], B2[1], B2[2]);
	for(int i=0; i<3; i++)dBz[i] = B1[i]-B2[i];
	
	dBdx = (dBx[0]+dBy[0]+dBz[0])/3.0;
	dBdy = (dBx[1]+dBy[1]+dBz[1])/3.0;
	dBdy = (dBx[2]+dBy[2]+dBz[2])/3.0;
}

