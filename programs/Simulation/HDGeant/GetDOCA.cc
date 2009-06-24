

//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// WARNING!!
//
// This file is DISABLED by the #if 0 block below. It was used
// to enforce more consistency between the simulation and reconstuction
// but that is no longer needed (for the moment). I am commenting
// out the centents rather than removing the whole file so that it
// can be more easily accessed if we need it again. The price of 
// having this ability is that it requires a strong connection
// between the reconstruction and simulation that is undesirable
// with all else being equal.
//
// To re-enable this functionality, you'll need to not only un-comment
// this file, but also uncomment lines in the hitCDC.c and Makefile.bms
// (Makefile.orig) files.
//
// 6/24/2009 DL
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#if 0

#include <string.h>

#include <particleType.h>


extern "C" {
	void GetDOCA(int ipart, float x[3], float p[5], float doca[3]);
	void gfpart_(int *ipart, char * chnpar, int *itrtyp, float *amass, float *charge, float *tlife, float *ubuf, int *nubuf);

#include "geant3.h"
}


#include <HDGEOMETRY/DMagneticFieldMap.h>
#include <TRACKING/DReferenceTrajectory.h>

extern DMagneticFieldMap *Bmap; // from clibDB.cc

//-----------------
// GetDOCA
//------------------
void GetDOCA(int ipart, float x[3], float p[5], float doca[3])
{
	/* Project the particle from the given position with the given momentum
	 * to the DOCA as determined by the origin of the local coordinate system.
	 *
	 * The values of x and p are in the global coordinate system while doca
	 * is in the local coordinate system.
	 */

	/* Get the charge of this particle (jeesh!) */
	int nubuf, itrtyp;
	float ubuf[10], amass,charge,tlife;
	char chnpar[99];
	gfpart_(&ipart, chnpar,&itrtyp,&amass,&charge,&tlife,ubuf,&nubuf);

	// We need to define the coordinate system of the wire. The only way
	// we have of doing this (that I know of) is via the transformCoord()
	// routine.
	float origin[3] = {0.0, 0.0, 0.0};
	float sdir[3] = {1.0, 0.0, 0.0};
	float tdir[3] = {0.0, 1.0, 0.0};
	float udir[3] = {0.0, 0.0, 1.0};
	float origin_global[3], sdir_global[3], tdir_global[3], udir_global[3];
	transformCoord(origin,"local",origin_global,"global");
	transformCoord(sdir,"local",sdir_global,"global");
	transformCoord(tdir,"local",tdir_global,"global");
	transformCoord(udir,"local",udir_global,"global");
	DCoordinateSystem wire;
	wire.origin.SetXYZ(origin_global[0], origin_global[1], origin_global[2]);
	wire.sdir.SetXYZ(sdir_global[0], sdir_global[1], sdir_global[2]);
	wire.tdir.SetXYZ(tdir_global[0], tdir_global[1], tdir_global[2]);
	wire.udir.SetXYZ(udir_global[0], udir_global[1], udir_global[2]);
	wire.sdir -= wire.origin;
	wire.tdir -= wire.origin;
	wire.udir -= wire.origin;wire.L=200.0;

	// Create a "short" reference trajectory that uses only local memory
	// and will be swum just a couple of steps.
	DReferenceTrajectory::swim_step_t steps[64];
	DReferenceTrajectory rt(Bmap , charge , steps , 64);

	DVector3 pos(x[0], x[1], x[2]);
	DVector3 mom(p[0], p[1], p[2]);
	mom *= p[4];
	rt.Swim(pos, mom, charge, 2.0); // swim for a maximum of 2cm
	
	// Get the DOCA
	rt.DistToRT(&wire);
	rt.GetLastDOCAPoint(pos, mom);

	pos -= wire.origin;
	
	doca[0] = pos.Dot(wire.sdir);
	doca[1] = pos.Dot(wire.tdir);
	doca[2] = pos.Dot(wire.udir);
}


#endif
