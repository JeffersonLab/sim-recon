
#include <iostream>
#include <iomanip>
using namespace std;

#include "DMagneticFieldMap.h"
#include "DMagneticFieldStepper.h"

//-----------------------
// DMagneticFieldStepper
//-----------------------
DMagneticFieldStepper::DMagneticFieldStepper(DMagneticFieldMap *map, double q, TVector3 *x, TVector3 *p)
{
	bfield = map;
	this->q = q;
	start_pos = pos = *x;
	start_mom = mom = *p;
	stepsize = 10.0; // in cm
}

//-----------------------
// ~DMagneticFieldStepper
//-----------------------
DMagneticFieldStepper::~DMagneticFieldStepper()
{

}
	
//-----------------------
// SetStartingParams
//-----------------------
derror_t DMagneticFieldStepper::SetStartingParams(double q, TVector3 *x, TVector3 *p)
{
	this->q = q;
	start_pos = pos = *x;
	start_mom = mom = *p;	

	return NOERROR;
}

//-----------------------
// SetMagneticFieldMap
//-----------------------
derror_t DMagneticFieldStepper::SetMagneticFieldMap(DMagneticFieldMap *map)
{
	this->bfield = map;

	return NOERROR;
}

//-----------------------
// SetStepSize
//-----------------------
derror_t DMagneticFieldStepper::SetStepSize(double step)
{
	this->stepsize = step;

	return NOERROR;
}

//-----------------------
// Step
//-----------------------
derror_t DMagneticFieldStepper::Step(TVector3 *newpos)
{
	/// Advance the track one step and return the new position
	
	// The idea here is to work in the coordinate system whose
	// axes point in directions defined in the following way:
	// z-axis is along direction of B-field
	// x-axis is in direction of perpendicular to both B and p (particle momentum)
	// y-axis is then just cross product of z and x axes.
	//
	// These coordinates are referred to as the natual coordinates below.
	// The step is calculated based on moving along a perfect helix a distance
	// of "stepsize". This means that the new position will actually be
	// closer that stepsize to the current position (unless the magnetic field
	// is zero).
	
	// Get B-field
	float hbarc = 197.326;
	D3Vector_t tmp = bfield->getQuick(pos.x()/2.54, pos.y()/2.54, (66.0-30.0*2.54+pos.z())/2.54);
	TVector3 B(tmp.x/hbarc, tmp.y/hbarc, tmp.z/hbarc);
	//TVector3 B(0.0, 0.0, -2.0/hbarc);

	// If the magnetic field is zero or the charge is zero, then our job is easy
	if(B.Mag()==0.0 || q==0.0){
		TVector3 pstep = mom;
		pstep.SetMag(stepsize);
		pos += pstep;
		if(newpos)*newpos = pos;

		return NOERROR;
	}

	// cross product of p and B (natural x-direction)
	TVector3 xdir = mom.Cross(B);
	xdir.SetMag(1.0);
	
	// cross product of B and pxB (natural y-direction)
	TVector3 ydir = B.Cross(xdir);
	ydir.SetMag(1.0);
	
	// B-field is natural z-direction
	TVector3 zdir = B;
	zdir.SetMag(1.0);
	
	// cosine of angle between p and B
	double theta = B.Angle(mom);
	double cos_theta = cos(theta);
	double sin_theta = sin(theta);
	
	// delta_z is step size in z direction
	TVector3 delta_z = zdir*stepsize*cos_theta;
	
	// The ratio p/qB appears in a few places.
	double Rp = mom.Mag()/(q*B.Mag()*0.5931); // 0.5931 is fudge factor for now
	double Ro = Rp*sin_theta;

	// delta_phi is angle of rotation in natural x/y plane
	double delta_phi = stepsize/Rp;

	// delta_x is magnitude of step in natural x direction
	TVector3 delta_x = Ro*(1.0-cos(delta_phi))*xdir;
	
	// delta_y is magnitude of step in natural y direction
	TVector3 delta_y = Ro*sin(delta_phi)*ydir;
	
	// Step to new position
	pos += delta_x + delta_y + delta_z;
	
	// Update momentum by rotating it by delta_phi about B
	mom.Rotate(delta_phi, B);
	
	// return new position 
	if(newpos)*newpos = pos;

	return NOERROR;
}

//-----------------------
// Dist
//-----------------------
TVector3 DMagneticFieldStepper::Closest(TVector3 *v)
{
	/// This method not yet implemented

	return TVector3();
}

//-----------------------
// DMagneticFieldStepper
//-----------------------
double DMagneticFieldStepper::Dist(TVector3 *v)
{
	/// This method not yet implemented

	return 0.0;
}
	

