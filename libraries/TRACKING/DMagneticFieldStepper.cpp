
#include <iostream>
#include <iomanip>
using namespace std;

#include "DMagneticFieldMap.h"
#include "DMagneticFieldStepper.h"

//-----------------------
// DMagneticFieldStepper
//-----------------------
DMagneticFieldStepper::DMagneticFieldStepper(const DMagneticFieldMap *bfield)
{
	this->bfield = bfield;
	q = 1.0;
	start_pos = pos = TVector3(0.0,0.0,0.0);
	start_mom = mom = TVector3(0.0,0.0,1.0);
	stepsize = 1.0; // in cm
}

//-----------------------
// DMagneticFieldStepper
//-----------------------
DMagneticFieldStepper::DMagneticFieldStepper(const DMagneticFieldMap *bfield, double q, TVector3 *x, TVector3 *p)
{
	this->bfield = bfield;
	this->q = q;
	start_pos = pos = *x;
	start_mom = mom = *p;
	stepsize = 1.0; // in cm
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
derror_t DMagneticFieldStepper::SetMagneticFieldMap(const DMagneticFieldMap *bfield)
{
	this->bfield = bfield;

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
	// x-axis is in direction perpendicular to both B and p (particle momentum)
	// y-axis is then just cross product of z and x axes.
	//
	// These coordinates are referred to as the natual coordinates below.
	// The step is calculated based on moving along a perfect helix a distance
	// of "stepsize". This means that the new position will actually be
	// closer than stepsize to the current position (unless the magnetic field
	// is zero).
	
	// Get B-field
	const DBfieldPoint_t* tmp = bfield->getQuick(pos.x(), pos.y(), pos.z());
	TVector3 B(tmp->Bx, tmp->By, tmp->Bz);

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
	double Rp = mom.Mag()/(q*B.Mag()*qBr2p); // qBr2p converts to GeV/c/cm so Rp will be in cm
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
	mom.Rotate(-delta_phi, B);
	
	// Energy loss for 1.0GeV pions in Air is roughly 2.4keV/cm
	double m = 0.13957;
	double p = mom.Mag();
	double E = sqrt(m*m + p*p) - 0.0000024*stepsize;
	mom.SetMag(sqrt(E*E-m*m));
	
	// return new position 
	if(newpos)*newpos = pos;

	return NOERROR;
}

//-----------------------
// GetBField
//-----------------------
const DBfieldPoint_t* DMagneticFieldStepper::GetDBfieldPoint(void)
{
	return bfield->getQuick(pos.x(), pos.y(), pos.z());
}

