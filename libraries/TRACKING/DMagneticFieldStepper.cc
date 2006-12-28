
#include <iostream>
#include <iomanip>
using namespace std;

#include "DMagneticFieldMap.h"
#include "DMagneticFieldStepper.h"

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c

//-----------------------
// DMagneticFieldStepper
//-----------------------
DMagneticFieldStepper::DMagneticFieldStepper(const DMagneticFieldMap *bfield)
{
	this->bfield = bfield;
	q = 1.0;
	start_pos = pos = TVector3(0.0,0.0,0.0);
	start_mom = mom = TVector3(0.0,0.0,1.0);
	last_stepsize = stepsize = 1.0; // in cm
	CalcDirs();
}

//-----------------------
// DMagneticFieldStepper
//-----------------------
DMagneticFieldStepper::DMagneticFieldStepper(const DMagneticFieldMap *bfield, double q, const TVector3 *x, const TVector3 *p)
{
	this->bfield = bfield;
	this->q = q;
	start_pos = pos = *x;
	start_mom = mom = *p;
	last_stepsize = stepsize = 1.0; // in cm
	CalcDirs();
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
jerror_t DMagneticFieldStepper::SetStartingParams(double q, const TVector3 *x, const TVector3 *p)
{
	this->q = q;
	start_pos = pos = *x;
	start_mom = mom = *p;
	this->last_stepsize = this->stepsize;
	CalcDirs();

	return NOERROR;
}

//-----------------------
// SetMagneticFieldMap
//-----------------------
jerror_t DMagneticFieldStepper::SetMagneticFieldMap(const DMagneticFieldMap *bfield)
{
	this->bfield = bfield;

	return NOERROR;
}

//-----------------------
// SetStepSize
//-----------------------
jerror_t DMagneticFieldStepper::SetStepSize(double step)
{
	this->stepsize = step;
	this->last_stepsize = step;

	return NOERROR;
}

//-----------------------
// CalcDirs
//-----------------------
void DMagneticFieldStepper::CalcDirs(void)
{
	/// Calculate the directions of the "natural coordinates"
	/// (aka reference trajectory coordinate system) in the
	/// lab frame using the current momentum and magnetic field
	/// at the current position. The results are left in the
	/// private member fields, copies of which may be obtained
	/// by a subsequent call to the GetDirs(...) and GetRo()
	/// methods.

	// Get B-field
	double Bx,By,Bz;
	bfield->GetField(pos.x(), pos.y(), pos.z(), Bx, By, Bz);
	B.SetXYZ(Bx, By, Bz);

	// If the B-field is zero, then default to lab system
	double B2 = B.Mag2();
	if(B2==0.0){
		xdir.SetXYZ(1.0,0.0,0.0);
		ydir.SetXYZ(0.0,1.0,0.0);
		zdir.SetXYZ(0.0,0.0,1.0);
		cos_theta = 1.0;
		sin_theta = 0.0;
		Rp = Ro = 1.0E20; // This shouldn't really be used ever
		return;
	}

	// cross product of p and B (natural x-direction)
	xdir = mom.Cross(B);
	xdir.SetMag(1.0);
	
	// cross product of B and pxB (natural y-direction)
	ydir = B.Cross(xdir);
	ydir.SetMag(1.0);
	
	// B-field is natural z-direction
	zdir = B;
	zdir.SetMag(1.0);

	// cosine of angle between p and B
	double theta = B.Angle(mom);
	cos_theta = cos(theta);
	sin_theta = sin(theta);

	// Calculate Rp and Ro for this momentum and position
	Rp = mom.Mag()/(q*sqrt(B2)*qBr2p); // qBr2p converts to GeV/c/cm so Rp will be in cm
	Ro = Rp*sin_theta;
}

#if 1
//extern "C" {
void grkuta_(double *CHARGE, double *STEP, double *VECT, double *VOUT,const DMagneticFieldMap *bfield);
//}

//-----------------------
// Step
//-----------------------
double DMagneticFieldStepper::Step(TVector3 *newpos, double stepsize)
{
	double VECT[7], VOUT[7];
	VECT[0] = pos.x();
	VECT[1] = pos.y();
	VECT[2] = pos.z();
	VECT[6] = mom.Mag();
	VECT[3] = mom.x()/VECT[6];
	VECT[4] = mom.y()/VECT[6];
	VECT[5] = mom.z()/VECT[6];
	if(stepsize==0.0)stepsize = this->stepsize;
	double &STEP = stepsize;
	double &CHARGE = q;

	// Call GEANT3 Runge Kutta step routine (see grkuta.c)
	grkuta_(&CHARGE, &STEP, VECT, VOUT, bfield);

	pos.SetXYZ(VOUT[0], VOUT[1], VOUT[2]);
	mom.SetXYZ(VOUT[3], VOUT[4], VOUT[5]);
	mom.SetMag(VOUT[6]);

	CalcDirs();

	if(newpos)*newpos = pos;

	return STEP;
}

#else

//-----------------------
// Step
//-----------------------
double DMagneticFieldStepper::Step(TVector3 *newpos, double stepsize)
{
	/// Advance the track one step. Copy the new position into the
	/// TVector3 pointer (if given). Returns distance along path
	/// traversed in step.
	
	// The idea here is to work in the coordinate system whose
	// axes point in directions defined in the following way:
	// z-axis is along direction of B-field
	// x-axis is in direction perpendicular to both B and p (particle momentum)
	// y-axis is then just cross product of z and x axes.
	//
	// These coordinates are referred to as the natural coordinates below.
	// The step is calculated based on moving along a perfect helix a distance
	// of "stepsize". This means that the new position will actually be
	// closer than stepsize to the current position (unless the magnetic field
	// is zero in which case they are equal).
	//
	// The values of the helix needed to project the step are
	// calculated for the current momentum and position in CalcDirs()
	// above. They should already be corect for the step we are
	// about to take. We call CalcDirs() after updating the position
	// and momentum in order to prepare for the next step.
	
	// If we weren't passed a step size, then calculate it using
	// info stored in the member data
	if(stepsize==0.0){
		// We always try growing the step size a little
		stepsize = 1.25*this->last_stepsize;
		
		// Don't grow past the set stepsize
		if(stepsize>this->stepsize)stepsize=this->stepsize;
	}

	// If the magnetic field is zero or the charge is zero, then our job is easy
	if(B.Mag2()==0.0 || q==0.0){
		TVector3 pstep = mom;
		pstep.SetMag(stepsize);
		pos += pstep;
		if(newpos)*newpos = pos;
		CalcDirs();

		return stepsize;
	}

	// Note that cos_theta, sin_theta, Rp, Ro, xdir, ydir, and zdir
	// should already be valid for the starting point of this step.
	// They are set via a call to CalcDirs()

	// delta_z is step size in z direction
	TVector3 delta_z = zdir*stepsize*cos_theta;

	// delta_phi is angle of rotation in natural x/y plane
	double delta_phi = stepsize/Rp;

	// delta_x is magnitude of step in natural x direction
	TVector3 delta_x = Ro*(1.0-cos(delta_phi))*xdir;
	
	// delta_y is magnitude of step in natural y direction
	TVector3 delta_y = Ro*sin(delta_phi)*ydir;
	
	// Calculate new position
	TVector3 mypos = pos + delta_x + delta_y + delta_z;

	// Adjust the step size if needed.
	// Check that the field hasn't changed too much here. If it has,
	// reduce the step size and try again.
	// We actually check three conditions:
	// 1.) That the step size is not already too small
	// 2.) that the radius of curvature is still not so big
	//     that it  is worth reducing the step size
	// 3.) that the difference of the field between the old
	//     and new points is small relative to the field size
	if(stepsize>0.1){
		if(fabs(stepsize/Ro) > 1.0E-4){
			double Bx,By,Bz;
			bfield->GetField(mypos.x(), mypos.y(), mypos.z(), Bx, By, Bz);
			TVector3 Bnew(Bx, By, Bz);
			TVector3 Bdiff = B-Bnew;
			double f = Bdiff.Mag()/B.Mag();
			if(f > 1.0E-4){
				return Step(newpos, stepsize/2.0);
			}
		}
	}

	// Step to new position
	pos = mypos;

	// Update momentum by rotating it by delta_phi about B
	mom.Rotate(-delta_phi, B);
	
	// Energy loss for 1.0GeV pions in Air is roughly 2.4keV/cm
	//double m = 0.13957;
	//double p = mom.Mag();
	//double E = sqrt(m*m + p*p) - 0.0000024*stepsize;
	//mom.SetMag(sqrt(E*E-m*m));

	// Calculate directions of natural coordinates at the new position
	CalcDirs();
	
	// return new position 
	if(newpos)*newpos = pos;
	
	// Record this step size for (possible) use on the next iteration
	last_stepsize = stepsize;

	return stepsize;
}
#endif

//-----------------------
// GetDirs
//-----------------------
void DMagneticFieldStepper::GetDirs(TVector3 &xdir, TVector3 &ydir, TVector3 &zdir)
{
	xdir = this->xdir;
	ydir = this->ydir;
	zdir = this->zdir;
}


