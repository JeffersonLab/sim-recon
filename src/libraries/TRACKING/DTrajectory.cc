// $Id$
//
//    File: DTrajectory.cc
// Created: Mon May 17 14:49:48 EDT 2010
// Creator: davidl (on Darwin eleanor.jlab.org 10.2.0 i386)
//

#include <cmath>
#include <iostream>
using namespace std;

#include "DTrajectory.h"

//---------------------------------
// DTrajectory    (Constructor)
//---------------------------------
DTrajectory::DTrajectory(const DMagneticFieldMap *bfield)
{
	Nswim_steps = 0;
	Max_swim_steps = 10000;
	swim_steps = new swim_step_t[Max_swim_steps];
	own_swim_steps = true;
	step_size = 2.0; // step size of zero means use adaptive step sizes
	
	ZMIN = -17.0;
	ZMAX = 650.0;
	RMAX = 88.0;
	R2MAX = RMAX*RMAX;
	mass = 0.0;
	
	this->bfield = bfield;
}

//---------------------------------
// ~DTrajectory    (Destructor)
//---------------------------------
DTrajectory::~DTrajectory()
{
	if(own_swim_steps)delete[] swim_steps;
}

//---------------------------------
// CalcDirs
//---------------------------------
void DTrajectory::CalcDirs(ThreeVector &pos, ThreeVector &p, RTdirs &dirs)
{
	// "x" and "p" give a position and momentum for which to
	// calculate the directions of the reference trajectory
	// coordinate system. The directions are copied into the
	// "dirs" structure along with the magnetic field at "x"
	// and Ro and the momentum components.
	//
	// Note that this is used by the Swim method which may pass
	// in a position not corresponding to the actual momentum.

	// Convenence references (these should cost nothing at run time)
	double &x = pos[0];
	double &y = pos[1];
	double &z = pos[2];
	double &px = dirs.p[0] = p[0];
	double &py = dirs.p[1] = p[1];
	double &pz = dirs.p[2] = p[2];
	double &Bx = dirs.B[0];
	double &By = dirs.B[1];
	double &Bz = dirs.B[2];
	double &xdir_x=dirs.xdir[0], &xdir_y=dirs.xdir[1], &xdir_z=dirs.xdir[2];
	double &ydir_x=dirs.ydir[0], &ydir_y=dirs.ydir[1], &ydir_z=dirs.ydir[2];
	double &zdir_x=dirs.zdir[0], &zdir_y=dirs.zdir[1], &zdir_z=dirs.zdir[2];
	double &Ro = dirs.Ro;
	double &sin_theta = dirs.sin_theta;
	double &cos_theta = dirs.cos_theta;
	
	// Get B-field at x
	bfield->GetField(x, y, z, Bx, By, Bz);
	//Bx = By = 0.0;
	//Bz = -2.0;
	
	// xdir = p cross B
	xdir_x = py*Bz - pz*Bx;
	xdir_y = pz*Bx - px*Bz;
	xdir_z = px*By - py*Bx;
	
	// zdir = Bdir
	zdir_x = Bx;
	zdir_y = By;
	zdir_z = Bz;
	
	// ydir = zdir cross xdir
	ydir_x = zdir_y*ydir_z - zdir_z*ydir_y;
	ydir_y = zdir_z*ydir_x - zdir_x*ydir_z;
	ydir_z = zdir_x*ydir_y - zdir_y*ydir_x;
	
	// Normalize all direction vectors. If any one of them has
	// a magnitude of zero or is not finite then default to lab
	// coordinate system. This may look ugly but it was written
	// to be optimized for speed.
	bool use_lab = false;
	double xmag, ymag, zmag;
	xmag = sqrt(xdir_x*xdir_x + xdir_y*xdir_y + xdir_z*xdir_z);
	use_lab |= (!finite(xmag) || xmag==0.0);
	if(!use_lab){
		ymag = sqrt(ydir_x*ydir_x + ydir_y*ydir_y + ydir_z*ydir_z);
		use_lab |= (!finite(ymag) || ymag==0.0);
		if(!use_lab){
			zmag = sqrt(zdir_x*zdir_x + zdir_y*zdir_y + zdir_z*zdir_z);
			use_lab |= (!finite(zmag) || zmag==0.0);
		}
	}
	if(use_lab){
		xdir_x=1.0;  xdir_y=0.0;  xdir_z=0.0;
		ydir_x=0.0;  ydir_y=1.0;  ydir_z=0.0;
		zdir_x=0.0;  zdir_y=0.0;  zdir_z=1.0;
		
		Ro = 1.0E20; // Effectively default to a straight line
	}else{
		xdir_x/=xmag;  xdir_y/=xmag;  xdir_z/=xmag;
		ydir_x/=ymag;  ydir_y/=ymag;  ydir_z/=ymag;
		zdir_x/=zmag;  zdir_y/=zmag;  zdir_z/=zmag;
		
		// Calculate Ro
		double &p_cross_B_mag=xmag; // p cross B already calculated above
		double B2 = Bx*Bx + By*By + Bz*Bz;
		Ro = p_cross_B_mag/B2/qBr2p;
	}
	
	// The values sin_theta and cos_theta represent the angle the 
	// momentum makes with the B-field. These are used to calculate
	// a helical step. Note that in the reference trajectory coordinate
	// system, the momentum should be in the y-z plane by definition
	double pmag = sqrt(px*px + py*py + pz*pz);
	cos_theta = (px*zdir_x + py*zdir_y + pz*zdir_z)/pmag;  // p_hat dot zdir
	sin_theta = (px*ydir_x + py*ydir_y + pz*ydir_z)/pmag;  // p_hat dot ydir
}

//---------------------------------
// CalcPosMom
//---------------------------------
void DTrajectory::CalcPosMom(double h, RTdirs &dirs, ThreeVector &pos, double *p)
{
	// Given the a step size "s" and a Reference Trajectory coordinate
	// system in "dirs", calculate the position after taking a step
	// along a perfect helix of the given size. The values are returned
	// in "x" and represent the delta of the step in lab coordinates. If
	// a non-NULL value is passed for "p", then the values in it are
	// updated to reflect the momentum at the end of the step. The values
	// at the begining of the step are contained in "dirs".
	//
	// NOTE: This is worth repeating: The values in "pos" will be a delta
	// in position which the caller will need to add to the starting
	// step position to get the actual location in lab coordinates. The
	// momentum values in "p" however, represent the actual momentum
	// (NOT the delta).
	
	// Convienence references (these should cost nothing at run time)
	double &x = pos[0];
	double &y = pos[1];
	double &z = pos[2];
	double &xdir_x=dirs.xdir[0], &xdir_y=dirs.xdir[1], &xdir_z=dirs.xdir[2];
	double &ydir_x=dirs.ydir[0], &ydir_y=dirs.ydir[1], &ydir_z=dirs.ydir[2];
	double &zdir_x=dirs.zdir[0], &zdir_y=dirs.zdir[1], &zdir_z=dirs.zdir[2];
	double &Ro = dirs.Ro;
	double &sin_theta = dirs.sin_theta;
	double &cos_theta = dirs.cos_theta;

	// Calculate step component parallel to B-field
	double delta_z = h*cos_theta;

	// Calculate step components perpendicular to B-field
	double delta_phi = h*sin_theta/Ro; // delta_phi is angle step makes in plane perpendicular to B
	//double cos_delta_phi = cos(delta_phi);
	//double sin_delta_phi = sin(delta_phi);
	double cos_delta_phi = 1.0 - (delta_phi*delta_phi/2.0);
	double sin_delta_phi = delta_phi;
	double delta_x = Ro*(1.0-cos_delta_phi);
	double delta_y = Ro*sin_delta_phi;
	
	// Calculate position (delta of) in lab coordinates
	x = delta_x*xdir_x + delta_y*ydir_x + delta_z*zdir_x;
	y = delta_x*xdir_y + delta_y*ydir_y + delta_z*zdir_y;
	z = delta_x*xdir_z + delta_y*ydir_z + delta_z*zdir_z;

//_DBG_<<"delta_x="<<delta_x<<" delta_y="<<delta_y<<" delta_z="<<delta_z<<" delta_phi="<<delta_phi<<" h="<<h<<endl;

	// If a non-NULL argument is passed in for "p" then the caller wants us also
	// to calculate the updated momentum at the end of the step. If they pass
	// NULL (default if argument is omitted), then they aren't interested in that
	// info and we can return now.
	if(p==NULL)return;
	
	//--------------- Calculate updated momentum vector --------------------
	// Calculate momentum by rotating it by delta_phi in the RT coordinate system
	// One would normally use a linear algebra package for this, but we try
	// and avoid some overhead by just coding it explicitly in order to optimize
	// speed.

	// Momentum stored in dirs is at starting point of step
	double &px = dirs.p[0];
	double &py = dirs.p[1];
	double &pz = dirs.p[2];
	
	// Momentum components in RT coordinate system
	double px_rt = px*xdir_x + py*xdir_y + pz*xdir_z;
	double py_rt = px*ydir_x + py*ydir_y + pz*ydir_z;
	double pz_rt = px*zdir_x + py*zdir_y + pz*zdir_z;
	
	// Rotate in X/Y plane of RT coordinate system
	double  px_rt_rotated =  px_rt*cos_delta_phi + py_rt*sin_delta_phi;
	double  py_rt_rotated = -py_rt*sin_delta_phi + py_rt*cos_delta_phi;
	double &pz_rt_rotated = pz_rt;

	// Convert back to lab coordinates
	p[0] = px_rt_rotated*xdir_x + py_rt_rotated*ydir_x + pz_rt_rotated*zdir_x;
	p[1] = px_rt_rotated*xdir_y + py_rt_rotated*ydir_y + pz_rt_rotated*zdir_y;
	p[2] = px_rt_rotated*xdir_z + py_rt_rotated*ydir_z + pz_rt_rotated*zdir_z;
}

//---------------------------------
// Swim
//---------------------------------
void DTrajectory::Swim(const DVector3 &pos, const DVector3 &mom, double q, double smax)
{
	swim_step_t *swim_step = swim_steps;
	
	// Copy starting parameters into first step
	swim_step->x = pos.x();
	swim_step->y = pos.y();
	swim_step->z = pos.z();
	swim_step->px = mom.x();
	swim_step->py = mom.y();
	swim_step->pz = mom.z();
	swim_step->s = 0.0;
	swim_step->t = 0.0;
	swim_step->dP = 0.0;
	swim_step->itheta02 = 0.0;
	swim_step->itheta02s = 0.0;
	swim_step->itheta02s2 = 0.0;
	
	bool DO_MATERIAL_LOSS = (mass!=0.0);

	double h = step_size;
	double h_2 = h/2.0;

	RTdirs dirs0, dirs1, dirs2, dirs3;
	
	double x0[3], x1[3], x2[3], x3[3];
	double k1[3], k2[3], k3[3], k4[3], k1_2[3], k2_2[3];
	double p1[3], p2[3], p3[3], p4[3];
	
	for(Nswim_steps=0; Nswim_steps<(Max_swim_steps-1); Nswim_steps++){
	
		// Setup references to make code more readable.
		// Note that these all point to the values at the start of the step
		double &x = swim_step->x;
		double &y = swim_step->y;
		double &z = swim_step->z;
		double &px = swim_step->px;
		double &py = swim_step->py;
		double &pz = swim_step->pz;
		double &s = swim_step->s;

		// Do 4th order Runge-Kutta
		//
		// Each of the following is a vector representing the full step in position.
		// k1 = projection of full step using Bfield at starting position
		// k2 = projection of full step using Bfield at mid-point of k1 step
		// k3 = projection of full step using Bfield at mid-point of k2 step
		// k4 = projection of full step using Bfield at end of k3 step
		//
		// The new position is given by the vector equation:
		//
		// x = x0 + k1/6 + k2/3 + k3/3 + k4/6
		//
		// The corresponding values for the momentum are simultaneously
		// calculated and stored in p1, p2, p3, and p4 such that
		//
		// p = p1/6 + p2/3 + p3/3 + p4/6
		//

		double p0[3] = {px, py, pz};

#if 0
		// k1
		x0[0]=x;  x0[1]=y;  x0[2]=z;
		CalcDirs(x0, p0, dirs0);
		CalcPosMom(h, dirs0, k1, p1);

		// k2
		CalcPosMom(h_2, dirs0, k1_2);
		x1[0]=x+k1_2[0];  x1[1]=y+k1_2[1];  x1[2]=z+k1_2[2];
		CalcDirs(x1, p0, dirs1);
		CalcPosMom(h, dirs1, k2, p2);
		
		// k3
		CalcPosMom(h_2, dirs1, k2_2);
		x2[0]=x+k2_2[0];  x2[1]=y+k2_2[1];  x2[2]=z+k2_2[2];
		CalcDirs(x2, p0, dirs2);
		CalcPosMom(h, dirs2, k3, p3);
		
		// k4
		x3[0]=x+k3[0];  x3[1]=y+k3[1];  x3[2]=z+k3[2];
		CalcDirs(x3, p0, dirs3);
		CalcPosMom(h, dirs3, k4, p4);
#endif

		// k1
		x0[0]=x;  x0[1]=y;  x0[2]=z;
		CalcDirs(x0, p0, dirs0);
		//CalcDirs(x0, p0, dirs0);

		// k2
		CalcPosMom(h_2, dirs0, k1_2);
		//CalcPosMom(h_2, dirs0, k1_2);
		x1[0]=x+k1_2[0];  x1[1]=y+k1_2[1];  x1[2]=z+k1_2[2];
		CalcDirs(x1, p0, dirs1);
		//CalcDirs(x1, p0, dirs1);
		CalcPosMom(h, dirs1, k2, p2);
		//CalcPosMom(h, dirs1, k2, p2);

		// Go to next swim step
		swim_step++;

#if 0
		// new_pos = pos + k1/6 + k2/3 + k3/3 + k4/6
		swim_step->x = x + (k1[0]+k4[0])/6.0 + (k2[0]+k3[0])/3.0;
		swim_step->y = y + (k1[1]+k4[1])/6.0 + (k2[1]+k3[1])/3.0;
		swim_step->z = z + (k1[2]+k4[2])/6.0 + (k2[2]+k3[2])/3.0;

		// new_mom = p1/6 + p2/3 + p3/3 + p4/6
		swim_step->px = (p1[0]+p4[0])/6.0 + (p2[0]+p3[0])/3.0;
		swim_step->py = (p1[1]+p4[1])/6.0 + (p2[1]+p3[1])/3.0;
		swim_step->pz = (p1[2]+p4[2])/6.0 + (p2[2]+p3[2])/3.0;
#endif
		// new_pos = pos + k2
		swim_step->x = x + k2[0];
		swim_step->y = y + k2[1];
		swim_step->z = z + k2[2];

		// new_mom = p2
		swim_step->px = p2[0];
		swim_step->py = p2[1];
		swim_step->pz = p2[2];

		
		// distance along trajectory
		swim_step->s = s+h;
		
		// Optionally account for material
		if(DO_MATERIAL_LOSS){
			_DBG__;
			bool particle_stopped = AdjustForMaterial(swim_step);
			if(particle_stopped) break;
		}else{
			// Set these values to reasonable defaults when not
			// adjusting for material. It's possible we could skip
			// this since these values should only be looked at
			// when material adjustments are in effect. That will
			// be left as a later optimization though.
			swim_step->t = 0.0;
			swim_step->dP = 0.0;
			swim_step->itheta02 = 0.0;
			swim_step->itheta02s = 0.0;
			swim_step->itheta02s2 = 0.0;
		}
		
		// Check if we have hit a boundary such that we should stop swimming
		if(swim_step->s >= smax)break; // max trajectory length
		if(swim_step->z <= ZMIN)break; // upstream z-boundary
		if(swim_step->z >= ZMAX)break; // downstream z-boundary
		double R2 = swim_step->x*swim_step->x + swim_step->y*swim_step->y;
		if(R2 >= R2MAX)break;  // r-boundary


//_DBG_<<"x,y,z="<<swim_step->x<<", "<<swim_step->y<<", "<<swim_step->z<<endl;
	}
	
	if(Nswim_steps >= (Max_swim_steps-1)){
		_DBG_<<"Maximum number of steps ("<<Max_swim_steps<<") reached. Swimming truncated."<<endl;
	}
}

//---------------------------------
// AdjustForMaterial
//---------------------------------
bool DTrajectory::AdjustForMaterial(swim_step_t *swim_step)
{
	swim_step_t *prev_swim_step = swim_step;
	prev_swim_step--; // potentially unsafe, but speedy
	
	// Total momentum at begining of step
	double P = sqrt(swim_step->px*swim_step->px + swim_step->py*swim_step->py + swim_step->pz*swim_step->pz);
	
	// Step size to integrate over
	double ds = swim_step->s - prev_swim_step->s;

	// Get Material properties for begining of step
	double radlen = 2.0E-5; // replace with call/calculation

	// Calculate momentum loss
	double dPdx = 0.0; // replace with call/calculation
	swim_step->dP = dPdx*ds;
	if(swim_step->dP >= P)swim_step->dP = P; // don't allow negative momentum

	// Apply momentum loss for the step
	double momentum_scale_factor = (P-swim_step->dP)/P;
	swim_step->px *= momentum_scale_factor;
	swim_step->py *= momentum_scale_factor;
	swim_step->pz *= momentum_scale_factor;

	// Update counters used to keep track of material for calculating MULS errors
	if(radlen>1.0E-5){ // PDG 2008 pg 271, second to last paragraph
		
		// n.b. we use P for previous step when should be that same as for the
		// current one before we adjusted it above.
		double beta=1./sqrt(1.+mass*mass/P/P);
		double theta0 = 0.0136/(P*beta)*sqrt(radlen)*(1.0+0.038*log(radlen)); // From PDG 2008 eq 27.12
		double theta02 = theta0*theta0;
		double &s = swim_step->s;
		swim_step->itheta02 = prev_swim_step->itheta02 + theta02;
		swim_step->itheta02s = prev_swim_step->itheta02s + s*theta02;
		swim_step->itheta02s2 = prev_swim_step->itheta02s2 + s*s*theta02;
	}
	
	// return true if particle stopped
	return (swim_step->dP >= P);
}

