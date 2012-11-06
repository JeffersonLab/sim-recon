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
DTrajectory::DTrajectory(const DMagneticFieldMap *bfield, const DGeometry *geom)
{
	Nswim_steps = 0;
	Max_swim_steps = 10000;
	swim_steps = new swim_step_t[Max_swim_steps];
	own_swim_steps = true;
	
	// Initialize some values from configuration parameters
	ZMIN = -17.0;
	ZMAX = 650.0;
	RMAX = 88.0;
	mass = 0.1396;
	BOUNDARY_STEP_FRACTION = 0.80;
	MIN_STEP_SIZE = 0.05;	// cm
	MAX_STEP_SIZE = 3.0;		// cm
	
	
	R2MAX = RMAX*RMAX;
	
	check_material_boundaries = true;
	
	this->bfield = bfield;
	this->geom = geom;
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
	bool   &straight_track = dirs.straight_track;
	
	// Get B-field at x
	bfield->GetField(x, y, z, Bx, By, Bz);
	//Bx = By = 0.0;
	//Bz = -2.0;
	
	// xdir = p cross B
	xdir_x = py*Bz - pz*By;
	xdir_y = pz*Bx - px*Bz;
	xdir_z = px*By - py*Bx;
	
	// zdir = Bdir
	zdir_x = Bx;
	zdir_y = By;
	zdir_z = Bz;
	
	// ydir = zdir cross xdir
	ydir_x = zdir_y*xdir_z - zdir_z*xdir_y;
	ydir_y = zdir_z*xdir_x - zdir_x*xdir_z;
	ydir_z = zdir_x*xdir_y - zdir_y*xdir_x;
	
	// Normalize all direction vectors. If any one of them has
	// a magnitude of zero or is not finite then set the 
	// "straight_track" flag and return immediately
	double xmag = sqrt(xdir_x*xdir_x + xdir_y*xdir_y + xdir_z*xdir_z);
	if(!finite(xmag) || xmag==0.0){straight_track=true; return;}
	
	double ymag = sqrt(ydir_x*ydir_x + ydir_y*ydir_y + ydir_z*ydir_z);
	if(!finite(ymag) || ymag==0.0){straight_track=true; return;}
	
	double zmag = sqrt(zdir_x*zdir_x + zdir_y*zdir_y + zdir_z*zdir_z);
	if(!finite(zmag) || zmag==0.0){straight_track=true; return;}
	
	// It looks like we're not a straight track
	straight_track = false;

	// Apply normalization to direction vectors
	xdir_x/=xmag;  xdir_y/=xmag;  xdir_z/=xmag;
	ydir_x/=ymag;  ydir_y/=ymag;  ydir_z/=ymag;
	zdir_x/=zmag;  zdir_y/=zmag;  zdir_z/=zmag;
		
	// Calculate Ro
	double &p_cross_B_mag=xmag; // p cross B already calculated above
	double B2 = Bx*Bx + By*By + Bz*Bz;
	Ro = p_cross_B_mag/B2/qBr2p;
	
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
	bool   &straight_track = dirs.straight_track;

	// First check if we're a straight track or not
	if(straight_track){
		double &dx = dirs.p[0];
		double &dy = dirs.p[1];
		double &dz = dirs.p[2];
		double f = h/sqrt(dx*dx + dy*dy + dz*dz);
		x = f*dx;
		y = f*dy;
		z = f*dz;
		if(p){
			p[0] = dirs.p[0];
			p[1] = dirs.p[1];
			p[2] = dirs.p[2];
		}
		return;
	}

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
	double  py_rt_rotated = -px_rt*sin_delta_phi + py_rt*cos_delta_phi;
	double &pz_rt_rotated =  pz_rt;

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
#if 0  // disable for now
	swim_step_t *swim_step = swim_steps;
	
	// Copy starting parameters into first step
	swim_step->x = pos.x();
	swim_step->y = pos.y();
	swim_step->z = pos.z();
	swim_step->px = mom.x();
	swim_step->py = mom.y();
	swim_step->pz = mom.z();
	swim_step->P = mom.Mag();
	swim_step->s = 0.0;
	swim_step->t = 0.0;
	swim_step->dP = 0.0;
	swim_step->itheta02 = 0.0;
	swim_step->itheta02s = 0.0;
	swim_step->itheta02s2 = 0.0;
	
	bool DO_MATERIAL_LOSS = (mass!=0.0) && (geom!=NULL);

	RTdirs dirs0, dirs1, dirs2, dirs3;
	
	double x0[3], x1[3], x2[3], x3[3];
	double k1[3], k2[3], k3[3], k4[3], k1_2[3], k2_2[3];
	double p1[3], p2[3], p3[3], p4[3];
	
	for(Nswim_steps=1; Nswim_steps<(Max_swim_steps-1); Nswim_steps++){
	
		// Setup references to point to begining of step
		// (and to make code more readable.)
		double &x  = swim_step->x;
		double &y  = swim_step->y;
		double &z  = swim_step->z;
		double &px = swim_step->px;
		double &py = swim_step->py;
		double &pz = swim_step->pz;
		double &s  = swim_step->s;
		
		double p0[3] = {px, py, pz}; // momentum vector at start of step
		
		// Increment to the end step. Note that the references
		// defined above will still point to the start step.
		swim_step_t *start_step=swim_step;
		swim_step++;
		swim_step_t *end_step = swim_step;

		// Initialize step size to maximum so it can be overwritten based
		// on material losses below
		double h = MAX_STEP_SIZE;

		// The total momentum of the particle is kept in local variable "P"
		// We set it here to the value at the beginning of the step. If
		// material losses are being included below, they will decrease it
		// by half of the momentum loss of the step so that the step is
		// actually calculated using the median momentum of the particle
		// over the step. Values used for calculating the MULS errors will
		// also use the median momentum for the step.
		double P = start_step->P;

		// Get the dPdx and radiation length of the material at the 
		// beginning of the step (if we are including material losses
		// in the swimming.)
		double momentum_scale_factor = 1.0;
		if(DO_MATERIAL_LOSS){
			double dP_dx, X0, s_to_boundary=MAX_STEP_SIZE;
			DVector3 pos(x, y, z);
			DVector3 mom(px, py, pz);
			bool particle_stopped = GetMaterialInfo(P, pos, mom, dP_dx, X0, s_to_boundary);
			if(particle_stopped) break;
			
			// Calculate step size based on 100keV/c momentum loss
			double step_size_ploss = 0.0001/dP_dx;
			
			// Keep the minimum step size
			if(s_to_boundary   < h) h =s_to_boundary;
			if(step_size_ploss < h) h =step_size_ploss;
			if( h < MIN_STEP_SIZE ) h = MIN_STEP_SIZE;
			
			// Use step size to calculate momentum loss and radiation lengths
			double &dP = end_step->dP;
			dP = dP_dx*h;
			double radlen = h/X0;
			
			// Adjust momentum vector to have a magnitude equal to the median
			// magnitude over the step
			momentum_scale_factor = (start_step->P - dP/2.0)/start_step->P;
			P *= momentum_scale_factor;
			p0[0] *= momentum_scale_factor;
			p0[1] *= momentum_scale_factor;
			p0[2] *= momentum_scale_factor;
			
			// Update counters used to keep track of material for calculating MULS errors
			if(radlen>1.0E-5){ // PDG 2008 pg 271, second to last paragraph
				// n.b. we use P for previous step when should be that same as for the
				// current one before we adjusted it above.
				double beta=1./sqrt(1.+mass*mass/P/P);
				double theta0 = 0.0136/(P*beta)*sqrt(radlen)*(1.0+0.038*log(radlen)); // From PDG 2008 eq 27.12
				double theta02 = theta0*theta0;
				double s = start_step->s + h;
				end_step->itheta02 = start_step->itheta02 + theta02;
				end_step->itheta02s = start_step->itheta02s + s*theta02;
				end_step->itheta02s2 = start_step->itheta02s2 + s*s*theta02;
			}
		}else{
			// Set these values to reasonable defaults when not
			// adjusting for material. It's possible we could skip
			// this since these values should only be looked at
			// when material adjustments are in effect. That will
			// be left as a later optimization though.
			end_step->t = 0.0;
			end_step->dP = 0.0;
			end_step->itheta02 = 0.0;
			end_step->itheta02s = 0.0;
			end_step->itheta02s2 = 0.0;
		}
		
		double h_2 = h/2.0;

		// Do either 2nd or 4th order Runge-Kutta
		//
		// For 4th order Runge-Kutta:
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


#if 1	// 1=midpoint method  0= 4th order Runge-Kutta

		// Midpoint (2nd order Runge-Kutta)

		// k1
		x0[0]=x;  x0[1]=y;  x0[2]=z;
		CalcDirs(x0, p0, dirs0);

		// k2
		CalcPosMom(h_2, dirs0, k1_2);
		x1[0]=x+k1_2[0];  x1[1]=y+k1_2[1];  x1[2]=z+k1_2[2];
		CalcDirs(x1, p0, dirs1);
		CalcPosMom(h, dirs1, k2, p2);

		// new_pos = pos + k2
		end_step->x = x + k2[0];
		end_step->y = y + k2[1];
		end_step->z = z + k2[2];

		// new_mom = p2
		end_step->px = p2[0];
		end_step->py = p2[1];
		end_step->pz = p2[2];

#else
		// 4th order Runge-Kutta

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

		// new_pos = pos + k1/6 + k2/3 + k3/3 + k4/6
		end_step->x = x + (k1[0]+k4[0])/6.0 + (k2[0]+k3[0])/3.0;
		end_step->y = y + (k1[1]+k4[1])/6.0 + (k2[1]+k3[1])/3.0;
		end_step->z = z + (k1[2]+k4[2])/6.0 + (k2[2]+k3[2])/3.0;

		// new_mom = p1/6 + p2/3 + p3/3 + p4/6
		end_step->px = (p1[0]+p4[0])/6.0 + (p2[0]+p3[0])/3.0;
		end_step->py = (p1[1]+p4[1])/6.0 + (p2[1]+p3[1])/3.0;
		end_step->pz = (p1[2]+p4[2])/6.0 + (p2[2]+p3[2])/3.0;
#endif

		
		// Update distance along trajectory
		end_step->s = start_step->s + h;
		
		// If including material losses, adjust momemtum
		end_step->px *= momentum_scale_factor;
		end_step->py *= momentum_scale_factor;
		end_step->pz *= momentum_scale_factor;
		end_step->P = P*momentum_scale_factor;

		// Check if we have hit a boundary such that we should stop swimming
		if(end_step->s >= smax)break; // max trajectory length
		if(end_step->z <= ZMIN)break; // upstream z-boundary
		if(end_step->z >= ZMAX)break; // downstream z-boundary
		double R2 = end_step->x*end_step->x + end_step->y*end_step->y;
		if(R2 >= R2MAX)break;  // r-boundary

	}

	// If we broke out of loop, then increment Nswim_steps
	if(Nswim_steps<Max_swim_steps)Nswim_steps++;
	
	if(Nswim_steps > Max_swim_steps){
		_DBG_<<"Maximum number of steps ("<<Max_swim_steps<<") reached. Swimming truncated."<<endl;
	}
#endif
}

//---------------------------------
// GetMaterialInfo
//---------------------------------
int DTrajectory::GetMaterialInfo(double P, DVector3 &pos, DVector3 &mom, double &dP_dx, double &X0, double &s_to_boundary)
{
	/// Determine the material properties for the given position and momentum
	/// and use them to calculate the rate of momemtum loss per cm, radiation
	/// length of material, and estimated distance to material boundary.
	///
	/// It may seem redundant to pass in both the magnitude of the momentum 
	/// as well as the momentum vector. Do this is just an optimization that
	/// saves a semi-expensive calculation of the magnitude here.
	///
	/// Momentum loss per cm of the material is returned in GeV/c per cm
	/// This follows the July 2008 PDG section 27.2 ppg 268-270.
	///
	/// Return value is value returned by call to DGeometry::FindMatALT1


	// Get Material properties for begining of step
	double KrhoZ_overA, LogI, rhoZ_overA;
	int err;
	if(check_material_boundaries){
		err = geom->FindMatALT1(pos, mom, KrhoZ_overA, rhoZ_overA,LogI, X0, &s_to_boundary);
	}else{
		err = geom->FindMatALT1(pos, mom, KrhoZ_overA, rhoZ_overA,LogI, X0);
	}

	// Calculate rate of momentum loss

	// Since we are dividing by mass below we should verify that mass is not zero
	if(mass==0.0){
		_DBG_<<"DTrajectory::GetMaterialInfo called for mass=0.0 particle!"<<endl;
		dP_dx = 0.0;
		return 0; // no ionization losses for neutrals
	}
	
	double gammabeta = P/mass;
	double gammabeta2=gammabeta*gammabeta;
	double gamma = sqrt(gammabeta2+1);
	double beta = gammabeta/gamma;
	double beta2=beta*beta;
	double me = 0.511E-3;
	double m_ratio=me/mass;
	double two_me_gammabeta2=2.*me*gammabeta2;

	double Tmax = two_me_gammabeta2/(1.0+2.0*gamma*m_ratio+m_ratio*m_ratio);

	// Density effect
	double delta=0.;	
	double X=log10(gammabeta);
	double X_0,X_1;
	double Cbar=2.*(LogI-log(28.816e-9*sqrt(rhoZ_overA)))+1.;
	if (rhoZ_overA>0.01){ // not a gas
	  if (LogI<-1.6118){ // I<100
	    if (Cbar<=3.681) X_0=0.2;
	    else X_0=0.326*Cbar-1.;
	    X_1=2.;
	  }
	  else{
	    if (Cbar<=5.215) X_0=0.2;
	    else X_0=0.326*Cbar-1.5;
	    X_1=3.;
	  }
	}
	else { // gases
	  X_1=4.;
	  if (Cbar<=9.5) X_0=1.6;
	  else if (Cbar>9.5 && Cbar<=10.) X_0=1.7;
	  else if (Cbar>10 && Cbar<=10.5) X_0=1.8;    
	  else if (Cbar>10.5 && Cbar<=11.) X_0=1.9;
	  else if (Cbar>11.0 && Cbar<=12.25) X_0=2.;
	  else if (Cbar>12.25 && Cbar<=13.804){
	    X_0=2.;
	    X_1=5.;
	  }
	  else {
	    X_0=0.326*Cbar-2.5;
	    X_1=5.;
	  } 
	}
	if (X>=X_0 && X<X_1)
	  delta=4.606*X-Cbar+(Cbar-4.606*X_0)*pow((X_1-X)/(X_1-X_0),3.);
	else if (X>=X_1)
	  delta= 4.606*X-Cbar;  	

	double dEdx = KrhoZ_overA/beta2*(log(two_me_gammabeta2*Tmax) 
					 -2.*LogI - 2.0*beta2 -delta);

	dP_dx = dEdx/beta;

	// Empirical correction for really low momentum particles. We
	// hardcode the ugly numbers just to save CPU cycles. The
	// more readable formulas are:
	//
	//  g = 0.350/sqrt(-log(0.06)) = 0.31663529 = 1/3.15820764
	//
	//  dP_dx *= 1.0 + exp(-pow(ptot/g,2.0))
	//
	// Note that for a momentum of greater than 450MeV/c this correction
	// is less than 1 percent so we apply that threshold here
	if(P<0.450){
		double a = P/0.31663529;
		dP_dx *= 1.0 + exp(-a*a);
	}
	return err;
}

