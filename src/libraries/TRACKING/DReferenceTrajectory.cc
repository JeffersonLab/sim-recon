// $Id$
//
//    File: DReferenceTrajectory.cc
// Created: Wed Jul 19 13:42:58 EDT 2006
// Creator: davidl (on Darwin swire-b241.jlab.org 8.7.0 powerpc)
//

#include <signal.h>
#include <memory>
#include <cmath>

#include <DVector3.h>
using namespace std;
#include <math.h>
#include <algorithm>

#include "DReferenceTrajectory.h"
#include "DTrackCandidate.h"
#include "DMagneticFieldStepper.h"
#include "HDGEOMETRY/DRootGeom.h"
#define ONE_THIRD 0.33333333333333333
#define TWO_THIRD 0.66666666666666667
#define EPS 1e-8
#define NaN std::numeric_limits<double>::quiet_NaN()

struct StepStruct {DReferenceTrajectory::swim_step_t steps[256];};

//---------------------------------
// DReferenceTrajectory    (Constructor)
//---------------------------------
DReferenceTrajectory::DReferenceTrajectory(const DMagneticFieldMap *bfield
														, double q
														, swim_step_t *swim_steps
														, int max_swim_steps
														, double step_size)
{
	// Copy some values into data members
	this->q = q;
	this->step_size = step_size;
	this->bfield = bfield;
	this->Nswim_steps = 0;
	this->dist_to_rt_depth = 0;
	this->mass = 0.13957; // assume pion mass until otherwise specified
	this->mass_sq=this->mass*this->mass;
	this->hit_cdc_endplate = false;
	this->RootGeom=NULL;
	this->geom = NULL;
	this->ploss_direction = kForward;
	this->check_material_boundaries = true;
	this->zmin_track_boundary = -100.0;  // boundary at which to stop swimming
	this->zmax_track_boundary = 670.0;   // boundary at which to stop swimming
	this->Rsqmax_interior = 65.0*65.0; // Maximum radius (in cm) corresponding to inside of BCAL
	this->Rsqmax_exterior = 88.0*88.0; // Maximum radius (in cm) corresponding to outside of BCAL
	
	this->last_phi = 0.0;
	this->last_swim_step = NULL;
	this->last_dist_along_wire = 0.0;
	this->last_dz_dphi = 0.0;
	
	this->debug_level = 0;
	
	// Initialize some values from configuration parameters
	BOUNDARY_STEP_FRACTION = 0.80;
	MIN_STEP_SIZE = 0.1;	// cm
	MAX_STEP_SIZE = 3.0;		// cm
	int MAX_SWIM_STEPS = 2500;
	
	gPARMS->SetDefaultParameter("TRK:BOUNDARY_STEP_FRACTION" , BOUNDARY_STEP_FRACTION, "Fraction of estimated distance to boundary to use as step size");
	gPARMS->SetDefaultParameter("TRK:MIN_STEP_SIZE" , MIN_STEP_SIZE, "Minimum step size in cm to take when swimming a track with adaptive step sizes");
	gPARMS->SetDefaultParameter("TRK:MAX_STEP_SIZE" , MAX_STEP_SIZE, "Maximum step size in cm to take when swimming a track with adaptive step sizes");
	gPARMS->SetDefaultParameter("TRK:MAX_SWIM_STEPS" , MAX_SWIM_STEPS, "Number of swim steps for DReferenceTrajectory to allocate memory for (when not using external buffer)");

	// It turns out that the greatest bottleneck in speed here comes from
	// allocating/deallocating the large block of memory required to hold
	// all of the trajectory info. The preferred way of calling this is 
	// with a pointer allocated once at program startup. This code block
	// though allows it to be allocated here if necessary.
	if(!swim_steps){
		own_swim_steps = true;
		this->max_swim_steps = MAX_SWIM_STEPS;
		this->swim_steps = new swim_step_t[this->max_swim_steps];
	}else{
		own_swim_steps = false;
		this->max_swim_steps = max_swim_steps;
		this->swim_steps = swim_steps;
	}
}

//---------------------------------
// DReferenceTrajectory    (Copy Constructor)
//---------------------------------
DReferenceTrajectory::DReferenceTrajectory(const DReferenceTrajectory& rt)
{
	/// The copy constructor will always allocate its own memory for the
	/// swim steps and set its internal flag to indicate that is owns them
	/// regardless of the owner of the source trajectory's.

	this->Nswim_steps = rt.Nswim_steps;
	this->q = rt.q;
	this->max_swim_steps = rt.max_swim_steps;
	this->own_swim_steps = true;
	this->step_size = rt.step_size;
	this->bfield = rt.bfield;
	this->last_phi = rt.last_phi;
	this->last_dist_along_wire = rt.last_dist_along_wire;
	this->last_dz_dphi = rt.last_dz_dphi;
	this->RootGeom = rt.RootGeom;
	this->geom = rt.geom;
	this->dist_to_rt_depth = 0;
	this->mass = rt.GetMass();
	this->mass_sq=this->mass*this->mass;
	this->ploss_direction = rt.ploss_direction;
	this->check_material_boundaries = rt.GetCheckMaterialBoundaries();
	this->BOUNDARY_STEP_FRACTION = rt.GetBoundaryStepFraction();
	this->MIN_STEP_SIZE = rt.GetMinStepSize();
	this->MAX_STEP_SIZE = rt.GetMaxStepSize();
	this->debug_level=rt.debug_level;
	this->zmin_track_boundary = -100.0;  // boundary at which to stop swimming
	this->zmax_track_boundary = 670.0;   // boundary at which to stop swimming
	this->Rsqmax_interior = 65.0*65.0; // Maximum radius (in cm) corresponding to inside of BCAL
	this->Rsqmax_exterior = 88.0*88.0; // Maximum radius (in cm) corresponding to outside of BCAL
	

	this->swim_steps = new swim_step_t[this->max_swim_steps];
	this->last_swim_step = NULL;
	for(int i=0; i<Nswim_steps; i++)
	{
		swim_steps[i] = rt.swim_steps[i];
		if(&(rt.swim_steps[i]) == rt.last_swim_step)
			this->last_swim_step = &(swim_steps[i]);
	}

}

//---------------------------------
// operator=               (Assignment operator)
//---------------------------------
DReferenceTrajectory& DReferenceTrajectory::operator=(const DReferenceTrajectory& rt)
{
	/// The assignment operator will always make sure the memory allocated
	/// for the swim_steps is owned by the object being copied into.
	/// If it already owns memory of sufficient size, then it will be
	/// reused. If it owns memory that is too small, it will be freed and
	/// a new block allocated. If it does not own its swim_steps coming
	/// in, then it will allocate memory so that it does own it on the
	/// way out.
	
	if(&rt == this)return *this; // protect against self copies

	// Free memory if block is too small
	if(own_swim_steps==true && max_swim_steps<rt.Nswim_steps){
		delete[] swim_steps;
		swim_steps=NULL;
	}
	
	// Forget memory block if we don't currently own it
	if(!own_swim_steps){
		swim_steps=NULL;
	}

	this->Nswim_steps = rt.Nswim_steps;
	this->q = rt.q;
	this->max_swim_steps = rt.max_swim_steps;
	this->own_swim_steps = true;
	this->step_size = rt.step_size;
	this->bfield = rt.bfield;
	this->last_phi = rt.last_phi;
	this->last_dist_along_wire = rt.last_dist_along_wire;
	this->last_dz_dphi = rt.last_dz_dphi;
	this->RootGeom = rt.RootGeom;
	this->geom = rt.geom;
	this->dist_to_rt_depth = rt.dist_to_rt_depth;
	this->mass = rt.GetMass();
	this->mass_sq=this->mass*this->mass;
	this->ploss_direction = rt.ploss_direction;
	this->check_material_boundaries = rt.GetCheckMaterialBoundaries();
	this->BOUNDARY_STEP_FRACTION = rt.GetBoundaryStepFraction();
	this->MIN_STEP_SIZE = rt.GetMinStepSize();
	this->MAX_STEP_SIZE = rt.GetMaxStepSize();

	// Allocate memory if needed
	if(swim_steps==NULL)this->swim_steps = new swim_step_t[this->max_swim_steps];

	// Copy swim steps
	this->last_swim_step = NULL;
	for(int i=0; i<Nswim_steps; i++)
	{
		swim_steps[i] = rt.swim_steps[i];
		if(&(rt.swim_steps[i]) == rt.last_swim_step)
			this->last_swim_step = &(swim_steps[i]);
	}

	
	return *this;
}

//---------------------------------
// ~DReferenceTrajectory    (Destructor)
//---------------------------------
DReferenceTrajectory::~DReferenceTrajectory()
{
	if(own_swim_steps){
		delete[] swim_steps;
	}
}

//---------------------------------
// CopyWithShift
//---------------------------------
void DReferenceTrajectory::CopyWithShift(const DReferenceTrajectory *rt, DVector3 shift)
{
	// First, do a straight copy
	*this = *rt;
	
	// Second, shift all positions
	for(int i=0; i<Nswim_steps; i++)swim_steps[i].origin += shift;
}


//---------------------------------
// Reset
//---------------------------------
void DReferenceTrajectory::Reset(void){
	//reset DReferenceTrajectory for re-use
	this->Nswim_steps = 0;
	this->ploss_direction = kForward;
	this->mass = 0.13957; // assume pion mass until otherwise specified
	this->mass_sq=this->mass*this->mass;
	this->hit_cdc_endplate = false;
	this->last_phi = 0.0;
	this->last_swim_step = NULL;
	this->last_dist_along_wire = 0.0;
	this->last_dz_dphi = 0.0;
	this->dist_to_rt_depth = 0;
	this->check_material_boundaries = true;
}

//---------------------------------
// FastSwim -- light-weight swim to a wire that does not treat multiple 
// scattering but does handle energy loss.
// No checks for distance to boundaries are done.
//---------------------------------
void DReferenceTrajectory::FastSwim(const DVector3 &pos, const DVector3 &mom,
				    DVector3 &last_pos,DVector3 &last_mom,
				    double q,double smax,
				    const DCoordinateSystem *wire){
  DVector3 mypos(pos);
  DVector3 mymom(mom);

  // Initialize the stepper
  DMagneticFieldStepper stepper(bfield, q, &pos, &mom);
  double s=0,doca=1000.,old_doca=1000.,dP_dx=0.;
  double mass=GetMass();
  while (s<smax){
    // Save old value of doca
    old_doca=doca;

    // Adjust step size to take smaller steps in regions of high momentum loss
    if(mass>0. && step_size<0.0 && geom){	
      double KrhoZ_overA=0.0;
      double rhoZ_overA=0.0;
      double LogI=0.0;
      double X0=0.0;
      if (geom->FindMatALT1(mypos,mymom,KrhoZ_overA,rhoZ_overA,LogI,X0)
	  ==NOERROR){ 
	// Calculate momentum loss due to ionization
	dP_dx = dPdx(mymom.Mag(), KrhoZ_overA, rhoZ_overA,LogI);
	double my_step_size = 0.0001/fabs(dP_dx);
		
	if(my_step_size>MAX_STEP_SIZE)my_step_size=MAX_STEP_SIZE; // maximum step size in cm
	if(my_step_size<MIN_STEP_SIZE)my_step_size=MIN_STEP_SIZE; // minimum step size in cm

	stepper.SetStepSize(my_step_size);
      }
    }
    // Swim to next
    double ds=stepper.Step(NULL);
    s+=ds;

    stepper.GetPosMom(mypos,mymom);
    if (mass>0 && dP_dx<0.){
      double ptot=mymom.Mag();
      if (ploss_direction==kForward) ptot+=dP_dx*ds;
      else ptot-=dP_dx*ds;
      mymom.SetMag(ptot);
      stepper.SetStartingParams(q, &mypos, &mymom);
    }
    
    // Break if we have passed the wire
    DVector3 wirepos=wire->origin;
    if (fabs(wire->udir.z())>0.){ // for CDC wires
      wirepos+=((mypos.z()-wire->origin.z())/wire->udir.z())*wire->udir;
    }
    doca=(wirepos-mypos).Mag();
    if (doca>old_doca) break;

    // Store the position and momentum for this step
    last_pos=mypos;
    last_mom=mymom;
  }  
}

// Faster version of the swimmer that uses an alternate stepper and does not
// check for material boundaries.
void DReferenceTrajectory::FastSwim(const DVector3 &pos, const DVector3 &mom, double q,double smax, double zmin,double zmax){
  
  /// (Re)Swim the trajectory starting from pos with momentum mom.
  /// This will use the charge and step size (if given) passed to
  /// the constructor when the object was created. It will also
  /// (re)use the swim_step buffer, replacing it's contents.
  
  // If the charged passed to us is greater that 10, it means use the charge
  // already stored in the class. Otherwise, use what was passed to us.
  if(fabs(q)>10)
    q = this->q;
  else
    this->q = q;
  
  DMagneticFieldStepper stepper(bfield, q, &pos, &mom);
  if(step_size>0.0)stepper.SetStepSize(step_size);

  // Step until we hit a boundary (don't track more than 20 meters)
  swim_step_t *swim_step = this->swim_steps;
  double t=0.;
  Nswim_steps = 0;
  double itheta02 = 0.0;
  double itheta02s = 0.0;
  double itheta02s2 = 0.0;
  double X0sum=0.0;
  swim_step_t *last_step=NULL;
  double old_radius_sq=1e6;

  // Variables used to tag the step at which the track passes into one one of
  // the outer detectors
  index_at_bcal=-1;
  index_at_tof=-1;
  index_at_fcal=-1;
  bool hit_bcal=false,hit_fcal=false,hit_tof=false;
	
  for(double s=0; fabs(s)<smax; Nswim_steps++, swim_step++){
       
    if(Nswim_steps>=this->max_swim_steps){
      if (debug_level>0){
	jerr<<__FILE__<<":"<<__LINE__<<" Too many steps in trajectory. Truncating..."<<endl;
      }
      break;
    }
    
    stepper.GetDirs(swim_step->sdir, swim_step->tdir, swim_step->udir);
    stepper.GetPosMom(swim_step->origin, swim_step->mom);
    swim_step->Ro = stepper.GetRo();
    swim_step->s = s;
    swim_step->t = t;
    
    // Magnetic field at current position
    bfield->GetField(swim_step->origin,swim_step->B);

    //magnitude of momentum and beta
    double p_sq=swim_step->mom.Mag2();
    double one_over_beta_sq=1.+mass_sq/p_sq;

    // Add material if geom or RootGeom is not NULL
    // If both are non-NULL, then use RootGeom
    double dP = 0.0;
    double dP_dx=0.0;
    if(RootGeom || geom){
      double KrhoZ_overA=0.0;
      double rhoZ_overA=0.0;
      double LogI=0.0;
      double X0=0.0;
      jerror_t err;
      if(RootGeom){
	double rhoZ_overA,rhoZ_overA_logI;
	err = RootGeom->FindMatLL(swim_step->origin,
				  rhoZ_overA, 
				  rhoZ_overA_logI, 
				  X0);
	KrhoZ_overA=0.1535e-3*rhoZ_overA;
	LogI=rhoZ_overA_logI/rhoZ_overA;
      }else{
	err = geom->FindMatALT1(swim_step->origin, swim_step->mom, KrhoZ_overA, rhoZ_overA,LogI, X0);
      }
      if(err == NOERROR){
	if(X0>0.0){
	  double p=sqrt(p_sq);
	  double delta_s = s;
	  if(last_step)delta_s -= last_step->s;
	  double radlen = delta_s/X0;
	  
	  if(radlen>1.0E-5){ // PDG 2008 pg 271, second to last paragraph
			      
	    //  double theta0 = 0.0136*sqrt(one_over_beta_sq)/p*sqrt(radlen)*(1.0+0.038*log(radlen)); // From PDG 2008 eq 27.12
	    //double theta02 = theta0*theta0;
	    double factor=1.0+0.038*log(radlen);
	    double theta02=1.8496e-4*factor*factor*radlen*one_over_beta_sq/p_sq;

	    itheta02 += theta02;
	    itheta02s += s*theta02;
	    itheta02s2 += s*s*theta02;
	    X0sum+=X0;
	  }
	  
	  // Calculate momentum loss due to ionization
	  dP_dx = dPdx(p, KrhoZ_overA, rhoZ_overA,LogI);
	}
      }
      last_step = swim_step;
    }
    swim_step->itheta02 = itheta02;
    swim_step->itheta02s = itheta02s;
    swim_step->itheta02s2 = itheta02s2;
    swim_step->invX0=Nswim_steps/X0sum;
    
    
    if(step_size<0.0){ // step_size<0 indicates auto-calculated step size
      // Adjust step size to take smaller steps in regions of high momentum loss
      double my_step_size = 0.0001/fabs(dP_dx);
      if(my_step_size>MAX_STEP_SIZE)my_step_size=MAX_STEP_SIZE; // maximum step size in cm
      if(my_step_size<MIN_STEP_SIZE)my_step_size=MIN_STEP_SIZE; // minimum step size in cm
      
      stepper.SetStepSize(my_step_size);
    }

    // Swim to next
    double ds=stepper.FastStep(swim_step->B);
  
    // Calculate momentum loss due to the step we're about to take
    dP = ds*dP_dx;
  
    // Adjust momentum due to ionization losses
    if(dP!=0.0){
      DVector3 pos, mom;
      stepper.GetPosMom(pos, mom);
      double ptot = mom.Mag() - dP; // correct for energy loss
      if (ptot<0) {Nswim_steps++; break;}
      mom.SetMag(ptot);
      stepper.SetStartingParams(q, &pos, &mom);
    }
		
    // update flight time
    t+=ds*sqrt(one_over_beta_sq)/SPEED_OF_LIGHT;
    s += ds;
  
    // Mark places along trajectory where we pass into one of the 
    // main detectors
    double Rsq=swim_step->origin.Perp2();
    double z=swim_step->origin.Z();
    if (hit_bcal==false && Rsq>Rsqmax_interior && z<407 &&z>0){
      index_at_bcal=Nswim_steps-1;
      hit_bcal=true;
    }
    if (hit_tof==false && z>606.){
      index_at_tof=Nswim_steps-1;
      hit_tof=true;
    }
    if (hit_fcal==false && z>625.){
      index_at_fcal=Nswim_steps-1;
      hit_fcal=true;     
    }
    
    // Exit the loop if we are already inside the volume of the BCAL
    // and the radius is decreasing
    if (Rsq<old_radius_sq && Rsq>Rsqmax_interior && z<407.0 && z>-100.0){
      Nswim_steps++; break;
    }
			
    // Exit loop if we leave the tracking volume
    if (z>zmax){Nswim_steps++; break;} 
    if(Rsq>Rsqmax_exterior && z<407.0){Nswim_steps++; break;} // ran into BCAL
    if (fabs(swim_step->origin.X())>160.0
	|| fabs(swim_step->origin.Y())>160.0)
      {Nswim_steps++; break;} // left extent of TOF + 31cm Buffer
    if(z>670.0){Nswim_steps++; break;} // ran into FCAL
    if(z<zmin){Nswim_steps++; break;} // exit upstream
    
    old_radius_sq=Rsq;
  }
  
  // OK. At this point the positions of the trajectory in the lab
  // frame have been recorded along with the momentum of the
  // particle and the directions of reference trajectory
  // coordinate system at each point.
}

//---------------------------------
// Swim
//---------------------------------
void DReferenceTrajectory::Swim(const DVector3 &pos, const DVector3 &mom, double q, const TMatrixFSym *cov,double smax, const DCoordinateSystem *wire)
{
        /// (Re)Swim the trajectory starting from pos with momentum mom.
	/// This will use the charge and step size (if given) passed to
	/// the constructor when the object was created. It will also
	/// (re)use the sim_step buffer, replacing it's contents.

	// If the charged passed to us is greater that 10, it means use the charge
	// already stored in the class. Otherwise, use what was passed to us.
	if(fabs(q)>10)
		q = this->q;
	else
		this->q = q;

	DMagneticFieldStepper stepper(bfield, q, &pos, &mom);
	if(step_size>0.0)stepper.SetStepSize(step_size);

	// Step until we hit a boundary (don't track more than 20 meters)
	swim_step_t *swim_step = this->swim_steps;
	double t=0.;
	Nswim_steps = 0;
	double itheta02 = 0.0;
	double itheta02s = 0.0;
	double itheta02s2 = 0.0;
	double X0sum=0.0;
	swim_step_t *last_step=NULL;
	double old_radius_sq=1e6;
	
	TMatrixFSym mycov(7);
	if (cov!=NULL){
	  mycov=*cov;
	}

	// Reset flag indicating whether we hit the CDC endplate
	// and get the parameters of the endplate so we can check
	// if we hit it while swimming.
	//hit_cdc_endplate = false;
	/*
#if 0 // The GetCDCEndplate call goes all the way back to the XML and slows down
      // overall tracking by a factor of 20. Therefore, we skip finding it
	  // and just hard-code the values instead.  1/28/2011  DL
	double cdc_endplate_z=150+17;  // roughly, from memory
	double cdc_endplate_dz=5.0;	   // roughly, from memory
	double cdc_endplate_rmin=10.0; // roughly, from memory
	double cdc_endplate_rmax=55.0; // roughly, from memory
	if(geom)geom->GetCDCEndplate(cdc_endplate_z, cdc_endplate_dz, cdc_endplate_rmin, cdc_endplate_rmax);
	double cdc_endplate_zmin = cdc_endplate_z - cdc_endplate_dz/2.0;
	double cdc_endplate_zmax = cdc_endplate_zmin + cdc_endplate_dz;
#else
	double cdc_endplate_rmin=10.0; // roughly, from memory
	double cdc_endplate_rmax=55.0; // roughly, from memory
	double cdc_endplate_zmin = 167.6;
	double cdc_endplate_zmax = 168.2;
#endif	
	*/

#if 0
	// Get Bfield from stepper to initialize Bz_old
	DVector3 B;
	stepper.GetBField(B);
	double Bz_old = B.z();
#endif

	// Variables used to tag the step at which the track passes into one 
	// one of the outer detectors
	index_at_bcal=-1;
	index_at_tof=-1;
	index_at_fcal=-1;
	bool hit_bcal=false,hit_fcal=false,hit_tof=false;
	
	for(double s=0; fabs(s)<smax; Nswim_steps++, swim_step++){
	
		if(Nswim_steps>=this->max_swim_steps){
		  if (debug_level>0){
			jerr<<__FILE__<<":"<<__LINE__<<" Too many steps in trajectory. Truncating..."<<endl;
		  }
			break;
		}

		stepper.GetDirs(swim_step->sdir, swim_step->tdir, swim_step->udir);
		stepper.GetPosMom(swim_step->origin, swim_step->mom);
		swim_step->Ro = stepper.GetRo();
		swim_step->s = s;
		swim_step->t = t;
	
		//magnitude of momentum and beta
		double p_sq=swim_step->mom.Mag2();
		double one_over_beta_sq=1.+mass_sq/p_sq;

		// Add material if geom or RootGeom is not NULL
		// If both are non-NULL, then use RootGeom
		double dP = 0.0;
		double dP_dx=0.0;
		double s_to_boundary=1.0E6; // initialize to "infinity" in case we don't set this below
		if(RootGeom || geom){
			double KrhoZ_overA=0.0;
			double rhoZ_overA=0.0;
			double LogI=0.0;
			double X0=0.0;
			jerror_t err;
			if(RootGeom){
			  double rhoZ_overA,rhoZ_overA_logI;
			  err = RootGeom->FindMatLL(swim_step->origin,
						    rhoZ_overA, 
						    rhoZ_overA_logI, 
						    X0);
			  KrhoZ_overA=0.1535e-3*rhoZ_overA;
			  LogI=rhoZ_overA_logI/rhoZ_overA;
			}else{
				if(check_material_boundaries){
					err = geom->FindMatALT1(swim_step->origin, swim_step->mom, KrhoZ_overA, rhoZ_overA,LogI, X0, &s_to_boundary);
				}else{
					err = geom->FindMatALT1(swim_step->origin, swim_step->mom, KrhoZ_overA, rhoZ_overA,LogI, X0);
				}
				
				// Check if we hit the CDC endplate
				//double z = swim_step->origin.Z();
				//if(z>=cdc_endplate_zmin && z<=cdc_endplate_zmax){
				// double r = swim_step->origin.Perp();
				// if(r>=cdc_endplate_rmin && r<=cdc_endplate_rmax){
				// hit_cdc_endplate = true;
				//}
				//}
			}

			if(err == NOERROR){
			  if(X0>0.0){
			    double p=sqrt(p_sq);
			    double delta_s = s;
			    if(last_step)delta_s -= last_step->s;
			    double radlen = delta_s/X0;
			    
			    if(radlen>1.0E-5){ // PDG 2008 pg 271, second to last paragraph
			      
			      //  double theta0 = 0.0136*sqrt(one_over_beta_sq)/p*sqrt(radlen)*(1.0+0.038*log(radlen)); // From PDG 2008 eq 27.12
			      //double theta02 = theta0*theta0;
			      double factor=1.0+0.038*log(radlen);
			      double theta02=1.8496e-4*factor*factor*radlen*one_over_beta_sq/p_sq;
			      
			      itheta02 += theta02;
			      itheta02s += s*theta02;
			      itheta02s2 += s*s*theta02;
			      X0sum+=X0;
			     
			      if (cov){
				
			      }
			    }

			    // Calculate momentum loss due to ionization
			    dP_dx = dPdx(p, KrhoZ_overA, rhoZ_overA,LogI);
			  }
			}
			last_step = swim_step;
		}
		swim_step->itheta02 = itheta02;
		swim_step->itheta02s = itheta02s;
		swim_step->itheta02s2 = itheta02s2;
		swim_step->invX0=Nswim_steps/X0sum;

		// Adjust step size to take smaller steps in regions of high momentum loss or field gradient
		if(step_size<0.0){ // step_size<0 indicates auto-calculated step size
			// Take step so as to change momentum by 100keV
			//double my_step_size=p/fabs(dP_dx)*0.01;
			double my_step_size = 0.0001/fabs(dP_dx);

			// Now check the field gradient
#if 0
			stepper.GetBField(B);
			double Bz = B.z();
			if (fabs(Bz-Bz_old)>EPS){
			  double my_step_size_B=0.01*my_step_size
			    *fabs(Bz/(Bz_old-Bz));
			  if (my_step_size_B<my_step_size) 
			    my_step_size=my_step_size_B;
			}
			Bz_old=Bz; // Save old z-component of B-field
#endif
			// Use the estimated distance to the boundary to make sure we don't overstep
			// into a high density region and miss some material. Use half the estimated
			// distance since it's only an estimate. Note that even though this would lead
			// to infinitely small steps, there is a minimum step size imposed below to
			// ensure the step size is reasonable.
			/*
			double step_size_to_boundary = BOUNDARY_STEP_FRACTION*s_to_boundary;
			if(step_size_to_boundary < my_step_size)my_step_size = step_size_to_boundary;
			*/

			if(my_step_size>MAX_STEP_SIZE)my_step_size=MAX_STEP_SIZE; // maximum step size in cm
			if(my_step_size<MIN_STEP_SIZE)my_step_size=MIN_STEP_SIZE; // minimum step size in cm

			stepper.SetStepSize(my_step_size);
		}

		// Swim to next
		double ds=stepper.Step(NULL,&swim_step->B);
		if (cov){
		  PropagateCovariance(ds,q,mass_sq,mom,pos,swim_step->B,mycov);
		  swim_step->cov_t_t=mycov(6,6);
		  swim_step->cov_px_t=mycov(6,0);
		  swim_step->cov_py_t=mycov(6,1);
		  swim_step->cov_pz_t=mycov(6,2);
		}

		// Calculate momentum loss due to the step we're about to take
		dP = ds*dP_dx;

		// Adjust momentum due to ionization losses
		if(dP!=0.0){
			DVector3 pos, mom;
			stepper.GetPosMom(pos, mom);
			double ptot = mom.Mag() - dP; // correct for energy loss
			bool ranged_out = false;
			/*
			if (ptot<0.05){
			  swim_step->origin.Print();
			  cout<<"N: " << Nswim_steps <<" x " << pos.x() <<" y " <<pos.y() <<" z " << pos.z() <<" r " << pos.Perp()<< " s " << s  << " p " << ptot << endl;
			}
			*/
			if(ptot<0.0)ranged_out=true;
			if(dP<0.0 && ploss_direction==kForward)ranged_out=true;
			if(dP>0.0 && ploss_direction==kBackward)ranged_out=true;
			if(mom.Mag()==0.0)ranged_out=true;
			if(ranged_out){
				Nswim_steps++; // This will at least allow for very low momentum particles to have 1 swim step
				break;
			}
			mom.SetMag(ptot);
			stepper.SetStartingParams(q, &pos, &mom);
		}
		
		// update flight time
		t+=ds*sqrt(one_over_beta_sq)/SPEED_OF_LIGHT;
		s += ds;

		// Mark places along trajectory where we pass into one of the 
		// main detectors
		double Rsq=swim_step->origin.Perp2();
		double z=swim_step->origin.Z();
		if (hit_bcal==false && Rsq>Rsqmax_interior && z<407 &&z>0){
		  index_at_bcal=Nswim_steps-1;
		  hit_bcal=true;
		}
		if (hit_tof==false && z>618.){
		  index_at_tof=Nswim_steps-1;
		  hit_tof=true;
		}
		if (hit_fcal==false && z>625.){
		  index_at_fcal=Nswim_steps-1;
		  hit_fcal=true;
		}

		// Exit the loop if we are already inside the volume of the BCAL
		// and the radius is decreasing
		if (Rsq<old_radius_sq && Rsq>Rsqmax_interior && z<407.0 && z>-100.0){
		  Nswim_steps++; break;
		}
	
		
		// Exit loop if we leave the tracking volume
		if(Rsq>Rsqmax_exterior && z<407.0){Nswim_steps++; break;} // ran into BCAL
		if (fabs(swim_step->origin.X())>160.0
		    || fabs(swim_step->origin.Y())>160.0)
		  {Nswim_steps++; break;} // left extent of TOF + 31cm Buffer
		if(z>zmax_track_boundary){Nswim_steps++; break;} // ran into FCAL
		if(z<zmin_track_boundary){Nswim_steps++; break;} // exit upstream
		if(wire && Nswim_steps>0){ // optionally check if we passed a wire we're supposed to be swimming to
			swim_step_t *closest_step = FindClosestSwimStep(wire);
			if(++closest_step!=swim_step){Nswim_steps++; break;}
		}

		old_radius_sq=Rsq;
	}

	// OK. At this point the positions of the trajectory in the lab
	// frame have been recorded along with the momentum of the
	// particle and the directions of reference trajectory
	// coordinate system at each point.
}

// Routine to find position on the trajectory where the track crosses a radial
// position R.  Also returns the path length to this position.
jerror_t DReferenceTrajectory::GetIntersectionWithRadius(double R,
							 DVector3 &mypos,
							 double *s,
							 double *t,
							 DVector3 *p_at_intersection) const{
  mypos.SetXYZ(NaN,NaN,NaN);
  if(p_at_intersection)
    p_at_intersection->SetXYZ(NaN,NaN,NaN);

  if(Nswim_steps<1){
    _DBG_<<"No swim steps! You must \"Swim\" the track before calling GetIntersectionWithRadius(...)"<<endl;
  }
  // Return early if the radius at the end of the reference trajectory is still less than R
  double outer_radius=swim_steps[Nswim_steps-1].origin.Perp();
  if (outer_radius<R){
    if (s) *s=0.;
    if (t) *t=0.;
    return VALUE_OUT_OF_RANGE;
  }
  // Return early if the radius at the beginning of the trajectory is outside
  // the radius to which we are trying to match
  double inner_radius=swim_steps[0].origin.Perp();
  if (inner_radius>R){
    if (s) *s=0.;
    if (t) *t=0.;
    return VALUE_OUT_OF_RANGE;
  }


  // Loop over swim steps and find the one that crosses the radius
  swim_step_t *swim_step = swim_steps;
  swim_step_t *step=NULL;
  swim_step_t *last_step=NULL;

  //  double inner_radius=swim_step->origin.Perp();
  for(int i=0; i<Nswim_steps; i++, swim_step++){
    if (swim_step->origin.Perp()>R){
      step=swim_step;
      break;
    }
    if (swim_step->origin.Z()>407.0) return VALUE_OUT_OF_RANGE;
    last_step=swim_step;
  }
  if (step==NULL||last_step==NULL) return VALUE_OUT_OF_RANGE;
  if (p_at_intersection!=NULL){
    *p_at_intersection=last_step->mom;
  }

  // At this point, the location where the track intersects the cyclinder 
  // is somewhere between last_step and step. For simplicity, we're going
  // to just find the intersection of the cylinder with the line that joins
  // the 2 positions. We do this by working in the X/Y plane only and
  // finding the value of "alpha" which is the fractional distance the
  // intersection point is between last_pos and mypos. We'll then apply
  // the alpha found in the 2D X/Y space to the 3D x/y/Z space to find
  // the actual intersection point.
  DVector2 x1(last_step->origin.X(), last_step->origin.Y());
  DVector2 x2(step->origin.X(), step->origin.Y());
  DVector2 dx = x2-x1;
  double A = dx.Mod2();
  double B = 2.0*(x1.X()*dx.X() + x1.Y()*dx.Y());
  double C = x1.Mod2() - R*R;
  
  double sqrt_D=sqrt(B*B-4.0*A*C);
  double one_over_denom=0.5/A;
  double alpha1 = (-B + sqrt_D)*one_over_denom;
  double alpha2 = (-B - sqrt_D)*one_over_denom;
  double alpha = alpha1;
  if(alpha1<0.0 || alpha1>1.0)alpha=alpha2;
  if(!isfinite(alpha))return VALUE_OUT_OF_RANGE;
	
  DVector3 delta = step->origin - last_step->origin;
  mypos = last_step->origin + alpha*delta;
  
  // The value of s actually represents the pathlength
  // to the outside point. Adjust it back to the
  // intersection point (approximately).
  if (s) *s = step->s-(1.0-alpha)*delta.Mag();

  // flight time
  if (t){	
    double p_sq=step->mom.Mag2();
    double one_over_beta=sqrt(1.+mass_sq/p_sq);
    *t = step->t-(1.0-alpha)*delta.Mag()*one_over_beta/SPEED_OF_LIGHT;
  }

  return NOERROR;
}

//---------------------------------
// GetIntersectionWithPlane
//---------------------------------
jerror_t DReferenceTrajectory::GetIntersectionWithPlane(const DVector3 &origin, const DVector3 &norm, DVector3 &pos, double *s,double *t,double *var_t,DetectorSystem_t detector) const{
  DVector3 dummy;
  return GetIntersectionWithPlane(origin,norm,pos,dummy,s,t,var_t,detector);
}
jerror_t DReferenceTrajectory::GetIntersectionWithPlane(const DVector3 &origin, const DVector3 &norm, DVector3 &pos, DVector3 &p_at_intersection, double *s,
							double *t,double *var_t,DetectorSystem_t detector) const
{
	/// Get the intersection point of this trajectory with a plane.
	/// The plane is specified by <i>origin</i> and <i>norm</i>. The
	/// <i>origin</i> vector should give the coordinates of any point
	/// on the plane and <i>norm</i> should give a vector normal to
	/// the plane. The <i>norm</i> vector will be copied and normalized
	/// so it can be of any magnitude upon entry.
	///
	/// The coordinates of the intersection point will copied into
	/// the supplied <i>pos</i> vector. If a non-NULL pointer for <i>s</i>
	/// is passed in, the pathlength of the trajectory from its begining
	/// to the intersection point is copied into location pointed to.
	
	// Set reasonable defaults
	pos.SetXYZ(0,0,0);
	if(s)*s=0.0;
	if(t)*t=0.0;
	p_at_intersection.SetXYZ(0,0,0);
	
	// Return early if the z-position of the plane with which we are 
	// intersecting is beyong the reference trajectory.
	if (origin.z()>swim_steps[Nswim_steps-1].origin.z()){
	  return VALUE_OUT_OF_RANGE;
	}
	// Find the closest swim step to the position where the track crosses
	// the plane
	int first_i=0;
	switch(detector){
	case SYS_FCAL:
	  if (index_at_fcal<0) return VALUE_OUT_OF_RANGE;
	  first_i=index_at_fcal;
	  break;
	case SYS_TOF:
	  if (index_at_tof<0) return VALUE_OUT_OF_RANGE;
	  first_i=index_at_tof;
	  break;
	default:	
	  break;
	}
	swim_step_t *step=FindPlaneCrossing(origin,norm,first_i,detector);
	if(!step){
		_DBG_<<"Could not find closest swim step!"<<endl;
		return RESOURCE_UNAVAILABLE;
	}

	// Here we follow a scheme described in more detail in the 
	// DistToRT(DVector3 hit) method below. The basic idea is to
	// express a point on the helix in terms of a single variable
	// and then solve for that variable by setting the distance 
	// to zero.
	//
	//   x = Ro*(cos(phi) - 1)
	//   y = Ro*sin(phi)
	//   z = phi*(dz/dphi)
	//
	// As is done below, we work in the RT coordinate system. Well,
	// sort of. The distance to the plane is given by:
	//
	//   d = ( x - c ).n
	//
	// where x is a point on the helix, c is the "origin" point
	// which lies somewhere in the plane and n is the "norm"
	// vector. Since we want a point in the plane, we set d=0
	// and solve for phi (with the components of x expressed in
	// terms of phi as given in the DistToRT method below).
	//
	// Thus, the equation we need to solve is:
	//
	// x.n - c.n = 0
	//
	// notice that "c" only gets dotted into "n" so that
	// dot product can occur in any coordinate system. Therefore,
	// we do that in the lab coordinate system to avoid the
	// overhead of transforming "c" to the RT system. The "n"
	// vector, however, still must be transformed.
	//
	// Expanding the trig functions to 2nd order in phi, performing
	// the x.n dot product, and gathering equal powers of phi
	// leads us to he following:
	//
	//  (-Ro*nx/2)*phi^2 + (Ro*ny+dz_dphi*nz)*phi - c.n = 0
	//
	// which is quadratic in phi. We solve for both roots, but use
	// the one with the smller absolute value (if both are finite).

	double &Ro = step->Ro;

	// OK, having said all of that, it turns out that the above 
	// mechanism will tend to fail in regions of low or no
	// field because the value of Ro is very large. Thus, we need to
	// use a straight line projection in such cases. We also
	// want to use a straight line projection if the helical intersection
	// fails for some other reason.
	//
	// The algorthim is then to only try the helical calculation
	// for small (<10m) values of Ro and then do the straight line
	// if R is larger than that OR the helical calculation fails.

	// Try helical calculation
	if(Ro<1000.0){
		double nx = norm.Dot(step->sdir);
		double ny = norm.Dot(step->tdir);
		double nz = norm.Dot(step->udir);
	
		double delta_z = step->mom.Dot(step->udir);
		double delta_phi = step->mom.Dot(step->tdir)/Ro;
		double dz_dphi = delta_z/delta_phi;
		
		double A = -Ro*nx/2.0;
		double B = Ro*ny + dz_dphi*nz;
		double C = norm.Dot(step->origin-origin);
		double sqroot=sqrt(B*B-4.0*A*C);
		double twoA=2.0*A;

		double phi_1 = (-B + sqroot)/(twoA);
		double phi_2 = (-B - sqroot)/(twoA);
		
		double phi = fabs(phi_1)<fabs(phi_2) ? phi_1:phi_2;
		if(!isfinite(phi_1))phi = phi_2;
		if(!isfinite(phi_2))phi = phi_1;
		if(isfinite(phi)){
		
			double my_s = -Ro/2.0 * phi*phi;
			double my_t = Ro * phi;
			double my_u = dz_dphi * phi;
			
			pos = step->origin + my_s*step->sdir + my_t*step->tdir + my_u*step->udir;
			p_at_intersection = step->mom;
			if(s){
				double delta_s = sqrt(my_t*my_t + my_u*my_u);
				*s = step->s + (phi>0 ? +delta_s:-delta_s);
			}	  
			// flight time
			if (t){
			  double delta_s = sqrt(my_t*my_t + my_u*my_u);
			  double ds=(phi>0 ? +delta_s:-delta_s);
			  double p_sq=step->mom.Mag2();
			  double one_over_beta=sqrt(1.+mass_sq/p_sq);
			  *t = step->t+ds*one_over_beta/SPEED_OF_LIGHT;
			}
			if (var_t){
			  *var_t=step->cov_t_t;
			}
			
			// Success. Go ahead and return
			return NOERROR;
		}
	}
	
	// If we got here then we need to try a straight line calculation
	double p_sq=step->mom.Mag2();
	double dz_over_pz=(origin.z()-step->origin.z())/step->mom.z();
	double ds=sqrt(p_sq)*dz_over_pz;
	pos.SetXYZ(step->origin.x()+dz_over_pz*step->mom.x(),
		   step->origin.y()+dz_over_pz*step->mom.y(),
		   origin.z());
	p_at_intersection = step->mom;

	if(s){
		*s = step->s + ds;
	}
	// flight time
	if (t){
	  double one_over_beta=sqrt(1.+mass_sq/p_sq);
	  *t = step->t+ds*one_over_beta/SPEED_OF_LIGHT;
	}
	// Flight time variance
	if (var_t){
	  *var_t=step->cov_t_t;
	}

	return NOERROR;
}

//---------------------------------
// InsertSteps
//---------------------------------
int DReferenceTrajectory::InsertSteps(const swim_step_t *start_step, double delta_s, double step_size)
{
	/// Insert additional steps into the reference trajectory starting
	/// at start_step and swimming for at least delta_s by step_size
	/// sized steps. Both delta_s and step_size are in centimeters.
	/// If the value of delta_s is negative then the particle's momentum
	/// and charge are reversed before swimming. This could be a problem
	/// if energy loss is implemented.
	
	if(!start_step)return -1;

	// We do this by creating another, temporary DReferenceTrajectory object
	// on the stack and swimming it.
	DVector3 pos = start_step->origin;
	DVector3 mom = start_step->mom;
	double my_q = q;
	int direction = +1;
	if(delta_s<0.0){
		mom *= -1.0;
		my_q = -q;
		direction = -1;
	}
	
	// Here I allocate the steps using an auto_ptr so I don't have to mess with
	// deleting them at all of the possible exits. The problem with auto_ptr
	// is it can't handle arrays so it has to be wrapped in a struct.
	auto_ptr<StepStruct> steps_aptr(new StepStruct);
	DReferenceTrajectory::swim_step_t *steps = steps_aptr->steps;
	DReferenceTrajectory rt(bfield , my_q , steps , 256);
	rt.SetStepSize(step_size);
	rt.Swim(pos, mom, my_q,NULL,fabs(delta_s));
	if(rt.Nswim_steps==0)return 1;

	// Check that there is enough space to add these points
	if((Nswim_steps+rt.Nswim_steps)>max_swim_steps){
		//_DBG_<<"Not enough swim steps available to add new ones! Max="<<max_swim_steps<<" had="<<Nswim_steps<<" new="<<rt.Nswim_steps<<endl;
		return 2;
	}
	
	// At this point, we may have swum forward or backwards so the points
	// will need to be added either before start_step or after it. We also
	// may want to replace an old step that overlaps our high density steps
	// since they are presumably more accurate. Find the indexes of the
	// existing steps that the new steps will be inserted between.
	double sdiff = rt.swim_steps[rt.Nswim_steps-1].s;
	double s1 = start_step->s;
	double s2 = start_step->s + (double)direction*sdiff;
	double smin = s1<s2 ? s1:s2;
	double smax = s1<s2 ? s2:s1;
	int istep_start = 0;
	int istep_end = 0;
	for(int i=0; i<Nswim_steps; i++){
		if(swim_steps[i].s <  smin)istep_start = i;
		if(swim_steps[i].s <= smax)istep_end = i+1;
	}
	
	// istep_start and istep_end now point to the steps we want to keep.
	// All steps between them (exclusive) will be overwritten. Note that 
	// the original start_step should be in the "overwrite" range since 
	// it is included already in the new trajectory.
	int steps_to_overwrite = istep_end - istep_start - 1;
	int steps_to_shift = rt.Nswim_steps - steps_to_overwrite;
	
	// Shift the steps down (or is it up?) starting with istep_end.
	for(int i=Nswim_steps-1; i>=istep_end; i--)swim_steps[i+steps_to_shift] = swim_steps[i];
	
	// Copy the new steps into this object
	double s_0 = start_step->s;
	double itheta02_0 = start_step->itheta02;
	double itheta02s_0 = start_step->itheta02s;
	double itheta02s2_0 = start_step->itheta02s2;
	for(int i=0; i<rt.Nswim_steps; i++){
		int index = direction>0 ? (istep_start+1+i):(istep_start+1+rt.Nswim_steps-1-i);
		swim_steps[index] = rt.swim_steps[i];
		swim_steps[index].s = s_0 + (double)direction*swim_steps[index].s;
		swim_steps[index].itheta02 = itheta02_0 + (double)direction*swim_steps[index].itheta02;
		swim_steps[index].itheta02s = itheta02s_0 + (double)direction*swim_steps[index].itheta02s;
		swim_steps[index].itheta02s2 = itheta02s2_0 + (double)direction*swim_steps[index].itheta02s2;
		if(direction<0.0){
			swim_steps[index].sdir *= -1.0;
			swim_steps[index].tdir *= -1.0;
		}
	}
	Nswim_steps += rt.Nswim_steps-steps_to_overwrite;

	// Note that the above procedure may leave us with "kinks" in the itheta0 
	// variables. It may be that we need to recalculate those for all of the 
	// new points and the ones after them by making one more pass. I'm hoping
	// it is a realitively small correction though so we can skip it here.
	return 0;
}

//---------------------------------
// DistToRTwithTime
//---------------------------------
double DReferenceTrajectory::DistToRTwithTime(DVector3 hit, double *s,double *t,
					      double *var_t,
					      DetectorSystem_t detector) const{
  double dist=DistToRT(hit,s,detector);
  if (s!=NULL && t!=NULL) 
  {
    if(last_swim_step==NULL)
    {
      *s = 1.0E6;
      *t = 1.0E6;
      if (var_t!=NULL){
	*var_t=1.0E6;
      }
    }
    else
    {
      double p_sq=last_swim_step->mom.Mag2();
      double one_over_beta=sqrt(1.+mass_sq/p_sq);
      *t=last_swim_step->t+(*s-last_swim_step->s)*one_over_beta/SPEED_OF_LIGHT;
      if (var_t!=NULL){
	*var_t=last_swim_step->cov_t_t;
      }
    }
  }
  return dist;
}

//---------------------------------
// DistToRT
//---------------------------------
double DReferenceTrajectory::DistToRT(DVector3 hit, double *s,
				      DetectorSystem_t detector) const
{
	last_swim_step=NULL;
	if(Nswim_steps<1)_DBG__;

	int start_index=0;
	switch(detector){
	case SYS_BCAL:
	  if (index_at_bcal<0) return numeric_limits<double>::quiet_NaN();
	  start_index=index_at_bcal;
	  break;	
	case SYS_FCAL:
	  if (index_at_fcal<0) return numeric_limits<double>::quiet_NaN();
	  start_index=index_at_fcal;
	  break;
	case SYS_TOF:
	  if (index_at_tof<0) return numeric_limits<double>::quiet_NaN();
	  start_index=index_at_tof;
	  break;
	default:
	  break;
	}

	// First, find closest step to point
	swim_step_t *swim_step = &swim_steps[start_index];
	swim_step_t *step=NULL;
	
	//double min_delta2 = 1.0E6;
	double old_delta2=10.e6,delta2=1.0e6;
	
	// Check if we should start at the end of the reference trajectory 
	// or the beginning...
	int last_index=Nswim_steps-1;
	double forward_delta2=(swim_step->origin - hit).Mag2();
	double backward_delta2=(swim_steps[last_index].origin-hit).Mag2();

	if (forward_delta2<backward_delta2){ // start at the beginning
	  for(int i=start_index; i<Nswim_steps; i++, swim_step++){
	    
	    DVector3 pos_diff = swim_step->origin - hit;
	    delta2 = pos_diff.Mag2();

	    if (delta2>old_delta2){
	      break;
	    }
	    
	    //if(delta2 < min_delta2){
	    //min_delta2 = delta2;
	    
	    step = swim_step;
	    old_delta2=delta2;
	    //}
	  }
	}
	else{// start at the end
	  for(int i=last_index; i>=start_index; i--){
	    swim_step=&swim_steps[i];
	    DVector3 pos_diff = swim_step->origin - hit;
	    delta2 = pos_diff.Mag2();
	    if (delta2>old_delta2) break;
	    
	    //if(delta2 < min_delta2){
	    //min_delta2 = delta2;
	    
	    step = swim_step;
	    old_delta2=delta2;
	    //}
	  }

	}

	if(step==NULL){
		// It seems to occasionally occur that we have 1 swim step
		// and it's values are invalid. Supress warning messages
		// for these as they are "known" (even if not fully understood!)
	  if(s != NULL)
	    *s = 1.0E6;
	  if(Nswim_steps>1){
	    _DBG_<<"\"hit\" passed to DistToRT(DVector3) out of range!"<<endl;
	    _DBG_<<"hit x,y,z = "<<hit.x()<<", "<<hit.y()<<", "<<hit.z()<<"  Nswim_steps="<<Nswim_steps<<"  min_delta2="<<delta2<<endl;
	  }
	  return 1.0E6;
	}
	
	// store last step
	last_swim_step=step;

	
	// Next, define a point on the helical segment defined by the
	// swim step it the RT coordinate system. The directions of
	// the RT coordinate system are defined by step->xdir, step->ydir,
	// and step->zdir. The coordinates of a point on the helix
	// in this coordinate system are:
	//
	//   x = Ro*(cos(phi) - 1)
	//   y = Ro*sin(phi)
	//   z = phi*(dz/dphi)
	//
	// where phi is the phi angle of the point in this coordinate system.
	// phi=0 corresponds to the swim step point itself
	//
	// Transform the given coordinates to the RT coordinate system
	// and call these x0,y0,z0. Then, the distance of point to a
	// point on the helical segment is given by:
	//
	//   d^2 = (x0-x)^2 + (y0-y)^2 + (z0-z)^2
	//
	// where x,y,z are all functions of phi as given above.
	//
	// writing out d^2 in terms of phi, but using the small angle
	// approximation for the trig functions, an equation for the
	// distance in only phi is obtained. Taking the derivative 
	// and setting it equal to zero leaves a 3rd order polynomial
	// in phi whose root corresponds to the minimum distance.
	// Skipping some math, this equation has the form:
	//
	// d(d^2)/dphi = 0 = Ro^2*phi^3 + 2*alpha*phi + beta
	//
	// where:
	//       alpha = x0*Ro + Ro^2 + (dz/dphi)^2
	//
	//        beta = -2*y0*Ro - 2*z0*(dz/dphi)
	//
	// The above 3rd order poly is convenient in that it does not
	// contain a phi^2 term. This means we can skip the step
	// done in the general case where a change of variables is
	// made such that the 2nd order term disappears.
	//
	// In general, an equation of the form
	//
	//  w^3 + 3.0*b*w + 2*c = 0 
	//
	// has one real root:
	//
	//  w0 = q - p
	//
	// where:
	//    q^3 = d - c
	//    p^3 = d + c
	//    d^2 = b^3 + c^2      (don't confuse with d^2 above!)
	//
	// So for us ...
	//
	//    3b = 2*alpha/(Ro^2)
	//    2c = beta/(Ro^2)

	hit -= step->origin;
	double x0 = hit.Dot(step->sdir);
	double y0 = hit.Dot(step->tdir);
	double z0 = hit.Dot(step->udir);
	double &Ro = step->Ro;
	double Ro2 = Ro*Ro;
	double delta_z = step->mom.Dot(step->udir);
	double delta_phi = step->mom.Dot(step->tdir)/Ro;
	double dz_dphi = delta_z/delta_phi;
	
	//  double alpha = x0*Ro + Ro2 + pow(dz_dphi,2.0);
	double alpha=x0*Ro + Ro2 +dz_dphi*dz_dphi;
	//  double beta = -2.0*y0*Ro - 2.0*z0*dz_dphi;
	double beta = -2.0*(y0*Ro + z0*dz_dphi);
	//  double b = (2.0*alpha/Ro2)/3.0;
	double b= TWO_THIRD*alpha/Ro2;
	//  double c = (beta/Ro2)/2.0;
	double c = 0.5*(beta/Ro2);
	//  double d = sqrt(pow(b,3.0) + pow(c,2.0));
	double d2=b*b*b+c*c;
	double phi=0.,dist2=1e8;
	if (d2>=0){
	  double d=sqrt(d2);
	  //double q = pow(d-c, ONE_THIRD);
	  //double p = pow(d+c, ONE_THIRD);
	  double p=cbrt(d+c);
	  double q=cbrt(d-c);
	  phi = q - p;
	  if (fabs(phi)<0.2){ // check small angle approximation
	    double phisq=phi*phi;
	  
	    dist2 = 0.25*Ro2*phisq*phisq + alpha*phisq + beta*phi 
	      + x0*x0 + y0*y0 + z0*z0;
	  }
	  else{
	    return numeric_limits<double>::quiet_NaN();
	  }
	}
	else{
	  // Use DeMoivre's theorem to find the cube root of a complex
	  // number.  In this case there are three real solutions.
	  double d=sqrt(-d2);
	  c*=-1.;
	  double temp=sqrt(cbrt(c*c+d*d));
	  double theta1=ONE_THIRD*atan2(d,c);
	  double sum_over_2=temp*cos(theta1);
	  double diff_over_2=-temp*sin(theta1);

	  double phi0=2.*sum_over_2;
	  double phi0sq=phi0*phi0;
	  double phi1=-sum_over_2+sqrt(3)*diff_over_2;
	  double phi1sq=phi1*phi1;
	  double phi2=-sum_over_2-sqrt(3)*diff_over_2;
	  double phi2sq=phi2*phi2;
	  double d2_2 = 0.25*Ro2*phi2sq*phi2sq + alpha*phi2sq + beta*phi2 + x0*x0 + y0*y0 + z0*z0;
	  double d2_1 = 0.25*Ro2*phi1sq*phi1sq + alpha*phi1sq + beta*phi1 + x0*x0 + y0*y0 + z0*z0; 
	  double d2_0 = 0.25*Ro2*phi0sq*phi0sq + alpha*phi0sq + beta*phi0 + x0*x0 + y0*y0 + z0*z0;
 
	  if (d2_0<d2_1 && d2_0<d2_2){
	    phi=phi0;
	    dist2=d2_0;
	  }
	  else if (d2_1<d2_0 && d2_1<d2_2){
	    phi=phi1;
	    dist2=d2_1;
	  }
	  else{
	    phi=phi2;
	    dist2=d2_2;
	  }
	  if (fabs(phi)<0.2){ // Check small angle approximation
	    return numeric_limits<double>::quiet_NaN();
	  }

	  if (std::isnan(Ro))
	    {
	  }
	}
	
	// Calculate distance along track ("s")
	if(s!=NULL){
		double dz = dz_dphi*phi;
		double Rodphi = Ro*phi;
		double ds = sqrt(dz*dz + Rodphi*Rodphi);
		*s = step->s + (phi>0.0 ? ds:-ds);
	}

	this->last_phi = phi;
	this->last_swim_step = step;
	this->last_dz_dphi = dz_dphi;

	return sqrt(dist2);
}

//---------------------------------
// FindClosestSwimStep
//---------------------------------
DReferenceTrajectory::swim_step_t* DReferenceTrajectory::FindClosestSwimStep(const DCoordinateSystem *wire, int *istep_ptr) const
{
	/// Find the closest swim step to the given wire. The value of
	/// "L" should be the active wire length. The coordinate system
	/// defined by "wire" should have its origin at the center of
	/// the wire with the wire running in the direction of udir.
	
	if(istep_ptr)*istep_ptr=-1;
	
	if(Nswim_steps<1){
		_DBG_<<"No swim steps! You must \"Swim\" the track before calling FindClosestSwimStep(...)"<<endl;
	}

	// Make sure we have a wire first!
	if(!wire)return NULL;
	
	// Loop over swim steps and find the one closest to the wire
	swim_step_t *swim_step = swim_steps;
	swim_step_t *step=NULL;
	//double min_delta2 = 1.0E6;
	double old_delta2=1.0e6;
	double L_over_2 = wire->L/2.0; // half-length of wire in cm
	int istep=-1;

	double dx, dy, dz;
  
  // w is a vector to the origin of the wire
  // u is a unit vector along the wire
  
  double wx, wy, wz;
  double ux, uy, uz;
  
  wx = wire->origin.X();
  wy = wire->origin.Y();
  wz = wire->origin.Z();
  
  ux = wire->udir.X();
  uy = wire->udir.Y();
  uz = wire->udir.Z();
  
  int i;
  for(i=0; i<Nswim_steps; i++, swim_step++){
		// Find the point's position along the wire. If the point
		// is past the end of the wire, calculate the distance
		// from the end of the wire.
    //		DVector3 pos_diff = swim_step->origin - wire->origin;

    dx = swim_step->origin.X() - wx;
    dy = swim_step->origin.Y() - wy;
    dz = swim_step->origin.Z() - wz;
		
    //    double u = wire->udir.Dot(pos_diff);
    double u = ux * dx + uy * dy + uz * dz;

		// Find distance perpendicular to wire
    //		double delta2 = pos_diff.Mag2() - u*u;
    double delta2 = dx*dx + dy*dy + dz*dz - u*u;

		// If point is past end of wire, calculate distance
		// from wire's end by adding on distance along wire direction.
		if( fabs(u)>L_over_2){
      //			delta2 += pow(fabs(u)-L_over_2, 2.0);
		  double u_minus_L_over_2=fabs(u)-L_over_2;
      delta2 += ( u_minus_L_over_2*u_minus_L_over_2 );
      // printf("step %d\n",i);
		}

		if(debug_level>3)_DBG_<<"delta2="<<delta2<<"  old_delta2="<<old_delta2<<endl;
		if (delta2>old_delta2) break;

		//if(delta2 < min_delta2){
		//	min_delta2 = delta2;
		step = swim_step;
		istep=i;

		//}
		//printf("%d delta %f min %f\n",i,delta2,min_delta2);
		old_delta2=delta2;
	}

	if(istep_ptr)*istep_ptr=istep;
	
	if(debug_level>3)_DBG_<<"found closest step at i="<<i<<"  istep_ptr="<<istep_ptr<<endl;

	return step;	
}

//---------------------------------
// FindClosestSwimStep
//---------------------------------
DReferenceTrajectory::swim_step_t* DReferenceTrajectory::FindClosestSwimStep(const DVector3 &origin, DVector3 norm, int *istep_ptr) const
{
	/// Find the closest swim step to the plane specified by origin
	/// and norm. origin should indicate any point in the plane and
	/// norm a vector normal to the plane.
	if(istep_ptr)*istep_ptr=-1;
	
	if(Nswim_steps<1){
		_DBG_<<"No swim steps! You must \"Swim\" the track before calling FindClosestSwimStep(...)"<<endl;
	}

	// Make sure normal vector is unit lenght
	norm.SetMag(1.0);

	// Loop over swim steps and find the one closest to the plane
	swim_step_t *swim_step = swim_steps;
	swim_step_t *step=NULL;
	//double min_dist = 1.0E6;
	double old_dist=1.0e6;
	int istep=-1;

	for(int i=0; i<Nswim_steps; i++, swim_step++){
	
		// Distance to plane is dot product of normal vector with any
		// vector pointing from the current step to a point in the plane
		double dist = fabs(norm.Dot(swim_step->origin-origin));

		if (dist>old_dist) break;

		// Check if we're the closest step
		//if(dist < min_dist){
		//min_dist = dist;

		step = swim_step;
		istep=i;
			//}
		old_dist=dist;

		// We should probably have a break condition here so we don't
		// waste time looking all the way to the end of the track after
		// we've passed the plane.
	}

	if(istep_ptr)*istep_ptr=istep;

	return step;	
}


//---------------------------------
// FindPlaneCrossing
//---------------------------------
DReferenceTrajectory::swim_step_t* DReferenceTrajectory::FindPlaneCrossing(const DVector3 &origin, DVector3 norm,int first_i,DetectorSystem_t detector) const
{
  /// Find the closest swim step to the position where the track crosses
  /// the plane specified by origin
  /// and norm. origin should indicate any point in the plane and
  /// norm a vector normal to the plane.
	
	if(Nswim_steps<1){
		_DBG_<<"No swim steps! You must \"Swim\" the track before calling FindPlaneCrossing(...)"<<endl;
		raise(SIGSEGV);// force seg. fault
	}

	// Make sure normal vector is unit length
	norm.SetMag(1.0);

	// Loop over swim steps and find the one closest to the plane
	swim_step_t *swim_step = &swim_steps[first_i];
	swim_step_t *step=NULL;
	double old_dist=1.0e6;

	// Check if we should start from the beginning of the reference 
	// trajectory or the end
	int last_index=Nswim_steps-1;
	double forward_dist= norm.Dot(swim_step->origin-origin);
	if( forward_dist == 0.0 ) return swim_step;
	double backward_dist= norm.Dot(swim_steps[last_index].origin-origin);
	if( backward_dist ==0.0 ) return &swim_steps[last_index];
	if (detector==SYS_START || fabs(forward_dist)<fabs(backward_dist)){ // start at beginning
	  for(int i=first_i; i<Nswim_steps; i++, swim_step++){
	      
	    // Distance to plane is dot product of normal vector with any
	    // vector pointing from the current step to a point in the plane
	    //double dist = fabs(norm.Dot(swim_step->origin-origin));
	    double dist = norm.Dot(swim_step->origin-origin);
	      
	    // We've crossed the plane when the sign of dist changes
	    if (dist*old_dist<0 && i>0) {
	      if (fabs(dist)<fabs(old_dist)){
		step=swim_step;
	      }
	      break;
	    }
	    step = swim_step;
	    old_dist=dist;
	  }
	}
	else{ // start at end
	  for(int i=last_index; i>=0; i--){
	    swim_step=&swim_steps[i];
	    double dist = norm.Dot(swim_step->origin-origin);
	    // We've crossed the plane when the sign of dist changes
	    if (dist*old_dist<0 && i<last_index) {
	      if (fabs(dist)<fabs(old_dist)){
		step=swim_step;
	      }
	      break;
	    }
	    step = swim_step;
	    old_dist=dist;
	  }

	}

	return step;	
}




//---------------------------------
// DistToRT
//---------------------------------
double DReferenceTrajectory::DistToRT(const DCoordinateSystem *wire, double *s) const
{
	/// Find the closest distance to the given wire in cm. The value of
	/// "L" should be the active wire length (in cm). The coordinate system
	/// defined by "wire" should have its origin at the center of
	/// the wire with the wire running in the direction of udir.
	swim_step_t *step=FindClosestSwimStep(wire);

	return (step && step->s>0) ? DistToRT(wire, step, s):std::numeric_limits<double>::quiet_NaN();
}

//---------------------------------
// DistToRTBruteForce
//---------------------------------
double DReferenceTrajectory::DistToRTBruteForce(const DCoordinateSystem *wire, double *s) const
{
	/// Find the closest distance to the given wire in cm. The value of
	/// "L" should be the active wire length (in cm). The coordinate system
	/// defined by "wire" should have its origin at the center of
	/// the wire with the wire running in the direction of udir.
	swim_step_t *step=FindClosestSwimStep(wire);

	return step ? DistToRTBruteForce(wire, step, s):std::numeric_limits<double>::quiet_NaN();
}

//------------------
// DistToRT
//------------------
double DReferenceTrajectory::DistToRT(const DCoordinateSystem *wire, const swim_step_t *step, double *s) const
{
	/// Calculate the distance of the given wire(in the lab
	/// reference frame) to the Reference Trajectory which the
	/// given swim step belongs to. This uses the momentum directions
	/// and positions of the swim step
	/// to define a curve and calculate the distance of the hit
	/// from it. The swim step should be the closest one to the wire.
	/// IMPORTANT: This approximates the helix locally by a parabola.
	/// This means the swim step should be fairly close
	/// to the wire so that this approximation is valid. If the
	/// reference trajectory from which the swim step came is too
	/// sparse, the results will not be nearly as good.
	
	// Interestingly enough, this is one of the harder things to figure
	// out in the tracking code which is why the explanations may be
	// a bit long.

	// The general idea is to define the helix in a coordinate system
	// in which the wire runs along the z-axis. The distance to the
	// wire is then defined just in the X/Y plane of this coord. system.
	// The distance is expressed as a function of the phi angle in the
	// natural coordinate system of the helix. This way, phi=0 corresponds
	// to the swim step point itself and the DOCA point should be
	// at a small phi angle.
	//
	// The minimum distance between the helical segment and the wire
	// will be a function of sin(phi), cos(phi) and phi. Approximating
	// sin(phi) by phi and cos(phi) by (1-phi^2) leaves a 4th order
	// polynomial in phi. Taking the derivative leaves a 3rd order
	// polynomial whose root is the phi corresponding to the 
	// Distance Of Closest Approach(DOCA) point on the helix. Plugging
	// that value of phi back into the distance formula gives
	// us the minimum distance between the track and the wire.

	// First, we need to define the coordinate system in which the 
	// wire runs along the z-axis. This is actually done already
	// in the CDC package for each wire once, at program start.
	// The directions of the axes are defined in wire->sdir,
	// wire->tdir, and wire->udir.
	
	// Next, define a point on the helical segment defined by the
	// swim step it the RT coordinate system. The directions of
	// the RT coordinate system are defined by step->xdir, step->ydir,
	// and step->zdir. The coordinates of a point on the helix
	// in this coordinate system are:
	//
	//   x = Ro*(cos(phi) - 1)
	//   y = Ro*sin(phi)
	//   z = phi*(dz/dphi)
	//
	// where phi is the phi angle of the point in this coordinate system.
	
	// Now, a vector describing the helical point in the LAB coordinate
	// system is:
	//
	//  h = x*xdir + y*ydir + z*zdir + pos
	//
	// where h,xdir,ydir,zdir and pos are all 3-vectors.
	// xdir,ydir,zdir are unit vectors defining the directions
	// of the RT coord. system axes in the lab coord. system.
	// pos is a vector defining the position of the swim step
	// in the lab coord.system 
	
	// Now we just need to find the extent of "h" in the wire's
	// coordinate system (period . means dot product):
	//
	// s = (h-wpos).sdir
	// t = (h-wpos).tdir
	// u = (h-wpos).udir
	//
	// where wpos is the position of the center of the wire in
	// the lab coord. system and is given by wire->wpos.

	// At this point, the values of s,t, and u repesent a point
	// on the helix in the coord. system of the wire with the
	// wire in the "u" direction and positioned at the origin.
	// The distance(squared) from the wire to the point on the helix
	// is given by:
	//
	// d^2 = s^2 + t^2
	//
	// where s and t are both functions of phi.

	// So, we'll define the values of "s" and "t" above as:
	//
	//   s = A*x + B*y + C*z + D
	//   t = E*x + F*y + G*z + H
	//
	// where A,B,C,D,E,F,G, and H are constants defined below
	// and x,y,z are all functions of phi defined above.
	// (period . means dot product)
	//
	// A = sdir.xdir
	// B = sdir.ydir
	// C = sdir.zdir
	// D = sdir.(pos-wpos)
	//
	// E = tdir.xdir
	// F = tdir.ydir
	// G = tdir.zdir
	// H = tdir.(pos-wpos)
	const DVector3 &xdir = step->sdir;
	const DVector3 &ydir = step->tdir;
	const DVector3 &zdir = step->udir;
	const DVector3 &sdir = wire->sdir;
	const DVector3 &tdir = wire->tdir;
	const DVector3 &udir = wire->udir;
	DVector3 pos_diff = step->origin - wire->origin;
	
	double A = sdir.Dot(xdir);
	double B = sdir.Dot(ydir);
	double C = sdir.Dot(zdir);
	double D = sdir.Dot(pos_diff);

	double E = tdir.Dot(xdir);
	double F = tdir.Dot(ydir);
	double G = tdir.Dot(zdir);
	double H = tdir.Dot(pos_diff);

	// OK, here is the dirty part. Using the approximations given above
	// to write the x and y functions in terms of phi^2 and phi (instead
	// of cos and sin) we put them into the equations for s and t above.
	// Then, inserting those into the equation for d^2 above that, we
	// get a very long equation in terms of the constants A,...H and
	// phi up to 4th order. Combining coefficients for similar powers
	// of phi yields an equation of the form:
	//
	// d^2 = Q*phi^4 + R*phi^3 + S*phi^2 + T*phi + U
	//
	// The dirty part is that it takes the better part of a sheet of
	// paper to work out the relations for Q,...U in terms of
	// A,...H, and Ro, dz/dphi. You can work it out yourself on
	// paper to verify that the equations below are correct.
	double Ro = step->Ro;
	double Ro2 = Ro*Ro;
	double delta_z = step->mom.Dot(step->udir);
	double delta_phi = step->mom.Dot(step->tdir)/Ro;
	double dz_dphi = delta_z/delta_phi;
	double dz_dphi2=dz_dphi*dz_dphi;
	double Ro_dz_dphi=Ro*dz_dphi;

	//	double Q = pow(A*Ro/2.0, 2.0) + pow(E*Ro/2.0, 2.0);
	double Q=0.25*Ro2*(A*A+E*E);
	//	double R = -(2.0*A*B*Ro2 + 2.0*A*C*Ro_dz_dphi + 2.0*E*F*Ro2 + 2.0*E*G*Ro_dz_dphi)/2.0;
	double R = -((A*B+E*F)*Ro2 + (A*C+E*G)*Ro_dz_dphi);
	//	double S = pow(B*Ro, 2.0) + pow(C*dz_dphi,2.0) + 2.0*B*C*Ro_dz_dphi - 2.0*A*D*Ro/2.0
	//+ pow(F*Ro, 2.0) + pow(G*dz_dphi,2.0) + 2.0*F*G*Ro_dz_dphi - 2.0*E*H*Ro/2.0;
	double S= (B*B+F*F)*Ro2+(C*C+G*G)*dz_dphi2+2.0*(B*C+F*G)*Ro_dz_dphi
	  -(A*D+E*H)*Ro;
	//	double T = 2.0*B*D*Ro + 2.0*C*D*dz_dphi + 2.0*F*H*Ro + 2.0*G*H*dz_dphi;	
	double T = 2.0*((B*D+F*H)*Ro + (C*D+G*H)*dz_dphi);
	double U = D*D + H*H;
	
	// Aaarghh! my fingers hurt just from typing all of that!
	//
	// OK, now we differentiate the above equation for d^2 to get:
	//
	// d(d^2)/dphi = 4*Q*phi^3 + 3*R*phi^2 + 2*S*phi + T
	//
	// NOTE: don't confuse "R" with "Ro" in the above equations!
	//
	// Now we have to solve the 3rd order polynomial for the phi value of
	// the point of closest approach on the RT. This is a well documented
	// procedure. Essentially, when you have an equation of the form:
	//
	//  x^3 + a2*x^2 + a1*x + a0 = 0;
	//
	// a change of variables is made such that w = x + a2/3 which leads
	// to a third order poly with no w^2 term:
	//
	//  w^3 + 3.0*b*w + 2*c = 0 
	//
	// where:
	//    b = a1/3 - (a2^2)/9
	//    c = a0/2 - a1*a2/6  + (a2^3)/27
	//
	// The one real root of this is:
	//
	//  w0 = q - p
	//
	// where:
	//    q^3 = d - c
	//    p^3 = d + c
	//    d^2 = b^3 + c^2      (don't confuse with d^2 above!)
	//
	// For us this means that:
	//    a2 = 3*R/(4*Q)
	//    a1 = 2*S/(4*Q)
	//    a0 =   T/(4*Q)
	//
	// A potential problem could occur if Q is at or very close to zero.
	// This situation occurs when both A and E are zero. This would mean
	// that both sdir and tdir are perpendicular to xdir which means
	// xdir is in the same direction as udir (got that?). Physically,
	// this corresponds to the situation when both the momentum and
	// the magnetic field are perpendicular to the wire (though not
	// necessarily perpendicular to each other). This situation can't
	// really occur in the CDC detector where the chambers are well
	// contained in a region where the field is essentially along z as
	// are the wires.
	//
	// Just to be safe, we check that Q is greater than
	// some minimum before solving for phi. If it is too small, we fall
	// back to solving the quadratic equation for phi.
	double phi =0.0;
	if(fabs(Q)>1.0E-6){
	  /*
		double fourQ = 4.0*Q;
		double a2 = 3.0*R/fourQ;
		double a1 = 2.0*S/fourQ;
		double a0 =     T/fourQ;
	  */
	  double one_over_fourQ=0.25/Q;
	  double a2=3.0*R*one_over_fourQ;
	  double a1=2.0*S*one_over_fourQ;
	  double a0=T*one_over_fourQ;
	  double a2sq=a2*a2;
	  /*
		double b = a1/3.0 - a2*a2/9.0;
		double c = a0/2.0 - a1*a2/6.0 + a2*a2*a2/27.0;
	  */
	  double b=ONE_THIRD*(a1-ONE_THIRD*a2sq);
	  double c=0.5*(a0-ONE_THIRD*a1*a2)+a2*a2sq/27.0;
	  double my_d2=b*b*b+c*c;
	  if (my_d2>0){
	    //double d = sqrt(pow(b, 3.0) + pow(c, 2.0)); // occasionally, this is zero. See below
	    double d=sqrt(my_d2);
	    //double q = pow(d - c, ONE_THIRD);
	    //double p = pow(d + c, ONE_THIRD);
	    double q=cbrt(d-c);
	    double p=cbrt(d+c);
	    
	    double w0 = q - p;
	    //phi = w0 - a2/3.0;
	    phi = w0 - ONE_THIRD*a2;
	  }
	  else{
	    // Use DeMoivre's theorem to find the cube root of a complex
	    // number.  In this case there are three real solutions.
	    double d=sqrt(-my_d2);
	    c*=-1.;
	    double temp=sqrt(cbrt(c*c+d*d));
	    double theta1=ONE_THIRD*atan2(d,c);
	    double sum_over_2=temp*cos(theta1);
	    double diff_over_2=-temp*sin(theta1);

	    double phi0=-a2/3+2.*sum_over_2;
	    double phi1=-a2/3-sum_over_2+sqrt(3)*diff_over_2;
	    double phi2=-a2/3-sum_over_2-sqrt(3)*diff_over_2;

	    double d2_0 = U + phi0*(T + phi0*(S + phi0*(R + phi0*Q)));
	    double d2_1 = U + phi1*(T + phi1*(S + phi1*(R + phi1*Q)));
	    double d2_2 = U + phi2*(T + phi2*(S + phi2*(R + phi2*Q)));

	    if (d2_0<d2_1 && d2_0<d2_2){
	      phi=phi0;
	    }
	    else if (d2_1<d2_0 && d2_1<d2_2){
	      phi=phi1;
	    }
	    else{
	      phi=phi2;
	    }
	  }
	}
	
	if(fabs(Q)<=1.0E-6 || !isfinite(phi)){
		double a = 3.0*R;
		double b = 2.0*S;
		double c = 1.0*T;
		phi = (-b + sqrt(b*b - 4.0*a*c))/(2.0*a); 
	}
	
	// The accuracy of this method is limited by how close the step is to the
	// actual minimum. If the value of phi is large then the step size is
	// not too close and we should add another couple of steps in the right
	// place in order to get a more accurate value. Note that while this will
	// increase the time it takes this round, presumably the fitter will be
	// calling this often for each wire and having a high density of points
	// near the wires will just make subsequent calls go quicker. This also
	// allows larger initial step sizes with the high density regions getting
	// filled in as needed leading to overall faster tracking.
#if 0
	if(isfinite(phi) && fabs(phi)>2.0E-4){
		if(dist_to_rt_depth>=3){
			_DBG_<<"3 or more recursive calls to DistToRT(). Something is wrong! bailing ..."<<endl;
			//for(int k=0; k<Nswim_steps; k++){
			//	DVector3 &v = swim_steps[k].origin;
			//	_DBG_<<"  "<<k<<": "<<v.X()<<", "<<v.Y()<<", "<<v.Z()<<endl;
			//}
			//exit(-1);
			return std::numeric_limits<double>::quiet_NaN();
		}
		double scale_step = 1.0;
		double s_range = 1.0*scale_step;
		double step_size = 0.02*scale_step;
		int err = InsertSteps(step, phi>0.0 ? +s_range:-s_range, step_size);			// Add new steps near this step by swimming in the direction of phi
		if(!err){
			step=FindClosestSwimStep(wire);									// Find the new closest step
			if(!step)return std::numeric_limits<double>::quiet_NaN();
			dist_to_rt_depth++;
			double doca = DistToRT(wire, step, s);								// re-call ourself with the new step
			dist_to_rt_depth--;
			return doca;
		}else{
			if(err<0)return std::numeric_limits<double>::quiet_NaN();
			
			// If InsertSteps() returns an error > 0 then it indicates that it
			// was unable to add additional steps (perhaps because there
			// aren't enough spaces available). In that case, we just go ahead
			// and use the phi we have and make the best estimate possible.
		}
	}
#endif
	
	// It is possible at this point that the value of phi corresponds to
	// a point past the end of the wire. We should check for this here and
	// recalculate, if necessary, the DOCA at the end of the wire. First,
	// calculate h (the vector defined way up above) and dot it into the
	// wire's u-direction to get the position of the DOCA point along the
	// wire.
	double x = -0.5*Ro*phi*phi;
	double y = Ro*phi;
	double z = dz_dphi*phi;
	DVector3 h = pos_diff + x*xdir + y*ydir + z*zdir;
	double u = h.Dot(udir);
	if(fabs(u) > wire->L/2.0){
		// Looks like our DOCA point is past the end of the wire.
		// Find phi corresponding to the end of the wire.
		double L_over_2 = u>0.0 ? wire->L/2.0:-wire->L/2.0;
		double a = -0.5*Ro*udir.Dot(xdir);
		double b = Ro*udir.Dot(ydir) + dz_dphi*udir.Dot(zdir);
		double c = udir.Dot(pos_diff) - L_over_2;
		double twoa=2.0*a;
		double sqroot=sqrt(b*b-4.0*a*c);
		double phi1 = (-b + sqroot)/(twoa); 
		double phi2 = (-b - sqroot)/(twoa);
		phi = fabs(phi1)<fabs(phi2) ? phi1:phi2;
		u=L_over_2;
	}
	this->last_dist_along_wire = u;

	// Use phi to calculate DOCA
	double d2 = U + phi*(T + phi*(S + phi*(R + phi*Q)));
	double d = sqrt(d2);

	// Calculate distance along track ("s")
	double dz = dz_dphi*phi;
	double Rodphi = Ro*phi;
	double ds = sqrt(dz*dz + Rodphi*Rodphi);
	if(s)*s=step->s + (phi>0.0 ? ds:-ds);
	if(debug_level>3){
	  _DBG_<<"distance to rt: "<< step->s + (phi>0.0 ? ds:-ds) <<" from step at "<<step->s<<" with ds="<<ds<<" d="<<d<<" dz="<<dz<<" Rodphi="<<Rodphi<<endl;
		_DBG_<<"phi="<<phi<<" U="<<U<<" u="<<u<<endl;
	}

	// Remember phi and step so additional info on the point can be obtained
	this->last_phi = phi;
	this->last_swim_step = step;
	this->last_dz_dphi = dz_dphi;

	return d; // WARNING: This could return nan!
}

//------------------
// DistToRTBruteForce
//------------------
double DReferenceTrajectory::DistToRTBruteForce(const DCoordinateSystem *wire, const swim_step_t *step, double *s) const
{
	/// Calculate the distance of the given wire(in the lab
	/// reference frame) to the Reference Trajectory which the
	/// given swim step belongs to. This uses the momentum directions
	/// and positions of the swim step
	/// to define a curve and calculate the distance of the hit
	/// from it. The swim step should be the closest one to the wire.
	/// IMPORTANT: This calculates the distance using a "brute force"
	/// method of taking tiny swim steps to find the minimum distance.
	/// It is vey SLOW and you should be using DistToRT(...) instead.
	/// This is only here to provide an independent check of DistToRT(...).
	
	const DVector3 &xdir = step->sdir;
	const DVector3 &ydir = step->tdir;
	const DVector3 &zdir = step->udir;
	const DVector3 &sdir = wire->sdir;
	const DVector3 &tdir = wire->tdir;
	DVector3 pos_diff = step->origin - wire->origin;
	
	double Ro = step->Ro;
	double delta_z = step->mom.Dot(step->udir);
	double delta_phi = step->mom.Dot(step->tdir)/Ro;
	double dz_dphi = delta_z/delta_phi;

	// Brute force
	double min_d2 = 1.0E6;
	double phi=M_PI;
	for(int i=-2000; i<2000; i++){
		double myphi=(double)i*0.000005;
		DVector3 d = Ro*(cos(myphi)-1.0)*xdir
	            	+ Ro*sin(myphi)*ydir
						+ dz_dphi*myphi*zdir
						+ pos_diff;

		double d2 = pow(d.Dot(sdir),2.0) + pow(d.Dot(tdir),2.0);
		if(d2<min_d2){
			min_d2 = d2;
			phi = myphi;
			this->last_phi = myphi;
		}
	}
	double d2 = min_d2;
	double d = sqrt(d2);
	this->last_phi = phi;
	this->last_swim_step = step;
	this->last_dz_dphi = dz_dphi;
	
	// Calculate distance along track ("s")
	double dz = dz_dphi*phi;
	double Rodphi = Ro*phi;
	double ds = sqrt(dz*dz + Rodphi*Rodphi);
	if(s)*s=step->s + (phi>0.0 ? ds:-ds);

	return d;
}

//------------------
// Straw_dx
//------------------
double DReferenceTrajectory::Straw_dx(const DCoordinateSystem *wire, double radius) const
{
	/// Find the distance traveled within the specified radius of the
	/// specified wire. This will give the "dx" component of a dE/dx
	/// measurement for cylindrical geometry as we have with straw tubes.
	///
	/// At this point, the estimate is done using a simple linear 
	/// extrapolation from the DOCA point in the direction of the momentum
	/// to the 2 points at which it itersects the given radius. Segments
	/// which extend past the end of the wire will be clipped to the end
	/// of the wire before calculating the total dx.
	
	// First, find the DOCA point for this wire
	double s;
	double doca = DistToRT(wire, &s);
	if(!isfinite(doca))
		return 0.0;

	// If doca is outside of the given radius, then we're done
	if(doca>=radius)return 0.0;
	
	// Get the location and momentum direction of the DOCA point
	DVector3 pos, momdir;
	GetLastDOCAPoint(pos, momdir);
	if(momdir.Mag()!=0.0)momdir.SetMag(1.0);
	
	// Get wire direction
	const DVector3 &udir = wire->udir;
	
	// Calculate vectors used to form quadratic equation for "alpha"
	// the distance along the mometum direction from the DOCA point
	// to the intersection with a cylinder of the given radius.
	DVector3 A = udir.Cross(pos-wire->origin);
	DVector3 B = udir.Cross(momdir);
	
	// If the magnitude of B is zero at this point, it means the momentum
	// direction is parallel to the wire. In this case, this method will
	// not work. Return NaN.
	if(B.Mag()<1.0E-10)return std::numeric_limits<double>::quiet_NaN();
	
	double a = B.Mag();
	double b = A.Dot(B);
	double c = A.Mag() - radius;
	double d = sqrt(b*b - 4.0*a*c);
	
	// The 2 roots should correspond to the 2 intersection points.
	double alpha1 = (-b + d)/(2.0*a);
	double alpha2 = (-b - d)/(2.0*a);
	
	DVector3 int1 = pos + alpha1*momdir;
	DVector3 int2 = pos + alpha2*momdir;
	
	// Check if point1 is past the end of the wire
	double q = udir.Dot(int1 - wire->origin);
	if(fabs(q) > wire->L/2.0){
		double gamma = udir.Dot(wire->origin - pos) + (q>0.0 ? +1.0:-1.0)*wire->L/2.0;
		gamma /= momdir.Dot(udir);
		int1 = pos + gamma*momdir;
	}

	// Check if point2 is past the end of the wire
	q = udir.Dot(int2 - wire->origin);
	if(fabs(q) > wire->L/2.0){
		double gamma = udir.Dot(wire->origin - pos) + (q>0.0 ? +1.0:-1.0)*wire->L/2.0;
		gamma /= momdir.Dot(udir);
		int2 = pos + gamma*momdir;
	}
	
	// Calculate distance
	DVector3 delta = int1 - int2;
	
	return delta.Mag();
}

//------------------
// GetLastDOCAPoint
//------------------
void DReferenceTrajectory::GetLastDOCAPoint(DVector3 &pos, DVector3 &mom) const
{
	/// Use values saved by the last call to one of the DistToRT functions
	/// to calculate the 3-D DOCA position in lab coordinates and momentum
	/// in GeV/c.
	
	if(last_swim_step==NULL){
		if(Nswim_steps>0){
			last_swim_step = &swim_steps[0];
			last_phi = 0.0;
		}else{
			pos.SetXYZ(NaN,NaN,NaN);
			mom.SetXYZ(NaN,NaN,NaN);
			return;
		}
	}

	// If last_phi is not finite, set it to 0 as a last resort
	if(!isfinite(last_phi))last_phi = 0.0;
	
	const DVector3 &xdir = last_swim_step->sdir;
	const DVector3 &ydir = last_swim_step->tdir;
	const DVector3 &zdir = last_swim_step->udir;

	double x = -(last_swim_step->Ro/2.0)*last_phi*last_phi;
	double y = last_swim_step->Ro*last_phi;
	double z = last_dz_dphi*last_phi;

	pos = last_swim_step->origin + x*xdir + y*ydir + z*zdir;
	mom = last_swim_step->mom;

	mom.Rotate(-last_phi, zdir);
}

//------------------
// GetLastDOCAPoint
//------------------
DVector3 DReferenceTrajectory::GetLastDOCAPoint(void) const
{
	/// Use values saved by the last call to one of the DistToRT functions
	/// to calculate the 3-D DOCA position in lab coordinates. This is
	/// mainly intended for debugging.
	if(last_swim_step==NULL){
		if(Nswim_steps>0){
			last_swim_step = &swim_steps[0];
			last_phi = 0.0;
		}else{
			return DVector3(NaN,NaN,NaN);
		}
	}
	const DVector3 &xdir = last_swim_step->sdir;
	const DVector3 &ydir = last_swim_step->tdir;
	const DVector3 &zdir = last_swim_step->udir;
	double Ro = last_swim_step->Ro;
	double delta_z = last_swim_step->mom.Dot(zdir);
	double delta_phi = last_swim_step->mom.Dot(ydir)/Ro;
	double dz_dphi = delta_z/delta_phi;

	double x = -(Ro/2.0)*last_phi*last_phi;
	double y = Ro*last_phi;
	double z = dz_dphi*last_phi;

	return last_swim_step->origin + x*xdir + y*ydir + z*zdir;
}

//------------------
// dPdx
//------------------
double DReferenceTrajectory::dPdx_from_A_Z_rho(double ptot, double A, double Z, double density) const
{
	double I = (Z*12.0 + 7.0)*1.0E-9; // From Leo 2nd ed. pg 25.
	if (Z>=13) I=(9.76*Z+58.8*pow(Z,-0.19))*1.0e-9;
	double rhoZ_overA=density*Z/A;
	double KrhoZ_overA = 0.1535e-3*rhoZ_overA;

	return dPdx(ptot, KrhoZ_overA,rhoZ_overA,log(I));
}

//------------------
// dPdx
//------------------
double DReferenceTrajectory::dPdx(double ptot, double KrhoZ_overA, 
				  double rhoZ_overA,double LogI) const
{
	/// Calculate the momentum loss per unit distance traversed of the material with
	/// the given A, Z, and density. Value returned is in GeV/c per cm
	/// This follows the July 2008 PDG section 27.2 ppg 268-270.
	if(mass==0.0)return 0.0; // no ionization losses for neutrals
	
	double gammabeta = ptot/mass;
	double gammabeta2=gammabeta*gammabeta;
	double gamma = sqrt(gammabeta2+1.);
	double beta = gammabeta/gamma;
	double beta2=beta*beta;
	double me = 0.511E-3;
	double m_ratio=me/mass;
	double two_me_gammabeta2=2.*me*gammabeta2;

	double Tmax = two_me_gammabeta2/(1.0+2.0*gamma*m_ratio+m_ratio*m_ratio);
	//double K = 0.307075E-3; // GeV gm^-1 cm^2
	// Density effect
	double delta=0.;	
	double X=log10(gammabeta);
	double X0,X1;
	double Cbar=2.*(LogI-log(28.816e-9*sqrt(rhoZ_overA)))+1.;
	if (rhoZ_overA>0.01){ // not a gas
	  if (LogI<-1.6118){ // I<100
	    if (Cbar<=3.681) X0=0.2;
	    else X0=0.326*Cbar-1.;
	    X1=2.;
	  }
	  else{
	    if (Cbar<=5.215) X0=0.2;
	    else X0=0.326*Cbar-1.5;
	    X1=3.;
	  }
	}
	else { // gases
	  X1=4.;
	  if (Cbar<=9.5) X0=1.6;
	  else if (Cbar>9.5 && Cbar<=10.) X0=1.7;
	  else if (Cbar>10 && Cbar<=10.5) X0=1.8;    
	  else if (Cbar>10.5 && Cbar<=11.) X0=1.9;
	  else if (Cbar>11.0 && Cbar<=12.25) X0=2.;
	  else if (Cbar>12.25 && Cbar<=13.804){
	    X0=2.;
	    X1=5.;
	  }
	  else {
	    X0=0.326*Cbar-2.5;
	    X1=5.;
	  } 
	}
	if (X>=X0 && X<X1)
	  delta=4.606*X-Cbar+(Cbar-4.606*X0)*pow((X1-X)/(X1-X0),3.);
	else if (X>=X1)
	  delta= 4.606*X-Cbar;  	

	double dEdx = KrhoZ_overA/beta2*(log(two_me_gammabeta2*Tmax) 
					 -2.*LogI - 2.0*beta2 -delta);

	double dP_dx = dEdx/beta;

	//double g = 0.350/sqrt(-log(0.06));
	//dP_dx *= 1.0 + exp(-pow(ptot/g,2.0)); // empirical for really low momentum particles

	
	if(ploss_direction==kBackward)dP_dx = -dP_dx;

	return dP_dx;
}

//------------------
// Dump
//------------------
void DReferenceTrajectory::Dump(double zmin, double zmax)
{	
	swim_step_t *step = swim_steps;
	for(int i=0; i<Nswim_steps; i++, step++){
		vector<pair<string,string> > item;
		double x = step->origin.X();
		double y = step->origin.Y();
		double z = step->origin.Z();
		if(z<zmin || z>zmax)continue;

		double px = step->mom.X();
		double py = step->mom.Y();
		double pz = step->mom.Z();
		
		cout<<i<<": ";
		cout<<"(x,y,z)=("<<x<<","<<y<<","<<z<<") ";
		cout<<"(px,py,pz)=("<<px<<","<<py<<","<<pz<<") ";
		cout<<"(Ro,s,t)=("<<step->Ro<<","<<step->s<<","<<step->t<<") ";
		cout<<endl;
	}
	
}

// Propagate the covariance matrix for {px,py,pz,x,y,z,t} along the step ds
jerror_t DReferenceTrajectory::PropagateCovariance(double ds,double q,
						   double mass_sq,
						   const DVector3 &mom,
						   const DVector3 &pos,
						   const DVector3 &B,
						   TMatrixFSym &C) const{
  DMatrix J(7,7);

  double one_over_p_sq=1./mom.Mag2();
  double one_over_p=sqrt(one_over_p_sq);
  double px=mom.X();
  double py=mom.Y();
  double pz=mom.Z();
  double Bx=B.x(),By=B.y(),Bz=B.z();

  double ds_over_p=ds*one_over_p;
  double factor=0.003*q*ds_over_p;
  double temp=(Bz*py-Bx*pz)*one_over_p_sq;
  J(0,0)=1-factor*px*temp;
  J(0,1)=factor*(Bz-py*temp);
  J(0,2)=-factor*(By+pz*temp);

  temp=(Bx*pz-Bz*px)*one_over_p_sq;
  J(1,0)=-factor*(Bz+px*temp);
  J(1,1)=1-factor*py*temp;
  J(1,2)=factor*(Bx-pz*temp);

  temp=(By*px-Bx*py)*one_over_p_sq;
  J(2,0)=factor*(By-px*temp);
  J(2,1)=-factor*(Bx+py*temp);
  J(2,2)=1-factor*pz*temp;

  J(3,3)=1.;
  double ds_over_p3=one_over_p_sq*ds_over_p;
  J(3,0)=ds_over_p*(1-px*px*one_over_p_sq);
  J(3,1)=-px*py*ds_over_p3;
  J(3,2)=-px*pz*ds_over_p3;

  J(4,4)=1.;
  J(4,0)=J(3,1);
  J(4,1)=ds_over_p*(1-py*py*one_over_p_sq);
  J(4,2)=-py*pz*ds_over_p3;

  J(5,5)=1.;
  J(5,0)=J(3,2);
  J(5,1)=J(4,2);
  J(5,2)=ds_over_p*(1-pz*pz*one_over_p_sq);
  
  J(6,6)=1.;
  
  double fac2=(-ds/SPEED_OF_LIGHT)*mass_sq*one_over_p_sq*one_over_p_sq
    /sqrt(1.+mass_sq*one_over_p_sq);
  J(6,0)=fac2*px;
  J(6,1)=fac2*py;
  J(6,2)=fac2*pz;
 
  C=C.Similarity(J);

  return NOERROR;
}

// Find the position along a reference trajectory closest to a line.
// The error matrix for the line can also be input via a pointer.  The error
// matrix is expected to be 7x7 with the order {Px,Py,Pz,X,Y,Z,T}.
// Outputs the kinematic data object (including the covariance) at this 
// position, and the doca and the variance on the doca.
jerror_t DReferenceTrajectory::FindPOCAtoLine(const DVector3 &origin,
					      const DVector3 &dir,
					      const DMatrixDSym *covline, 
					      DKinematicData *track_kd,
					      DVector3 &commonpos, double &doca, double &var_doca) const{ 
  const swim_step_t *swim_step=this->swim_steps;

  shared_ptr<TMatrixFSym> cov = (track_kd!=NULL) ? dResourcePool_TMatrixFSym->Get_SharedResource() : nullptr;
  if(track_kd!=NULL)
  {
	  cov->ResizeTo(7, 7);
	  *cov = *(track_kd->errorMatrix());
  }
  doca=1000.;
  double tflight=0.;
  double mass_sq=this->mass_sq;
  double q=this->q;
  double step_size=1.0,s=-step_size;
  DVector3 oldpos,oldmom;
  DVector3 point=origin;

  // Find the magnitude of the direction vector
  double pscale=dir.Mag();
  // If the magnitude of the direction vector is zero, don't bother to propagate
  // along a line from the input origin...
  bool move_along_line=(pscale>0)?true:false;

  // Propagate along the reference trajectory, comparing to the line at each
  // step
  for (int i=0;i<this->Nswim_steps-1; i++, swim_step++){
    DVector3 pos=swim_step->origin;     
    DVector3 diff=pos-point;
    double new_doca=diff.Mag();
    if (new_doca>doca){
      if (i==1){  // backtrack to find the true doca
	tflight=0.;
	
	swim_step=this->swim_steps;
   if(track_kd!=NULL)
	   *cov=*track_kd->errorMatrix();
	
	pos=swim_step->origin;
	DVector3 mom=swim_step->mom;
	DMagneticFieldStepper stepper(this->bfield, this->q, &pos, &mom);

	int inew=0;
	while (inew<100){
	  DVector3 B;
	  double ds=stepper.Step(&pos,&B,-0.5);
	  // Compute the revised estimate for the doca
	  diff=pos-point;
	  new_doca=diff.Mag();
	  
	  if(new_doca > doca) break;	

	  // Propagate the covariance matrix of the track along the trajectory
	  if(track_kd!=NULL){
	    this->PropagateCovariance(ds,q,mass_sq,mom,oldpos,B,*cov);
	  }
	  
	  // Store the current positions, doca and adjust flight times
	  oldpos=pos;
	  doca=new_doca;
	  
	  double one_over_p_sq=1./mom.Mag2();
	  tflight+=ds*sqrt(1.+mass_sq*one_over_p_sq)/SPEED_OF_LIGHT;
	  
	  // New momentum
	  stepper.GetMomentum(mom);

	  oldmom=/*(-1.)*/mom;
	  inew++;

	  // New point on line
	  if (move_along_line){
	    point-=(step_size/pscale)*dir;
	    s-=step_size;
	  }
	}
      }
      if(track_kd!=NULL)
      {
        track_kd->setErrorMatrix(cov);
        track_kd->setMomentum(oldmom);
        track_kd->setPosition(oldpos);
        track_kd->setTime(track_kd->time() + tflight);
      }
      
      // Compute the variance on the doca
      diff=oldpos-point;
      double dx=diff.x();
      double dy=diff.y();
      double dz=diff.z();
      
      if(track_kd==NULL)
        break;
      //calculate var_doca
      if (covline==NULL){
	var_doca=(dx*dx*((*cov)(kX,kX))+dy*dy*((*cov)(kY,kY))
		  +dz*dz*((*cov)(kZ,kZ))+2.*dx*dy*((*cov)(kX,kY))
		  +2.*dx*dz*((*cov)(kX,kZ))+2.*dy*dz*((*cov)(kY,kZ)))
	  /(doca*doca);
      }
      else{
	DMatrixDSym cov2(*covline);
	if (move_along_line){
	  double two_s=2.*s;
	  double s_sq=s*s;
	  cov2(kX,kX)+=two_s*cov2(kPx,kX)+s_sq*cov2(kPx,kPx); 
	  cov2(kY,kY)+=two_s*cov2(kPy,kY)+s_sq*cov2(kPy,kPy);
	  cov2(kZ,kZ)+=two_s*cov2(kPz,kZ)+s_sq*cov2(kPz,kPz);
	}
	var_doca=(dx*dx*((*cov)(kX,kX)+cov2(kX,kX))
		  +dy*dy*((*cov)(kY,kY)+cov2(kY,kY))
		  +dz*dz*((*cov)(kZ,kZ)+cov2(kZ,kZ))
		  +2.*dx*dy*((*cov)(kX,kY)+cov2(kX,kY))
		  +2.*dx*dz*((*cov)(kX,kZ)+cov2(kX,kZ))
		  +2.*dy*dz*((*cov)(kY,kZ)+cov2(kY,kZ)))
	  /(doca*doca);
      }
      break;
    }
    // New point on line
    if (move_along_line){
      point+=(step_size/pscale)*dir;
      s+=step_size;
    }

    // Propagate the covariance matrix of the track along the trajectory
    if(track_kd!=NULL)
    	this->PropagateCovariance(this->swim_steps[i+1].s-swim_step->s,q,mass_sq,swim_step->mom,swim_step->origin,swim_step->B,*cov);

    // Store the current position and doca
    oldpos=pos;
    oldmom=swim_step->mom;
    tflight=swim_step->t;
    doca=new_doca;
  }

	// "Vertex" is mid-point of line connecting the positions of closest
	// approach of the two tracks
	commonpos = 0.5*(oldpos + point);

  return NOERROR;
}

// Find the position along a reference trajectory closest to a given point.
// The error matrix for the point can also be input via a pointer. The error
// matrix is expected to be 7x7, with the order {Px,Py,Pz,X,Y,Z,T}.
//  Outputs the kinematic data object (including the covariance) at this 
// position,and the doca and the variance on the doca.
jerror_t DReferenceTrajectory::FindPOCAtoPoint(const DVector3 &point,
					       const DMatrixDSym *covpoint, 
					       DKinematicData *track_kd,
					       double &doca, double &var_doca) const{ 
  if (track_kd==NULL) return RESOURCE_UNAVAILABLE;
  
  DVector3 dir, commonpos;
  return FindPOCAtoLine(point,dir,covpoint,track_kd,commonpos,doca,var_doca);
}

// Find the mid-point of the line connecting the points of closest approach of the
// trajectories of two tracks.  Return the positions, momenta, and error matrices 
// at these points for the two tracks.
jerror_t DReferenceTrajectory::IntersectTracks(const DReferenceTrajectory *rt2, DKinematicData *track1_kd, DKinematicData *track2_kd, DVector3 &pos, double &doca, double &var_doca) const {
  const swim_step_t *swim_step1=this->swim_steps;
  const swim_step_t *swim_step2=rt2->swim_steps;
  
  TMatrixFSym cov1(7), cov2(7);
  shared_ptr<TMatrixFSym> locCovarianceMatrix1 = (track1_kd != NULL) ? dResourcePool_TMatrixFSym->Get_SharedResource() : nullptr;
  shared_ptr<TMatrixFSym> locCovarianceMatrix2 = (track2_kd != NULL) ? dResourcePool_TMatrixFSym->Get_SharedResource() : nullptr;

  if((track1_kd != NULL) && (track2_kd != NULL)){
	  locCovarianceMatrix1->ResizeTo(7, 7);
	  locCovarianceMatrix2->ResizeTo(7, 7);
    cov1=*track1_kd->errorMatrix();
    cov2=*track2_kd->errorMatrix();
  }

  double q1=this->q;
  double q2=rt2->q;
  double mass_sq1=this->mass_sq;
  double mass_sq2=rt2->mass_sq;
  
  // Initialize the doca and traverse both particles' trajectories
  doca=1000.;
  double tflight1=0.,tflight2=0.;
  for (int i=0;i<this->Nswim_steps-1&&i<rt2->Nswim_steps-1; i++, swim_step1++, swim_step2++){
    DVector3 pos1=swim_step1->origin;
    DVector3 pos2=swim_step2->origin;
    DVector3 diff=pos1-pos2;
    double new_doca=diff.Mag();
    
    if (new_doca>doca){
      int prev_i=i-1;		   
      // positions and momenta of tracks at the center of the 
      // bracketed region
      pos1=this->swim_steps[prev_i].origin;
      DVector3 mom1=this->swim_steps[prev_i].mom;
      pos2=rt2->swim_steps[prev_i].origin;
      DVector3 mom2=rt2->swim_steps[prev_i].mom;

      // If we break out of the loop immediately, we have not bracketed the 
      // doca yet...
      if (i==1) {  // backtrack to find the true doca
	tflight1=tflight2=0.;
	if((track1_kd != NULL) && (track2_kd != NULL)){
	  cov1=*track1_kd->errorMatrix();
	  cov2=*track2_kd->errorMatrix();
	}
	// Initialize the steppers
	DMagneticFieldStepper stepper1(this->bfield, q1, &pos1, &mom1);
	DMagneticFieldStepper stepper2(this->bfield, q2, &pos2, &mom2);

	// Do the backtracking...
	int inew=0;
	DVector3 oldpos1=pos1;
	DVector3 oldpos2=pos2;
	while (inew<20){
	  if (pos1.z()<0. || pos2.z()<0. || pos1.z()>400. || pos2.z()>400.
	      || pos1.Perp()>65. || pos2.Perp()>65.){
	    break;
	  }
	  DVector3 B1,B2;
	  double ds1=stepper1.Step(&pos1,&B1,-0.5);
	  double ds2=stepper2.Step(&pos2,&B2,-0.5);
	  
	  // Compute the revised estimate for the doca
	  diff=pos1-pos2;
	  new_doca=diff.Mag();
	  
	  if(new_doca > doca){
	    pos1=oldpos1;
	    pos2=oldpos2;
	    break;
	  }
	  
	  // Propagate the covariance matrices along the trajectories
	  if((track1_kd != NULL) && (track2_kd != NULL)){
	    this->PropagateCovariance(ds1,q1,mass_sq1,mom1,oldpos1,B1,cov1);
	    rt2->PropagateCovariance(ds2,q2,mass_sq2,mom2,oldpos2,B2,cov2);
	  }
				  
	  // Store the current positions, doca and adjust flight times
	  oldpos1=pos1;
	  oldpos2=pos2;
	  doca=new_doca;
	  
	  double one_over_p1_sq=1./mom1.Mag2();
	  tflight1+=ds1*sqrt(1.+mass_sq1*one_over_p1_sq)/SPEED_OF_LIGHT;
				  
	  double one_over_p2_sq=1./mom2.Mag2();
	  tflight2+=ds2*sqrt(1.+mass_sq2*one_over_p2_sq)/SPEED_OF_LIGHT;
	  
	  // New momenta
	  stepper1.GetMomentum(mom1);
	  stepper2.GetMomentum(mom2);
	}
      }
      	  
      // Use Brent's algorithm to find a better approximation for 
      // the poca of the two tracks
      double ds=0.5;
      BrentsAlgorithm(pos1,mom1,pos2,mom2,ds,q2,doca);

      // "Vertex" is mid-point of line connecting the positions of closest
      // approach of the two tracks
      pos=0.5*(pos1+pos2);

      if((track1_kd != NULL) && (track2_kd != NULL)){
	// Adjust flight times
	double one_over_p1_sq=1./mom1.Mag2();
	tflight1+=ds*sqrt(1.+mass_sq1*one_over_p1_sq)/SPEED_OF_LIGHT;
				  
	double one_over_p2_sq=1./mom2.Mag2();
	tflight2+=ds*sqrt(1.+mass_sq2*one_over_p2_sq)/SPEED_OF_LIGHT;

    *locCovarianceMatrix1 = cov1;
    track1_kd->setErrorMatrix(locCovarianceMatrix1);
	track1_kd->setMomentum(mom1);
	track1_kd->setPosition(pos1);
	track1_kd->setTime(track1_kd->time() + tflight1);

    *locCovarianceMatrix2 = cov2;
	track2_kd->setErrorMatrix(locCovarianceMatrix2);
	track2_kd->setMomentum(mom2);
	track2_kd->setPosition(pos2);
	track2_kd->setTime(track2_kd->time() + tflight2);
	
	// Compute the variance on the doca
	diff=pos1-pos2;
	double dx=diff.x();
	double dy=diff.y();
	double dz=diff.z();
	var_doca=(dx*dx*(cov1(kX,kX)+cov2(kX,kX))
		  +dy*dy*(cov1(kY,kY)+cov2(kY,kY))
		  +dz*dz*(cov1(kZ,kZ)+cov2(kZ,kZ))
		  +2.*dx*dy*(cov1(kX,kY)+cov2(kX,kY))
				    +2.*dx*dz*(cov1(kX,kZ)+cov2(kX,kZ))
		  +2.*dy*dz*(cov1(kY,kZ)+cov2(kY,kZ)))
	  /(doca*doca);
      }      
      break;
    }

    // Propagate the covariance matrices along the trajectories
    if((track1_kd != NULL) && (track2_kd != NULL)){
      this->PropagateCovariance(this->swim_steps[i+1].s-swim_step1->s,q1,mass_sq1,swim_step1->mom,swim_step1->origin,swim_step1->B,cov1);
      rt2->PropagateCovariance(rt2->swim_steps[i+1].s-swim_step2->s,q2,mass_sq2,swim_step2->mom,swim_step2->origin,swim_step2->B,cov2);
    }
    
    // Store the current positions and doca
    tflight1=swim_step1->t;
    tflight2=swim_step2->t;
    doca=new_doca;
  }
  
  return NOERROR;
}


// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
#define ITMAX 20
#define CGOLD 0.3819660
#define EPS2  1.e-4
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))
jerror_t DReferenceTrajectory::BrentsAlgorithm(DVector3 &pos1,DVector3 &mom1,
					       DVector3 &pos2,DVector3 &mom2,
					       double ds,double q2,
					       double &doca) const{
  double d=0.,u=0.;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-ds;
  double cx=-2.*ds;

  double a=(ax<cx?ax:cx);
  double b=(ax>cx?ax:cx);
  double x=bx,w=bx,v=bx;

  // initialization
  double fw=doca;
  double fv=fw;
  double fx=fw;
  double u_old=x;
  DMagneticFieldStepper stepper1(this->bfield, this->q, &pos1, &mom1);
  DMagneticFieldStepper stepper2(this->bfield, q2, &pos2, &mom2);

  // main loop
  for (unsigned int iter=1;iter<=ITMAX;iter++){
    double xm=0.5*(a+b);
    double tol1=EPS2*fabs(x)+ZEPS;
    double tol2=2.0*tol1;
    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      doca=(pos1-pos2).Mag();

      // New momenta
      stepper1.GetMomentum(mom1);
      stepper2.GetMomentum(mom2);

      return NOERROR;
    }
    // trial parabolic fit
    if (fabs(e)>tol1){
      double x_minus_w=x-w;
      double x_minus_v=x-v;
      double r=x_minus_w*(fx-fv);
      double q=x_minus_v*(fx-fw);
      double p=x_minus_v*q-x_minus_w*r;
      q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(q);
      double etemp=e;
      e=d;
      if (fabs(p)>=fabs(0.5*q*etemp) || p<=q*(a-x) || p>=q*(b-x))
	// fall back on the Golden Section technique
	d=CGOLD*(e=(x>=xm?a-x:b-x));
      else{
	// parabolic step
	d=p/q;
	u=x+d;
	if (u-a<tol2 || b-u <tol2)
	  d=SIGN(tol1,xm-x);
      }						
    } else{
      d=CGOLD*(e=(x>=xm?a-x:b-x));
    }
    u=(fabs(d)>=tol1 ? x+d: x+SIGN(tol1,d));

    // Function evaluation
    double du=u_old-u;
    stepper1.Step(&pos1,NULL,du);
    stepper2.Step(&pos2,NULL,du);
    DVector3 diff=pos1-pos2;
    double fu=diff.Mag();
    u_old=u;

    if (fu<=fx){
      if (u>=x) a=x; else b=x;
      SHFT(v,w,x,u);
      SHFT(fv,fw,fx,fu);      
    }
    else {
      if (u<x) a=u; else b=u;
      if (fu<=fw || w==x){
	v=w;
	w=u;
	fv=fw;
	fw=fu;
      }
      else if (fu<=fv || v==x || v==w){
	v=u;
	fv=fu;
      }
    }
  }
  
  // We only get here if there is a convergence issue...
  doca=(pos1-pos2).Mag();
  stepper1.GetMomentum(mom1);
  stepper2.GetMomentum(mom2);

  return NOERROR;
}


