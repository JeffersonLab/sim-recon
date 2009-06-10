//************************************************************************
// DKalmanFilter.cc
//************************************************************************

#include "DKalmanFilter.h"
#include "CDC/DCDCTrackHit.h"
#include "HDGEOMETRY/DLorentzDeflections.h"
#include "HDGEOMETRY/DMaterialMap.h"
#include "HDGEOMETRY/DRootGeom.h"
#include "DANA/DApplication.h"
#include <TDecompLU.h>

#include <TH2F.h>
#include <TROOT.h>

#include <iomanip>
#include <math.h>

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define EPS 3.0e-8
#define EPS2 1.e-4
#define BEAM_RADIUS  0.1 
#define MAX_ITER 25
#define STEP_SIZE 0.25  // 0.25
#define CDC_STEP_SIZE 0.25 // 0.25
#define CDC_FORWARD_STEP_SIZE 0.25
#define NUM_ITER 10
#define Z_MIN 0.
#define Z_MAX 175.
#define R_MAX 64.5
#define SPEED_OF_LIGHT 29.98
#define CDC_DRIFT_SPEED 55e-4
#define VAR_S 0.09

#define ONE_THIRD 0.33333333333333333
#define ONE_SIXTH 0.16666666666666667

// Local boolean routines for sorting
bool DKalmanHit_cmp(DKalmanHit_t *a, DKalmanHit_t *b){
  return a->z<b->z;
}
bool DKalmanFDCHit_cmp(DKalmanFDCHit_t *a, DKalmanFDCHit_t *b){
  return a->z<b->z;
}
bool DKalmanCDCHit_cmp(DKalmanCDCHit_t *a, DKalmanCDCHit_t *b){
  if (a==NULL || b==NULL){
    cout << "Null pointer in CDC hit list??" << endl;
    return false;
  }
  return (b->ring>a->ring);
}
int grkuta_(double *CHARGE, double *STEP, double *VECT, double *VOUT,const DMagneticFieldMap *bfield);

// Variance for position along wire using PHENIX angle dependence, transverse
// diffusion, and an intrinsic resolution of 127 microns.
#define DIFFUSION_COEFF     1.1e-6 // cm^2/s --> 200 microns at 1 cm
#define DRIFT_SPEED           .0055

inline double fdc_y_variance(double alpha,double x){
  double diffusion=2.*DIFFUSION_COEFF*fabs(x)/DRIFT_SPEED;
  return diffusion+0.00040+0.0064*tan(alpha)*tan(alpha);
}


DKalmanFilter::DKalmanFilter(const DMagneticFieldMap *bfield,
			     const DGeometry *dgeom,
			     const DLorentzDeflections *lorentz_def,
			     //const DMaterialMap *material,
			     const DRootGeom *RootGeom){
  this->bfield=bfield;
  this->geom=dgeom;
  this->lorentz_def=lorentz_def;
  //this->material=material;
  this->RootGeom=RootGeom;
  hits.clear();

  // Get the position of the CDC downstream endplate from DGeometry
  geom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);
  endplate_z-=endplate_dz;

  // Beginning of the cdc
  geom->Get("//posXYZ[@volume='CentralDC']/@X_Y_Z",cdc_origin);
  
  // Beginning of the FDC
  geom->Get("//posXYZ[@volume='ForwardDC']/@X_Y_Z",fdc_origin);
  vector<double>fdc_z1;
  geom->Get("//composition[@name='ForwardDC']/posXYZ[@volume='forwardDC']/@X_Y_Z", fdc_z1);
  fdc_origin[2]+=fdc_z1[2];
  geom->Get("//posXYZ[@volume='forwardDC_package_1']/@X_Y_Z",fdc_z1);
  fdc_origin[2]+=fdc_z1[2]; 
  geom->Get("//posXYZ[@volume='forwardDC_chamber_1']/@X_Y_Z/layer[@value='1']", fdc_z1);
  fdc_origin[2]+=fdc_z1[2]-1.; 

  // Number degrees of freedom
  ndf=0;
  
  // Mass hypothesis
  MASS=0.13957; //charged pion

  DEBUG_HISTS=true;
  // DEBUG_HISTS=false;
  DEBUG_LEVEL=0;
  //DEBUG_LEVEL=2;
}

// Initialize the state vector
jerror_t DKalmanFilter::SetSeed(double q,DVector3 pos, DVector3 mom){
  // Forward parameterization 
  x_=pos(0);
  y_=pos(1);
  z_=pos(2);
  tx_= mom(0)/mom(2);
  ty_= mom(1)/mom(2);
  q_over_p_=q/mom.Mag();
  
  // Central parameterization
  double pt=mom.Perp();
  phi_=mom.Phi();
  tanl_=tan(M_PI/2.-mom.Theta());
  q_over_pt_=q/pt;
  
  return NOERROR;
}

// Return the momentum at the distance of closest approach to the origin.
void DKalmanFilter::GetMomentum(DVector3 &mom){
  double pt=1./fabs(q_over_pt_);
  mom.SetXYZ(pt*cos(phi_),pt*sin(phi_),pt*tanl_);
}

// Return the "vertex" position (position at which track crosses beam line)
void DKalmanFilter::GetPosition(DVector3 &pos){
  pos.SetXYZ(x_,y_,z_);
}

// Add FDC hits
jerror_t DKalmanFilter::AddFDCHit(const DFDCPseudo *fdchit){
  DKalmanFDCHit_t *hit= new DKalmanFDCHit_t;
  
  hit->t=fdchit->time;
  hit->uwire=fdchit->w;
  hit->vstrip=fdchit->s;
  hit->z=fdchit->wire->origin.z();
  hit->cosa=fdchit->wire->udir(1);
  hit->sina=fdchit->wire->udir(0);
  hit->nr=0.;
  hit->nz=0.;
  hit->covu=hit->covv=0.0004;

  fdchits.push_back(hit);
  
  return NOERROR;
}

/// Add a hit to the list of hits using Cartesian coordinates
jerror_t DKalmanFilter::AddHit(double x,double y, double z,double covx,
                                double covy, double covxy, double dE){

  DKalmanHit_t *hit = new DKalmanHit_t;
  hit->x = x;
  hit->y = y;
  hit->z = z;
  hit->covx=covx;
  hit->covy=covy;
  hit->covxy=covxy;
  hit->dE=dE;
  
  hits.push_back(hit);
  
  return NOERROR;
}

//  Add CDC hits
jerror_t DKalmanFilter::AddCDCHit (const DCDCTrackHit *cdchit){
  DKalmanCDCHit_t *hit= new DKalmanCDCHit_t;
  
  hit->origin=cdchit->wire->origin;
  hit->t=cdchit->tdrift;
  hit->d=cdchit->dist;
  hit->stereo=cdchit->wire->stereo;
  //hit->dir=(1./cdchit->wire->udir.z())*cdchit->wire->udir;
  hit->dir=cdchit->wire->udir;
  hit->ring=cdchit->wire->ring;
  
  cdchits.push_back(hit);
  
  return NOERROR;
}

jerror_t DKalmanFilter::AddVertex(DVector3 vertex){
  DKalmanCDCHit_t *hit= new DKalmanCDCHit_t;
  
  hit->origin=vertex;
  hit->t=0.;
  hit->d=0.;
  hit->stereo=0.;
  hit->dir.SetXYZ(0.,0.,1.);
  
  cdchits.push_back(hit);
  
  return NOERROR;
}



// Calculate the derivative of the state vector with respect to z
jerror_t DKalmanFilter::CalcDeriv(double z,double dz,DMatrix S, double dEdx, 
				  DMatrix &D){
  double x=S(state_x,0), y=S(state_y,0),tx=S(state_tx,0),ty=S(state_ty,0);
  double q_over_p=S(state_q_over_p,0);
  double factor=sqrt(1.+tx*tx+ty*ty);

  //B-field at (x,y,z)
  double Bx,By,Bz;
  bfield->GetField(x,y,z, Bx, By, Bz);
  
  D(state_x,0)=tx
    +0.5*dz*qBr2p*q_over_p*factor*(ty*Bz+tx*ty*Bx-(1.+tx*tx)*By);
  D(state_y,0)=ty
    +0.5*dz*qBr2p*q_over_p*factor*(Bx*(1.+ty*ty)-tx*ty*By-tx*Bz);
  D(state_tx,0)=qBr2p*q_over_p*factor*(tx*ty*Bx-(1.+tx*tx)*By+ty*Bz);
  D(state_ty,0)=qBr2p*q_over_p*factor*((1.+ty*ty)*Bx-tx*ty*By-tx*Bz);
  
  if (fabs(dEdx)>0.){
    double E=sqrt(1./q_over_p/q_over_p+MASS*MASS); 
    D(state_q_over_p,0)=-q_over_p*q_over_p*q_over_p*E*dEdx*factor;  
  }
  return NOERROR;
}


// Calculate the derivative of the state vector with respect to z and the 
// Jacobian matrix relating the state vector at z to the state vector at z+dz.
jerror_t DKalmanFilter::CalcDerivAndJacobian(double z,double dz,DMatrix S,
					     double dEdx,
					     DMatrix &J,DMatrix &D){
  double x=S(state_x,0), y=S(state_y,0),tx=S(state_tx,0),ty=S(state_ty,0);
  double q_over_p=S(state_q_over_p,0);
  double factor=sqrt(1.+tx*tx+ty*ty);

  //B-field and field gradient at (x,y,z)
  double Bx,By,Bz,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bfield->GetField(x,y,z, Bx, By, Bz);
  bfield->GetFieldGradient(x,y,z,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,
			   dBzdy,dBzdz);

  // Derivative of S with respect to z
  D(state_x,0)=tx
    +0.5*dz*qBr2p*q_over_p*factor*(ty*Bz+tx*ty*Bx-(1.+tx*tx)*By);
  D(state_y,0)=ty
    +0.5*dz*qBr2p*q_over_p*factor*(Bx*(1.+ty*ty)-tx*ty*By-tx*Bz);
  D(state_tx,0)=qBr2p*q_over_p*factor*(tx*ty*Bx-(1.+tx*tx)*By+ty*Bz);
  D(state_ty,0)=qBr2p*q_over_p*factor*((1.+ty*ty)*Bx-tx*ty*By-tx*Bz);

  // Jacobian
  J(state_x,state_tx)=J(state_y,state_ty)=1.;
  J(state_tx,state_q_over_p)=qBr2p*factor*(Bx*tx*ty-By*(1.+tx*tx)+Bz*ty);
  J(state_ty,state_q_over_p)=qBr2p*factor*(Bx*(1.+ty*ty)-By*tx*ty-Bz*tx);
  J(state_tx,state_tx)=qBr2p*q_over_p*(Bx*ty*(1.+2.*tx*tx+ty*ty)
				       -By*tx*(3.+3.*tx*tx+2.*ty*ty)
				       +Bz*tx*ty)/factor;
  J(state_tx,state_x)=qBr2p*q_over_p*factor*(ty*dBzdx+tx*ty*dBxdx
					     -(1.+tx*tx)*dBydx);
  J(state_ty,state_ty)=qBr2p*q_over_p*(Bx*ty*(3.+2.*tx*tx+3.*ty*ty)
				       -By*tx*(1.+tx*tx+2.*ty*ty)
				       -Bz*tx*ty)/factor;
  J(state_ty,state_y)= qBr2p*q_over_p*factor*((1.+ty*ty)*dBxdy
					      -tx*ty*dBydy-tx*dBzdy);
  J(state_tx,state_ty)=qBr2p*q_over_p*((Bx*tx+Bz)*(1.+tx*tx+2.*ty*ty)
				       -By*ty*(1.+tx*tx))/factor;
  J(state_tx,state_y)= qBr2p*q_over_p*factor*(tx*dBzdy+tx*ty*dBxdy
					      -(1.+tx*tx)*dBydy);
  J(state_ty,state_tx)=-qBr2p*q_over_p*((By*ty+Bz)*(1.+2.*tx*tx+ty*ty)
					-Bx*tx*(1.+ty*ty))/factor;
  J(state_ty,state_x)=qBr2p*q_over_p*factor*((1.+ty*ty)*dBxdx-tx*ty*dBydx
					     -tx*dBzdx);
  J(state_q_over_p,state_tx)=D(state_q_over_p,0)*tx/factor/factor;
  J(state_q_over_p,state_ty)=D(state_q_over_p,0)*ty/factor/factor;

  if (fabs(dEdx)>0.){
    double E=sqrt(1./q_over_p/q_over_p+MASS*MASS); 
    D(state_q_over_p,0)=-q_over_p*q_over_p*q_over_p*E*dEdx*factor;
    J(state_q_over_p,state_q_over_p)=-dEdx*factor/E
      *(2.+3.*MASS*MASS*q_over_p*q_over_p);
  }
   
  return NOERROR;
}

// Reference trajectory for forward tracks in CDC region
// At each point we store the state vector and the Jacobian needed to get to this state 
// along z from the previous state.
jerror_t DKalmanFilter::SetCDCForwardReferenceTrajectory(DMatrix &S,DMatrix &C){
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  DMatrix S0(5,1);

  // Initialize some variables
  temp.h_id=0;
  temp.num_hits=0;
  int i=0,my_i=0,forward_traj_cdc_length=forward_traj_cdc.size();
  double dEdx=0.;
  double z=z_,newz,ds=0.;
  //double my_endz=endplate_z+STEP_SIZE;
  double my_endz=fdc_origin[2];
  double beta2=1.,q_over_p=1.,varE=0.;

  // doca variables
  double doca,old_doca;

  // Loop over CDC hits
  for (unsigned int m=0;m<cdchits.size();m++){
    // wire origin, direction
    DVector3 origin=cdchits[m]->origin;
    DVector3 dir=cdchits[m]->dir;

    // Position along wire
    DVector3 wirepos=origin+((z-origin.z())/dir.z())*dir;
    
    // doca
    old_doca=sqrt((S(state_x,0)-wirepos.x())*(S(state_x,0)-wirepos.x())
		  +(S(state_y,0)-wirepos.y())*(S(state_y,0)-wirepos.y()));
 
    while(z<my_endz){
      newz=z+STEP_SIZE;

      // State at current position 
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
      temp.s= len;
      i++;
      my_i=forward_traj_cdc_length-i;
      if (i<=forward_traj_cdc_length){
	forward_traj_cdc[my_i].s=temp.s;
	forward_traj_cdc[my_i].h_id=temp.h_id;
	forward_traj_cdc[my_i].pos=temp.pos;
	for (unsigned int j=0;j<5;j++){
	  forward_traj_cdc[my_i].S->operator()(j,0)=S(j,0);
	}
      } 
      else{
	temp.S= new DMatrix(S);
      }
   
      // get material properties from the Root Geometry
      if (RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
	  !=NOERROR){
	_DBG_<< "Material error!"<<endl;
	break;
      }
       
      // Get dEdx for the upcoming step
      dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 

      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,z,temp.X0,S,Q);

      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }	
      
      // Compute the Jacobian matrix for swimming away from target
      StepJacobian(z,newz,S,dEdx,J);

      // Update covariance matrix for swimming away from target
      C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;

      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
      
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);

      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }	
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory data
      if (i<=forward_traj_cdc_length){
	for (unsigned int j=0;j<5;j++){
	  for (unsigned int k=0;k<5;k++){
	    forward_traj_cdc[my_i].Q->operator()(j,k)=Q(j,k);
	    forward_traj_cdc[my_i].J->operator()(j,k)=J(j,k);
	  }
	}
      }
      else{	
	temp.Q= new DMatrix(Q);
	temp.J= new DMatrix(J);
	forward_traj_cdc.push_front(temp);    
      }
   
      z=newz;
      
      // update position along wire
      wirepos=origin+((z-cdchits[m]->origin.z())/dir.z())*dir;
    
      // new doca 
      doca=sqrt((S(state_x,0)-wirepos.x())*(S(state_x,0)-wirepos.x())
		+(S(state_y,0)-wirepos.y())*(S(state_y,0)-wirepos.y()));
      
      // Check if we've passed the true minimum doca...
      if (doca>old_doca){	
	// if the doca is increasing, mark the previous point as being close to the true minumum doca.
	// I am using this tag to indicate where the start the search in the kalman routine.
	if (i<=forward_traj_cdc_length){
	  forward_traj_cdc[my_i].h_id=m+1;
	}
	else{
	  forward_traj_cdc[0].h_id=m+1;
	}

	break;
      }
      old_doca=doca;
    }
  }
  // Continue adding to the trajectory until we have passed the endplate
  while(z<my_endz){
     newz=z+STEP_SIZE;
     
     // State at current position 
     temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
     temp.s= len;  
     i++;
     my_i=forward_traj_cdc_length-i;
     if (i<=forward_traj_cdc_length){
       forward_traj_cdc[my_i].s=temp.s;
       forward_traj_cdc[my_i].h_id=temp.h_id;
       forward_traj_cdc[my_i].pos=temp.pos;
       for (unsigned int j=0;j<5;j++){
	 forward_traj_cdc[my_i].S->operator()(j,0)=S(j,0);
       }
     } 
     else{
       temp.S= new DMatrix(S);
     }
     
     // get material properties from the Root Geometry
     if (RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
	 !=NOERROR){
       _DBG_<<"Material error!"<<endl;
       break;
     }
     
     // Get dEdx for the upcoming step
     dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 

     // Get the contribution to the covariance matrix due to multiple 
     // scattering
     if (do_multiple_scattering)
       GetProcessNoise(MASS,ds,z,temp.X0,S,Q);

     // Energy loss straggling in the approximation of thick absorbers
     if (temp.density>0. && do_energy_loss){
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
     }	
         
     // Compute the Jacobian matrix for swimming away from target
     StepJacobian(z,newz,S,dEdx,J);
     
     // Update covariance matrix for swimming away from target
     C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;

     // Step through field
     ds=Step(z,newz,dEdx,S);
     len+=ds;
     
     // Get the contribution to the covariance matrix due to multiple 
     // scattering
     if (do_multiple_scattering)
       GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);
      
     // Energy loss straggling in the approximation of thick absorbers
     if (temp.density>0. && do_energy_loss){
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
     }	
    
     // Compute the Jacobian matrix
     StepJacobian(newz,z,S,dEdx,J);
     
     // update the trajectory
     if (i<=forward_traj_cdc_length){
       for (unsigned int j=0;j<5;j++){
	 for (unsigned int k=0;k<5;k++){
	   forward_traj_cdc[my_i].Q->operator()(j,k)=Q(j,k);
	   forward_traj_cdc[my_i].J->operator()(j,k)=J(j,k);
	 }
       }
     }
     else{	
       temp.Q= new DMatrix(Q);
       temp.J= new DMatrix(J);
       forward_traj_cdc.push_front(temp);    
     }
 
     z=newz;
   }
  /*
   printf("================== %d %d\n",forward_traj_cdc_length,forward_traj_cdc.size());
   for (unsigned int m=0;m<forward_traj_cdc.size();m++){
     printf("id %d x %f y %f z %f s %f \n",
     forward_traj_cdc[m].h_id,forward_traj_cdc[m].pos.x(),
     forward_traj_cdc[m].pos.y(),forward_traj_cdc[m].pos.z(),
     forward_traj_cdc[m].s);
     }    	   
   */
   
   // Current state vector
   S=*(forward_traj_cdc[0].S);

   // position at the end of the swim
   z_=forward_traj_cdc[0].pos.Z();
   x_=forward_traj_cdc[0].pos.X();
   y_=forward_traj_cdc[0].pos.Y();
   
   return NOERROR;
}
// Same as above without outward propagation of C
jerror_t DKalmanFilter::SetCDCForwardReferenceTrajectory(DMatrix &S){
  int i=0,forward_traj_cdc_length=forward_traj_cdc.size();
  double z=z_;

  // Continue adding to the trajectory until we have reached the endplate
  while(z<endplate_z-CDC_FORWARD_STEP_SIZE){
    double step_size=CDC_FORWARD_STEP_SIZE;
    double r2=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);
    if (r2<81.) step_size=CDC_FORWARD_STEP_SIZE/2.;
    if (PropagateForwardCDC(forward_traj_cdc_length,i,z,step_size,S)!=NOERROR)
      return UNRECOVERABLE_ERROR;   
    z+=step_size;
  }  
  if (PropagateForwardCDC(forward_traj_cdc_length,i,z,endplate_z-z,S)!=NOERROR)
    return UNRECOVERABLE_ERROR;
  z=endplate_z;

  // Propagate through endplate with smaller steps
  while(z<endplate_z+2.*endplate_dz){
    if (PropagateForwardCDC(forward_traj_cdc_length,i,z,0.1,S)!=NOERROR)
      return UNRECOVERABLE_ERROR;
    z+=0.1;
  }

  /*
   printf("================== %d %d\n",forward_traj_cdc_length,forward_traj_cdc.size());
   for (unsigned int m=0;m<forward_traj_cdc.size();m++){
     printf("id %d x %f y %f z %f s %f p %f\n",
     forward_traj_cdc[m].h_id,forward_traj_cdc[m].pos.x(),
     forward_traj_cdc[m].pos.y(),forward_traj_cdc[m].pos.z(),
	    forward_traj_cdc[m].s,fabs(1./forward_traj_cdc[m].S->operator()(state_q_over_p,0)));
     }    	   
  */
   
   // Current state vector
   S=*(forward_traj_cdc[0].S);

   // position at the end of the swim
   z_=forward_traj_cdc[0].pos.Z();
   x_=forward_traj_cdc[0].pos.X();
   y_=forward_traj_cdc[0].pos.Y();
   
   return NOERROR;
}

jerror_t DKalmanFilter::PropagateForwardCDC(int length,int &index,double z,
					    double step,DMatrix &S){
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  int my_i=0;
  double newz=z+step;

  // Initialize some variables
  temp.h_id=0;
  temp.num_hits=0;
  double dEdx=0.;
  double ds=0.;
  double beta2=1.,q_over_p=1.,varE=0.;

  // State at current position 
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
  temp.s=len;  
  
  // get material properties from the Root Geometry
  if (RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
      !=NOERROR){
    _DBG_<<"Material error!"<<endl; 
    return UNRECOVERABLE_ERROR;
  }
  
  index++; 
  if (index<=length){
    my_i=length-index;
    forward_traj_cdc[my_i].s=temp.s;
    forward_traj_cdc[my_i].h_id=temp.h_id;
    forward_traj_cdc[my_i].pos=temp.pos;
    forward_traj_cdc[my_i].X0=temp.X0;
    forward_traj_cdc[my_i].A=temp.A;
    forward_traj_cdc[my_i].Z=temp.Z;
    forward_traj_cdc[my_i].density=temp.density;
    for (unsigned int j=0;j<5;j++){
      forward_traj_cdc[my_i].S->operator()(j,0)=S(j,0);
    }
  } 
  else{
    temp.S= new DMatrix(S);
  }
  
  // Get dEdx for the upcoming step
  dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
  
     // Step through field
  ds=Step(z,newz,dEdx,S);
  len+=ds;
  
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  if (do_multiple_scattering)
    GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);
  
  // Energy loss straggling in the approximation of thick absorbers
  if (temp.density>0. && do_energy_loss){
    q_over_p=S(state_q_over_p,0);
    beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
    varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
    
    Q(state_q_over_p,state_q_over_p)
      =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
  }	
  
  // Compute the Jacobian matrix
  StepJacobian(newz,z,S,dEdx,J);
  
  // update the trajectory
  if (index<=length){
    for (unsigned int j=0;j<5;j++){
      for (unsigned int k=0;k<5;k++){
	forward_traj_cdc[my_i].Q->operator()(j,k)=Q(j,k);
	forward_traj_cdc[my_i].J->operator()(j,k)=J(j,k);
      }
    }
  }
  else{	
    temp.Q= new DMatrix(Q);
    temp.J= new DMatrix(J);
    forward_traj_cdc.push_front(temp);    
  }

  return NOERROR;
}



// Swim the state vector through the field from the start of the reference
// trajectory to the end
jerror_t DKalmanFilter::SwimCentral(DVector3 &pos,DMatrix &Sc){
  double ds=CDC_STEP_SIZE;
  for (int m=central_traj.size()-1;m>0;m--){
    double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
    double dedx=GetdEdx(MASS,q_over_p,central_traj[m].Z,
			central_traj[m].A,central_traj[m].density);
    // Variables for the computation of D at the doca to the wire
    double D=Sc(state_D,0);
    double q=(Sc(state_q_over_pt,0)>0)?1.:-1.;
    double qpt=1./Sc(state_q_over_pt,0);
    double sinphi=sin(Sc(state_phi,0));
    double cosphi=cos(Sc(state_phi,0));
    
    // Magnetic field
    DVector3 B;
    bfield->GetField(pos.x(),pos.y(),pos.z(), B(0), B(1), B(2));
    double Bz_=fabs(B(2));
    DVector3 dpos1=central_traj[m-1].pos-central_traj[m].pos;

    // Propagate the state through the field
    FixedStep(pos,ds,Sc,dedx);

    // update D
    double rc=sqrt(dpos1.Perp2()
		   +2.*(D+qpt/qBr2p/Bz_)*(dpos1.x()*sinphi-dpos1.y()*cosphi)
		   +(D+qpt/qBr2p/Bz_)*(D+qpt/qBr2p/Bz_));
    Sc(state_D,0)=q*rc-qpt/qBr2p/Bz_;
  
  }

  return NOERROR;
}

// Reference trajectory for central tracks
// At each point we store the state vector and the Jacobian needed to get to this state 
// along s from the previous state.
// The tricky part is that we swim out from the target to find Sc and pos along the trajectory 
// but we need the Jacobians for the opposite direction, because the filter proceeds from 
// the outer hits toward the target.
jerror_t DKalmanFilter::SetCDCReferenceTrajectory(DVector3 pos,DMatrix &Sc,
						  DMatrix &Cc){
  DKalmanState_t temp;
  DMatrix J(5,5),JT(5,5);  // State vector Jacobian matrix and its transpose
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix C(5,5);  // covariance matrix 
  
  // Position, step, radius, etc. variables
  DVector3 oldpos; 
  double dedx=0;
  double ds=CDC_STEP_SIZE;
  double beta2=1.,varE=0.; 
  len=0.; 
  
  // Coordinates for outermost cdc hit
  unsigned int id=cdchits.size()-1;
  DVector3 origin=cdchits[id]->origin;
  DVector3 dir=cdchits[id]->dir;

  if (central_traj.size()>0){  // reuse existing deque
    // Reset D to zero
    Sc(state_D,0)=0.;
    for (int m=central_traj.size()-1;m>=0;m--){      
      central_traj[m].s=len;
      central_traj[m].pos=pos;
      for (unsigned int i=0;i<5;i++){
	central_traj[m].S->operator()(i,0)=Sc(i,0);
      }
      
      // update path length
      len+=ds;
      
      // get material properties from the Root Geometry
      if(RootGeom->FindMat(pos,central_traj[m].density,central_traj[m].A,central_traj[m].Z,central_traj[m].X0)
	 !=NOERROR){
	/*
	for (unsigned int m=0;m<central_traj.size();m++){
	  DMatrix S=*(central_traj[m].S);
	  double cosphi=cos(S(state_phi,0));
	  double sinphi=sin(S(state_phi,0));
	  double pt=fabs(1./S(state_q_over_pt,0));
	  double tanl=S(state_tanl,0);
	  
	  cout
	    << setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
	    << central_traj[m].pos.x() << ", "
	    << central_traj[m].pos.y() << ", "
	    << central_traj[m].pos.z() << "  mom: "
	    << pt*cosphi << ", " << pt*sinphi << ", " 
	<< pt*tanl << " -> " << pt/cos(atan(tanl))
	    <<"  s: " << setprecision(3) 	   
	    << central_traj[m].s << endl;
	}
	*/
	_DBG_<< "Material error!" <<endl;
	break;
      }

      // Energy loss	
      double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
      if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz){
	if (do_energy_loss){
	  dedx=GetdEdx(MASS,q_over_p,central_traj[m].Z,
		       central_traj[m].A,central_traj[m].density);
	}
	
	// Multiple scattering    
	if (do_multiple_scattering && Sc(state_z,0)>cdc_origin[2])
	  GetProcessNoiseCentral(MASS,ds,pos,central_traj[m].X0,Sc,Q);
      
	// Energy loss straggling in the approximation of thick absorbers
	if (central_traj[m].density>0. && do_energy_loss){
	  beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	  varE=GetEnergyVariance(ds,q_over_p,central_traj[m].Z,
				 central_traj[m].A,central_traj[m].density);

	  Q(state_q_over_pt,state_q_over_pt)
	    =varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	    *q_over_p*q_over_p;
	}
      }

      // Compute the Jacobian matrix for forward-tracking away from target
      StepJacobian(pos,origin,dir,ds,Sc,dedx,J);
      // Propagate the coviarance matrix
      Cc=J*(Cc*DMatrix(DMatrix::kTransposed,J))+Q;

      // Propagate the state through the field
      FixedStep(pos,ds,Sc,dedx);

      if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz){
	// Energy loss straggling in the approximation of thick absorbers
	if (central_traj[m].density>0. && do_energy_loss){
	  q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	  beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	  varE=GetEnergyVariance(ds,q_over_p,central_traj[m].Z,
				 central_traj[m].A,central_traj[m].density);

	  Q(state_q_over_pt,state_q_over_pt)
	    =varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	    *q_over_p*q_over_p;
	}

	// Multiple scattering    
	if (do_multiple_scattering)
	  GetProcessNoiseCentral(MASS,ds,pos,central_traj[m].X0,Sc,Q);
      }

	// Compute the Jacobian matrix for back-tracking towards target
      StepJacobian(pos,origin,dir,-ds,Sc,dedx,J);
      
      //central_traj[m].s=len;
      //central_traj[m].pos=pos;
      for (unsigned int i=0;i<5;i++){
	//central_traj[m].S->operator()(i,0)=Sc(i,0);
	for (unsigned int j=0;j<5;j++){
	  central_traj[m].J->operator()(i,j)=J(i,j);
	  central_traj[m].Q->operator()(i,j)=Q(i,j);	  
	}
      }
    }
  }
  
  // Swim out
  double r=pos.Perp();
  while(r<R_MAX && pos.z()<Z_MAX && pos.z()>Z_MIN){
    // Reset D to zero
    Sc(state_D,0)=0.;

    // store old position
    oldpos=pos;
    
    temp.pos=pos;	
    temp.s=len;
    temp.S= new DMatrix(Sc);	
 
    // update path length
    len+=ds;
    
    // get material properties from the Root Geometry
    if (RootGeom->FindMat(pos,temp.density,temp.A,temp.Z,temp.X0)!=NOERROR){
      _DBG_<<"Material error!"<<endl;
      break;
    }
    
    // Get dEdx for the upcoming step 
    double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
    if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz){
      if (do_energy_loss){
	dedx=GetdEdx(MASS,q_over_p,temp.Z,temp.A,temp.density);
      }
      
      // Multiple scattering    
      if (do_multiple_scattering)
	GetProcessNoiseCentral(MASS,ds,pos,temp.X0,Sc,Q);
    
      // Energy loss straggling in the approximation of thick absorbers  
      if (temp.density>0. && do_energy_loss){
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);

	Q(state_q_over_pt,state_q_over_pt)
	  =varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	  *q_over_p*q_over_p;
      }
    }

    // Compute the Jacobian matrix for forward-tracking away from target
    StepJacobian(pos,origin,dir,ds,Sc,dedx,J);
    // Propagate the covariance matrix
    Cc=J*(Cc*DMatrix(DMatrix::kTransposed,J))+Q;
    
    // Propagate the state through the field
    FixedStep(pos,ds,Sc,dedx);

    if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz){
      // Multiple scattering    
      if (do_multiple_scattering)
	GetProcessNoiseCentral(MASS,ds,pos,temp.X0,Sc,Q);
      
      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){ 
	q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);

	Q(state_q_over_pt,state_q_over_pt)
	  =varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	  *q_over_p*q_over_p;
      }
    }
    // Compute the Jacobian matrix
    StepJacobian(pos,origin,dir,-ds,Sc,dedx,J);
    
    // update the radius relative to the beam line
    r=pos.Perp();
    
    // Update the trajectory info
    temp.Q= new DMatrix(Q);
    temp.J= new DMatrix(J);
    temp.C= new DMatrix(C);

    central_traj.push_front(temp);    
  }
  
  if (DEBUG_LEVEL>1)
  {
    for (unsigned int m=0;m<central_traj.size();m++){
      DMatrix S=*(central_traj[m].S);
      double cosphi=cos(S(state_phi,0));
      double sinphi=sin(S(state_phi,0));
      double pt=fabs(1./S(state_q_over_pt,0));
      double tanl=S(state_tanl,0);
      
      cout
	<< setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
	<< central_traj[m].pos.x() << ", "
	<< central_traj[m].pos.y() << ", "
	<< central_traj[m].pos.z() << "  mom: "
	<< pt*cosphi << ", " << pt*sinphi << ", " 
	<< pt*tanl << " -> " << pt/cos(atan(tanl))
	<<"  s: " << setprecision(3) 	   
	<< central_traj[m].s << endl;
    }
  }
  // State at end of swim
  Sc=*(central_traj[0].S);

  // Position at the end of the swim
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  return NOERROR;
}
// Same as above without outward projection of covariance matrix
jerror_t DKalmanFilter::SetCDCReferenceTrajectory(DVector3 pos,DMatrix &Sc){
  DKalmanState_t temp;
  DMatrix J(5,5);  // State vector Jacobian matrix 
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix C(5,5);  // covariance matrix 
  
  // Position, step, radius, etc. variables
  DVector3 oldpos; 
  double dedx=0;
  double ds=CDC_STEP_SIZE;
  double beta2=1.,varE=0.,q_over_p=1.; 
  len=0.; 
  
  // Coordinates for outermost cdc hit
  unsigned int id=cdchits.size()-1;
  DVector3 origin=cdchits[id]->origin;
  DVector3 dir=cdchits[id]->dir;

  if (central_traj.size()>0){  // reuse existing deque
    // Reset D to zero
    Sc(state_D,0)=0.;
    for (int m=central_traj.size()-1;m>=0;m--){      
      central_traj[m].s=len;
      central_traj[m].pos=pos;
      for (unsigned int i=0;i<5;i++){
	central_traj[m].S->operator()(i,0)=Sc(i,0);
      }
      // Make sure D is zero
      central_traj[m].S->operator()(state_D,0)=0.;
      
      // update path length
      len+=ds;

      // Initialize energy loss 
      dedx=0.;

      // get material properties from the Root Geometry
      if(RootGeom->FindMat(pos,central_traj[m].density,central_traj[m].A,
			   central_traj[m].Z,central_traj[m].X0)!=NOERROR){
	_DBG_ << "Material error! " << endl;
	break;
      }

      //if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz)
      {
      	// Energy loss	
	if (do_energy_loss){	
	  q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	  dedx=GetdEdx(MASS,q_over_p,temp.Z,temp.A,temp.density);
	}
      }

      // Propagate the state through the field
      FixedStep(pos,ds,Sc,dedx);

      //if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz)
	{
	// Energy loss straggling in the approximation of thick absorbers
	if (central_traj[m].density>0. && do_energy_loss){
	  q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	  beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	  varE=GetEnergyVariance(ds,q_over_p,central_traj[m].Z,
				 central_traj[m].A,central_traj[m].density);

	  Q(state_q_over_pt,state_q_over_pt)
	    =varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	    *q_over_p*q_over_p;
	}

	// Multiple scattering    
	if (do_multiple_scattering)
	  GetProcessNoiseCentral(MASS,ds,pos,central_traj[m].X0,Sc,Q);
      }

	// Compute the Jacobian matrix for back-tracking towards target
      StepJacobian(pos,origin,dir,-ds,Sc,dedx,J);
    
      // Fill the deque with the Jacobian and Process Noise matrices
      for (unsigned int i=0;i<5;i++){
	for (unsigned int j=0;j<5;j++){
	  central_traj[m].J->operator()(i,j)=J(i,j);
	  central_traj[m].Q->operator()(i,j)=Q(i,j);	  
	}
      }
    }
  }
  
  // Swim out
  double r=pos.Perp();
  while(r<R_MAX && pos.z()<Z_MAX && pos.z()>Z_MIN){
    // Reset D to zero
    Sc(state_D,0)=0.;

    // store old position
    oldpos=pos;
    
    temp.pos=pos;	
    temp.s=len;
    temp.S= new DMatrix(Sc);	
 
    // update path length
    len+=ds;
    
    // Initialize energy loss 
    dedx=0.; 

    // get material properties from the Root Geometry
    if(RootGeom->FindMat(pos,temp.density,temp.A,temp.Z,temp.X0)!=NOERROR){
      _DBG_ << "Material error! " << endl;
      break;
    }
       
    //if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz)
    {
      // Energy loss	
      if (do_energy_loss){
	q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	dedx=GetdEdx(MASS,q_over_p,temp.Z,temp.A,temp.density);
      }
    }
     
    // Propagate the state through the field
    FixedStep(pos,ds,Sc,dedx);

    //if (Sc(state_z,0)>cdc_origin[2] && Sc(state_z,0)<endplate_z-endplate_dz)
    {
      // Multiple scattering    
      if (do_multiple_scattering)
	GetProcessNoiseCentral(MASS,ds,pos,temp.X0,Sc,Q);
      
      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){ 
	q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);

	Q(state_q_over_pt,state_q_over_pt)
	  =varE*Sc(state_q_over_pt,0)*Sc(state_q_over_pt,0)/beta2
	  *q_over_p*q_over_p;
      }
    }
    // Compute the Jacobian matrix
    StepJacobian(pos,origin,dir,-ds,Sc,dedx,J);
    
    // update the radius relative to the beam line
    r=pos.Perp();
    
    // Update the trajectory info
    temp.Q= new DMatrix(Q);
    temp.J= new DMatrix(J);
    temp.C= new DMatrix(C);

    central_traj.push_front(temp);    
  }
  
  if (DEBUG_LEVEL>1)
  {
    for (unsigned int m=0;m<central_traj.size();m++){
      DMatrix S=*(central_traj[m].S);
      double cosphi=cos(S(state_phi,0));
      double sinphi=sin(S(state_phi,0));
      double pt=fabs(1./S(state_q_over_pt,0));
      double tanl=S(state_tanl,0);
      
      cout
	<< setiosflags(ios::fixed)<< "pos: " << setprecision(4) 
	<< central_traj[m].pos.x() << ", "
	<< central_traj[m].pos.y() << ", "
	<< central_traj[m].pos.z() << "  mom: "
	<< pt*cosphi << ", " << pt*sinphi << ", " 
	<< pt*tanl << " -> " << pt/cos(atan(tanl))
	<<"  s: " << setprecision(3) 	   
	<< central_traj[m].s << endl;
    }
  }
  // State at end of swim
  Sc=*(central_traj[0].S);

  // Position at the end of the swim
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  return NOERROR;
}

  // Reference trajectory for trajectories with hits in the forward direction
// At each point we store the state vector and the Jacobian needed to get to this state 
// along z from the previous state.
jerror_t DKalmanFilter::SetReferenceTrajectory(DMatrix &S,DMatrix &C){   
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  double ds=0.; // path length increment

  // Initialize some variables
  temp.h_id=0;
  temp.num_hits=0;
  double dEdx=0.;
  double beta2=1.,q_over_p=1.,varE=0.;

   // progress in z from hit to hit
  double z=z_;
  double newz=z;
  int i=0,my_i=0;
  int forward_traj_length=forward_traj.size();
  for (unsigned int m=0;m<fdchits.size();m++){
    int num=int((fdchits[m]->z-z)/STEP_SIZE);
    newz=fdchits[m]->z-STEP_SIZE*double(num);
    
    if (newz-z>0.){
      temp.s=len;
      temp.h_id=0;
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
      i++;
      my_i=forward_traj_length-i;
      if (i<=forward_traj_length){
	forward_traj[my_i].s=temp.s;
	forward_traj[my_i].h_id=temp.h_id;
	forward_traj[my_i].pos=temp.pos;
	for (unsigned int j=0;j<5;j++){
	  forward_traj[my_i].S->operator()(j,0)=S(j,0);
	}
      } 
      else{
	temp.S=new DMatrix(S);
      }

      // get material properties from the Root Geometry
      if (RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
	  !=NOERROR){
	_DBG_ << "Material error! " << endl;
	break;
      }
      
      // Get dEdx for the upcoming step
      if (do_energy_loss){
	dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
      }

      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,z,temp.X0,S,Q);
      
      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){	
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }
      
      // Compute the Jacobian matrix for swimming away from target
      StepJacobian(z,newz,S,dEdx,J);
      
      // Update covariance matrix for swimming away from target
      C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;
      
      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
   
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);
      
      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }			      
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory data
      if (i<=forward_traj_length){
	for (unsigned int j=0;j<5;j++){
	  for (unsigned int k=0;k<5;k++){
	    forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	    forward_traj[my_i].J->operator()(j,k)=J(j,k);
	  }
	}
      }
      else{
	temp.Q=new DMatrix(Q);
	temp.J=new DMatrix(J);
	forward_traj.push_front(temp);
      }
      // update z
      z=newz;
    }
    
    for (int k=0;k<num;k++){
      i++;
      my_i=forward_traj_length-i;
      newz=z+STEP_SIZE;
      
      temp.s=len;
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
      if (i<=forward_traj_length){
	forward_traj[my_i].s=temp.s;
	forward_traj[my_i].h_id=temp.h_id;
	forward_traj[my_i].pos=temp.pos;
	for (unsigned int j=0;j<5;j++){
	  forward_traj[my_i].S->operator()(j,0)=S(j,0);
	}
      } 
      else{
	temp.S=new DMatrix(S);
      }
      
      // get material properties from the Root Geometry
      if(RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
	 !=NOERROR){
	_DBG_<<"Material error!"<<endl; 
	break;
      }
      
      // Get dEdx for the upcoming step
      if (do_energy_loss){
	dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
      }
 
      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,z,temp.X0,S,Q);
      
      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){	
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }
      
      // Compute the Jacobian matrix for swimming away from target
      StepJacobian(z,newz,S,dEdx,J);
      
      // Update covariance matrix for swimming away from target
      C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;
      
      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
      
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);

      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){	
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory data
      if (i<=forward_traj_length){
	for (unsigned int j=0;j<5;j++){
	  for (unsigned int k=0;k<5;k++){
	    forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	    forward_traj[my_i].J->operator()(j,k)=J(j,k);
	  }
	}
      }
      else{
	temp.Q=new DMatrix(Q);
	temp.J=new DMatrix(J);
	forward_traj.push_front(temp);
      }
      temp.h_id=0;
      
      //update z
      z=newz;
    }
    
    // Lorentz correction slope parameters
    double tanr=0.,tanz=0.;
    lorentz_def->GetLorentzCorrectionParameters(S(state_x,0),S(state_y,0),z,
						tanz,tanr);
    fdchits[m]->nr=tanr;
    fdchits[m]->nz=tanz;
    
    temp.h_id=m+1;
  }
  
  // Store final point
  i++;
  my_i=forward_traj_length-i;
  temp.s=len;
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),newz); 
  if (i<=forward_traj_length){
    forward_traj[my_i].s=temp.s;
    forward_traj[my_i].h_id=temp.h_id;
    forward_traj[my_i].pos=temp.pos;
    for (unsigned int j=0;j<5;j++){
      forward_traj[my_i].S->operator()(j,0)=S(j,0);
    }
  }
  else{
    temp.S=new DMatrix(S);
  }

  // get material properties from the Root Geometry
  if(RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)==NOERROR){
  
    // Get dEdx for the upcoming step
    if (do_energy_loss){
      dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
    }

    // multiple scattering
    if (do_multiple_scattering)
      GetProcessNoise(MASS,ds,z,temp.X0,S,Q);

    // Energy loss straggling in the approximation of thick absorbers
    if (temp.density>0. && do_energy_loss){	
      q_over_p=S(state_q_over_p,0);
      beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
      varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
      
      Q(state_q_over_p,state_q_over_p)
	=varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
    }
      
    // Compute the Jacobian matrix for swimming away from target
    StepJacobian(z,newz,S,dEdx,J);
    
    // Update covariance matrix for swimming away from target
    C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;
  
    // One more step after last hit point
    newz=z+STEP_SIZE;
    ds=Step(z,newz,dEdx,S);
    len+=ds;
    
    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    if (do_multiple_scattering)
      GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);
 
    // Energy loss straggling in the approximation of thick absorbers
    if (temp.density>0. && do_energy_loss){	
      q_over_p=S(state_q_over_p,0);
      beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
      varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
      
      Q(state_q_over_p,state_q_over_p)
	=varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
    }
    
    // Compute the Jacobian matrix
    StepJacobian(newz,z,S,dEdx,J);
    
    // update the trajectory data
    if (i<=forward_traj_length){
      for (unsigned int j=0;j<5;j++){
	for (unsigned int k=0;k<5;k++){
	  forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	  forward_traj[my_i].J->operator()(j,k)=J(j,k);
	}
      }
    }
    else{
      temp.Q=new DMatrix(Q);
      temp.J=new DMatrix(J);
      forward_traj.push_front(temp);
    }
  }
  /*
  printf("--------------- %d %d\n",forward_traj_length,forward_traj.size());
    for (unsigned int m=0;m<forward_traj.size();m++){
    printf("id %d x %f y %f z %f s %f \n",
    forward_traj[m].h_id,forward_traj[m].pos.x(),
    forward_traj[m].pos.y(),forward_traj[m].pos.z(),forward_traj[m].s);
    }    
  */

  // position at the end of the swim
  z_=z;
  x_=S(state_x,0);
  y_=S(state_y,0);

  return NOERROR;
}
// Same as above without outward propagation of C
jerror_t DKalmanFilter::SetReferenceTrajectory(DMatrix &S){   
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  double ds=0.; // path length increment

  // Initialize some variables
  temp.h_id=0;
  temp.num_hits=0;
  double dEdx=0.;
  double beta2=1.,q_over_p=1.,varE=0.;

   // progress in z from hit to hit
  double z=z_;
  double newz=z;
  int i=0,my_i=0;
  int forward_traj_length=forward_traj.size();
  for (unsigned int m=0;m<fdchits.size();m++){
    int num=int((fdchits[m]->z-z)/STEP_SIZE);
    newz=fdchits[m]->z-STEP_SIZE*double(num);
    
    if (newz-z>0.){
      temp.s=len;
      temp.h_id=0;
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);
      
      // get material properties from the Root Geometry
      if (RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
	  !=NOERROR){
	_DBG_<<"Material error!"<<endl; 
	break;
      }

      i++;
      my_i=forward_traj_length-i;
      if (i<=forward_traj_length){
	forward_traj[my_i].s=temp.s;
	forward_traj[my_i].h_id=temp.h_id;
	forward_traj[my_i].pos=temp.pos;
	forward_traj[my_i].A=temp.A;
	forward_traj[my_i].Z=temp.Z;
	forward_traj[my_i].density=temp.density;
	forward_traj[my_i].X0=temp.X0;
	for (unsigned int j=0;j<5;j++){
	  forward_traj[my_i].S->operator()(j,0)=S(j,0);
	}
      } 
      else{
	temp.S=new DMatrix(S);
      }

     
      
      // Get dEdx for the upcoming step
      if (do_energy_loss){
	dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
      }

      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
   
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);
      
      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }			      
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory data
      if (i<=forward_traj_length){
	for (unsigned int j=0;j<5;j++){
	  for (unsigned int k=0;k<5;k++){
	    forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	    forward_traj[my_i].J->operator()(j,k)=J(j,k);
	  }
	}
      }
      else{
	temp.Q=new DMatrix(Q);
	temp.J=new DMatrix(J);
	forward_traj.push_front(temp);
      }
      // update z
      z=newz;
    }
    
    for (int k=0;k<num;k++){
      i++;
      my_i=forward_traj_length-i;
      newz=z+STEP_SIZE;
      
      temp.s=len;
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z);

      // get material properties from the Root Geometry
      if(RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)
	 !=NOERROR){
	_DBG_<<"Material error!"<<endl;
	break;
      }

      if (i<=forward_traj_length){
	forward_traj[my_i].s=temp.s;
	forward_traj[my_i].h_id=temp.h_id;
	forward_traj[my_i].pos=temp.pos;
	forward_traj[my_i].A=temp.A;
	forward_traj[my_i].Z=temp.Z;
	forward_traj[my_i].density=temp.density;
	forward_traj[my_i].X0=temp.X0;
	for (unsigned int j=0;j<5;j++){
	  forward_traj[my_i].S->operator()(j,0)=S(j,0);
	}
      } 
      else{
	temp.S=new DMatrix(S);
      }    
      
      // Get dEdx for the upcoming step
      if (do_energy_loss){
	dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
      }
      
      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
      
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      if (do_multiple_scattering)
	GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);

      // Energy loss straggling in the approximation of thick absorbers
      if (temp.density>0. && do_energy_loss){	
	q_over_p=S(state_q_over_p,0);
	beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
	varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
	
	Q(state_q_over_p,state_q_over_p)
	  =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
      }
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory data
      if (i<=forward_traj_length){
	for (unsigned int j=0;j<5;j++){
	  for (unsigned int k=0;k<5;k++){
	    forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	    forward_traj[my_i].J->operator()(j,k)=J(j,k);
	  }
	}
      }
      else{
	temp.Q=new DMatrix(Q);
	temp.J=new DMatrix(J);
	forward_traj.push_front(temp);
      }
      temp.h_id=0;
      
      //update z
      z=newz;
    }
    
    // Lorentz correction slope parameters
    double tanr=0.,tanz=0.;
    lorentz_def->GetLorentzCorrectionParameters(S(state_x,0),S(state_y,0),z,
						tanz,tanr);
    fdchits[m]->nr=tanr;
    fdchits[m]->nz=tanz;
    
    temp.h_id=m+1;
  }

  // Final point 
  i++;
  my_i=forward_traj_length-i;
  temp.s=len;
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),newz); 

  // get material properties from the Root Geometry
  if(RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)==NOERROR){
    
    // Get dEdx for the upcoming step
    if (do_energy_loss){
      dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
    }
    
    // Store next-to-final point
    if (i<=forward_traj_length){
      forward_traj[my_i].s=temp.s;
      forward_traj[my_i].h_id=temp.h_id;
      forward_traj[my_i].pos=temp.pos;
      forward_traj[my_i].A=temp.A;
      forward_traj[my_i].Z=temp.Z;
      forward_traj[my_i].density=temp.density;
      forward_traj[my_i].X0=temp.X0;
      for (unsigned int j=0;j<5;j++){
	forward_traj[my_i].S->operator()(j,0)=S(j,0);
      }
    }
    else{
      temp.S=new DMatrix(S);
    }
  
    // One more step to last hit point
    newz=z+STEP_SIZE;
    ds=Step(z,newz,dEdx,S);
    len+=ds;
    
    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    if (do_multiple_scattering)
      GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);
    
    // Energy loss straggling in the approximation of thick absorbers
    if (temp.density>0. && do_energy_loss){	
      q_over_p=S(state_q_over_p,0);
      beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
      varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
      
      Q(state_q_over_p,state_q_over_p)
	=varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
    }
    
    // Compute the Jacobian matrix
    StepJacobian(newz,z,S,dEdx,J);
    
    // update the trajectory data
    if (i<=forward_traj_length){
      for (unsigned int j=0;j<5;j++){
	for (unsigned int k=0;k<5;k++){
	  forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	  forward_traj[my_i].J->operator()(j,k)=J(j,k);
	}
      }
    }
    else{
      temp.Q=new DMatrix(Q);
      temp.J=new DMatrix(J);
      forward_traj.push_front(temp);
    }
  }

  // Make sure the ref trajectory goes one step beyond the most downstream hit
  i++;
  my_i=forward_traj_length-i;
  temp.s=len;
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),newz); 
  temp.h_id=0;
  
  // get material properties from the Root Geometry
  if(RootGeom->FindMat(temp.pos,temp.density,temp.A,temp.Z,temp.X0)==NOERROR){
    
    // Get dEdx for the upcoming step
    if (do_energy_loss){
      dEdx=GetdEdx(MASS,S(state_q_over_p,0),temp.Z,temp.A,temp.density); 
    }
    
    // Store final point
    if (i<=forward_traj_length){
      forward_traj[my_i].s=temp.s;
      forward_traj[my_i].h_id=temp.h_id;
      forward_traj[my_i].pos=temp.pos;
      forward_traj[my_i].A=temp.A;
      forward_traj[my_i].Z=temp.Z;
      forward_traj[my_i].density=temp.density;
      forward_traj[my_i].X0=temp.X0;
      for (unsigned int j=0;j<5;j++){
	forward_traj[my_i].S->operator()(j,0)=S(j,0);
      }
    }
    else{
      temp.S=new DMatrix(S);
    }
    
    // One more step after last hit point
    newz=z+STEP_SIZE;
    ds=Step(z,newz,dEdx,S);
    len+=ds;
    
    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    if (do_multiple_scattering)
      GetProcessNoise(MASS,ds,newz,temp.X0,S,Q);
  
    // Energy loss straggling in the approximation of thick absorbers
    if (temp.density>0. && do_energy_loss){	
      q_over_p=S(state_q_over_p,0);
      beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
    varE=GetEnergyVariance(ds,q_over_p,temp.Z,temp.A,temp.density);
    
    Q(state_q_over_p,state_q_over_p)
      =varE*q_over_p*q_over_p*q_over_p*q_over_p/beta2;
    }
    
    // Compute the Jacobian matrix
    StepJacobian(newz,z,S,dEdx,J);
    
    // update the trajectory data
    if (i<=forward_traj_length){
      for (unsigned int j=0;j<5;j++){
	for (unsigned int k=0;k<5;k++){
	  forward_traj[my_i].Q->operator()(j,k)=Q(j,k);
	  forward_traj[my_i].J->operator()(j,k)=J(j,k);
	}
    }
    }
    else{
      temp.Q=new DMatrix(Q);
      temp.J=new DMatrix(J);
      forward_traj.push_front(temp);
    }
  }

  if (DEBUG_LEVEL==2){
    printf("--------------- %d %d\n",forward_traj_length,forward_traj.size());
    for (unsigned int m=0;m<forward_traj.size();m++){
      printf("id %d x %f y %f z %f s %f \n",
	     forward_traj[m].h_id,forward_traj[m].pos.x(),
	     forward_traj[m].pos.y(),forward_traj[m].pos.z(),forward_traj[m].s);
    }    
  }
  

  // position at the end of the swim
  z_=z;
  x_=S(state_x,0);
  y_=S(state_y,0);

  return NOERROR;
}


// Step the state vector through the field from oldz to newz.
// Uses the 4th-order Runga-Kutte algorithm.
double DKalmanFilter::Step(double oldz,double newz, double dEdx,DMatrix &S){
  double delta_z=newz-oldz;
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);

  double s=sqrt(1.+S(state_tx,0)*S(state_tx,0)+S(state_ty,0)*S(state_ty,0))
    *delta_z;
  CalcDeriv(oldz,delta_z,S,dEdx,D1);
  CalcDeriv(oldz+delta_z/2.,delta_z/2.,S+0.5*delta_z*D1,dEdx,D2);
  CalcDeriv(oldz+delta_z/2.,delta_z/2.,S+0.5*delta_z*D2,dEdx,D3);
  CalcDeriv(oldz+delta_z,delta_z,S+delta_z*D3,dEdx,D4);
	
  S+=delta_z*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);

  return s;
}

// Step the state vector through the magnetic field and compute the Jacobian
// matrix.  Uses the 4th-order Runga-Kutte algorithm.
jerror_t DKalmanFilter::StepJacobian(double oldz,double newz,DMatrix S,
				   double dEdx,DMatrix &J){
   // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix J1(5,5),J2(5,5),J3(5,5),J4(5,5);  
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  
  double delta_z=newz-oldz;
  CalcDerivAndJacobian(oldz,delta_z,S,dEdx,J1,D1);
  CalcDerivAndJacobian(oldz+delta_z/2.,delta_z/2.,S+0.5*delta_z*D1,dEdx,J2,D2);
  J2=J2+0.5*(J2*J1);
  CalcDerivAndJacobian(oldz+delta_z/2.,delta_z/2.,S+0.5*delta_z*D2,dEdx,J3,D3);
  J3=J3+0.5*(J3*J2);
  CalcDerivAndJacobian(oldz+delta_z,delta_z,S+delta_z*D3,dEdx,J4,D4);
  J4=J4+J4*J3;

  J+=delta_z*(ONE_SIXTH*J1+ONE_THIRD*J2+ONE_THIRD*J3+ONE_SIXTH*J4);
  
  return NOERROR;
}

// Calculate the derivative for the alternate set of parameters {q/pT, phi, 
// tan(lambda),D,z}
jerror_t DKalmanFilter::CalcDeriv(double ds,DVector3 pos,DVector3 &dpos,
				  DVector3 B,DMatrix S,double dEdx,
				  DMatrix &D1){
   //Direction at current point
  double tanl=S(state_tanl,0);
  double phi=S(state_phi,0);
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  // Other parameters
  double q_over_pt=S(state_q_over_pt,0);
  double pt=fabs(1./q_over_pt);

  // Derivative of S with respect to s
  double temp1=B.y()*cosphi-B.x()*sinphi;
  D1(state_q_over_pt,0)
    =qBr2p*q_over_pt*q_over_pt*sinl*temp1;
  if (fabs(dEdx)>0){    
    double p=pt/cosl;
    double E=sqrt(p*p+MASS*MASS);
    D1(state_q_over_pt,0)+=-q_over_pt*E/p/p*dEdx;
  }
  D1(state_phi,0)
    =qBr2p*q_over_pt*(B.x()*cosphi*sinl+B.y()*sinphi*sinl-B.z()*cosl);
  D1(state_tanl,0)=qBr2p*q_over_pt*temp1/cosl;
  D1(state_z,0)=sinl;
  // Second order correction
  D1(state_z,0)+=qBr2p*q_over_pt*ds/2.*cosl*cosl*temp1;
 
  // New direction
  dpos.SetXYZ(cosl*cosphi,cosl*sinphi,sinl);

  // second order correction
  dpos(0)+=qBr2p*q_over_pt*cosl*ds/2.*(B.z()*cosl*sinphi-B.y()*sinl);
  dpos(1)+=qBr2p*q_over_pt*cosl*ds/2.*(B.x()*sinl-B.z()*cosl*cosphi);
  dpos(2)+=qBr2p*q_over_pt*ds/2.*cosl*cosl*temp1;

  return NOERROR;
}

// Calculate the derivative and Jacobian matrices for the alternate set of 
// parameters {q/pT, phi, tan(lambda),D,z}
jerror_t DKalmanFilter::CalcDerivAndJacobian(double ds,DVector3 pos,
					     DVector3 &dpos,
					     DVector3 B,
					     DMatrix S,double dEdx,
					     DMatrix &J1,DMatrix &D1){  
  //Direction at current point
  double tanl=S(state_tanl,0);
  double phi=S(state_phi,0);
  double cosphi=cos(phi);
  double sinphi=sin(phi);
  double lambda=atan(tanl);
  double sinl=sin(lambda);
  double cosl=cos(lambda);
  // Other parameters
  double q_over_pt=S(state_q_over_pt,0);
  double pt=fabs(1./q_over_pt);
  double q=pt*q_over_pt;
 
  //field gradient at (x,y,z)
  double dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bfield->GetFieldGradient(pos.x(),pos.y(),pos.z(),dBxdx,dBxdy,dBxdz,dBydx,
			   dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  // Derivative of S with respect to s
  double temp1=B.y()*cosphi-B.x()*sinphi;
  D1(state_q_over_pt,0)=qBr2p*q_over_pt*q_over_pt*sinl*temp1;
  D1(state_phi,0)=qBr2p*q_over_pt*(B.x()*cosphi*sinl+B.y()*sinphi*sinl-B.z()*cosl);
  D1(state_tanl,0)=qBr2p*q_over_pt*temp1/cosl;
  D1(state_z,0)=sinl+q_over_pt*qBr2p*ds/2.*cosl*cosl*temp1;

  // New direction
  dpos.SetXYZ(cosl*cosphi
	      +q_over_pt*cosl*qBr2p*ds/2.*(B.z()*cosl*sinphi-B.y()*sinl),
	      cosl*sinphi
	      +q_over_pt*cosl*qBr2p*ds/2.*(B.x()*sinl-B.z()*cosl*cosphi),
	      D1(state_z,0));
 
  // Jacobian matrix elements
  J1(state_phi,state_phi)=qBr2p*q_over_pt*sinl*(B.y()*cosphi-B.x()*sinphi);
  J1(state_phi,state_q_over_pt)=qBr2p*(B.x()*cosphi*sinl+B.y()*sinphi*sinl
				       -B.z()*cosl);
  J1(state_phi,state_tanl)=qBr2p*q_over_pt*(B.x()*cosphi*cosl+B.y()*sinphi*cosl
					    +B.z()*sinl)/(1.+tanl*tanl);
  J1(state_phi,state_z)
    =qBr2p*q_over_pt*(dBxdz*cosphi*sinl+dBydz*sinphi*sinl-dBzdz*cosl);

  J1(state_tanl,state_phi)=-qBr2p*q_over_pt*(B.y()*sinphi+B.x()*cosphi)/cosl;
  J1(state_tanl,state_q_over_pt)=qBr2p*(B.y()*cosphi-B.x()*sinphi)/cosl;
  J1(state_tanl,state_tanl)=qBr2p*q_over_pt*sinl*(B.y()*cosphi-B.x()*sinphi);
  J1(state_tanl,state_z)=qBr2p*q_over_pt*(dBydz*cosphi-dBxdz*sinphi)/cosl;
  
  J1(state_q_over_pt,state_phi)
    =-qBr2p*q_over_pt*q_over_pt*sinl*(B.y()*sinphi+B.x()*cosphi);

  double temp=sqrt(1.+cosl*cosl*q_over_pt*q_over_pt*MASS*MASS);
  J1(state_q_over_pt,state_q_over_pt)
    =2.*qBr2p*q_over_pt*sinl*(B.y()*cosphi-B.x()*sinphi);

  if (fabs(dEdx)>0){  
    double p=pt/cosl;
    double E=sqrt(p*p+MASS*MASS);
    D1(state_q_over_pt,0)+=-q_over_pt*E/p/p*dEdx;
    J1(state_q_over_pt,state_q_over_pt)+=-q*dEdx*cosl*q_over_pt*(2.+3.*cosl*cosl*q_over_pt*q_over_pt*MASS*MASS)
      /temp;
  }
  J1(state_q_over_pt,state_tanl)
    =qBr2p*q_over_pt*q_over_pt*cosl*cosl*cosl*(B.y()*cosphi-B.x()*sinphi)
    +q*dEdx*sinl*cosl*cosl/temp*q_over_pt*q_over_pt
    *(1.+2.*cosl*cosl*MASS*MASS*q_over_pt*q_over_pt);
  J1(state_q_over_pt,state_z)
    =qBr2p*q_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);

  J1(state_z,state_tanl)=cosl/(1.+tanl*tanl);

  return NOERROR;
}

// Convert between the forward parameter set {x,y,tx,ty,q/p} and the central
// parameter set {q/pT,phi,tan(lambda),D,z}
jerror_t DKalmanFilter::ConvertStateVector(double z,double wire_x, 
                                           double wire_y,
                                           DMatrix S, DMatrix C,
                                           DMatrix &Sc, DMatrix &Cc){
  double x=S(state_x,0),y=S(state_y,0);
  double tx=S(state_tx,0),ty=S(state_ty,0),q_over_p=S(state_q_over_p,0);
  double factor=1./sqrt(1.+tx*tx+ty*ty);
  double tanl=1./sqrt(tx*tx+ty*ty);
  double cosl=cos(atan(tanl));
  Sc(state_q_over_pt,0)=q_over_p/cosl;
  Sc(state_phi,0)=atan2(ty,tx);
  Sc(state_tanl,0)=tanl;
  //  Sc(state_D,0)=sqrt((x-wire_x)*(x-wire_x)+(y-wire_y)*(y-wire_y));
  Sc(state_D,0)=sqrt(x*x+y*y);
  Sc(state_z,0)=z;

  DMatrix J(5,5);
  J(state_tanl,state_tx)=-tx*tanl*tanl*tanl;
  J(state_tanl,state_ty)=-ty*tanl*tanl*tanl;
  J(state_z,state_x)=1./tx;
  J(state_z,state_y)=1./ty;
  J(state_q_over_pt,state_q_over_p)=1./cosl;
  J(state_q_over_pt,state_tx)=-tx*q_over_p*tanl*tanl*tanl*factor;
  J(state_q_over_pt,state_ty)=-ty*q_over_p*tanl*tanl*tanl*factor;
  J(state_phi,state_tx)=-ty/(tx*tx+ty*ty);
  J(state_phi,state_ty)=tx/(tx*tx+ty*ty);
  //J(state_D,state_x)=(x-wire_x)/S(state_D,0);
  //J(state_D,state_y)=(y-wire_y)/S(state_D,0);
  J(state_D,state_x)=x/Sc(state_D,0);
  J(state_D,state_y)=y/Sc(state_D,0);

  Cc=J*(C*DMatrix(DMatrix::kTransposed,J));

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DKalmanFilter::FixedStep(DVector3 &pos,double ds,DMatrix &S,
				  double dEdx){
  // Matrices for intermediate steps
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  DMatrix S1(5,1),S2(5,1),S3(5,1),S4(5,1);
  DVector3 dpos1,dpos2,dpos3,dpos4;
  
  // Magnetic field
  DVector3 B;
  bfield->GetField(pos.x(),pos.y(),pos.z(), B(0), B(1), B(2));

  CalcDeriv(0.,pos,dpos1,B,S,dEdx,D1);

  DVector3 mypos=pos+(ds/2.)*dpos1;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(), B(0), B(1), B(2));
  S1=S+(0.5*ds)*D1; 

  CalcDeriv(ds/2.,mypos,dpos2,B,S1,dEdx,D2);

  mypos=pos+(ds/2.)*dpos2;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(), B(0), B(1), B(2));
  S2=S+(0.5*ds)*D2; 

  CalcDeriv(ds/2.,mypos,dpos3,B,S2,dEdx,D3);

  mypos=pos+ds*dpos3;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(), B(0), B(1), B(2));
  S3=S+ds*D3;

  CalcDeriv(ds,mypos,dpos4,B,S3,dEdx,D4);
   
  // New state vector
  S+=ds*(ONE_SIXTH*D1+ONE_THIRD*D2+ONE_THIRD*D3+ONE_SIXTH*D4);

  // New position
  pos+=
    (dpos1=ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4));

  return NOERROR;
}


jerror_t DKalmanFilter::CalcDirMom(double ds,double dEdx,double q,DVector3 mom,
				   DVector3 pos,DVector3 &dmom,DVector3 &dpos){
  DVector3 pdir=mom.Unit();

  //B-field and field gradient at (x,y,z)
  double Bx,By,Bz;
  bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
  
  // Change in momentum
  dmom.SetXYZ(q*qBr2p*(Bz*pdir.y()-By*pdir.z()),
	      q*qBr2p*(Bx*pdir.z()-Bz*pdir.x()),
	      q*qBr2p*(By*pdir.x()-Bx*pdir.y()));  
  
  // change in position 
  dpos.SetXYZ(pdir.x()
	      +q*qBr2p*ds/2.*(Bz*pdir.y()-By*pdir.z()),
	      pdir.y()
	      +q*qBr2p*ds/2.*(Bx*pdir.z()-Bz*pdir.x()),
	      pdir.z()
	      +q*qBr2p*ds/2.*(By*pdir.x()-Bx*pdir.y())
	      );

  return NOERROR;
}



// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
double DKalmanFilter::Step(DVector3 &pos, DVector3 wire_orig,DVector3 wiredir,
				     double ds,DMatrix &S,double dEdx){
  DVector3 dpos1,dpos2,dpos3,dpos4;
  DVector3 dmom1,dmom2,dmom3,dmom4;

  //Direction at current point
  double cosphi=cos(S(state_phi,0));
  double sinphi=sin(S(state_phi,0));
  // Other parameters
  double q=S(state_q_over_pt,0)>0?1.:-1.;
  double pt=fabs(1./S(state_q_over_pt,0));
  DVector3 mom(pt*cosphi,pt*sinphi,pt*S(state_tanl,0));
  DVector3 pdir=mom.Unit();

  /*
  CalcDirMom(0.,dEdx,q,mom,pos,dmom1,dpos1);
  CalcDirMom(ds/2.,dEdx,q,mom+(ds/2.)*dmom1,
	     pos+ds/2.*pdir+ds*ds/8.*dmom1,
	     dmom2,dpos2);
  CalcDirMom(ds/2.,dEdx,q,mom+(ds/2.)*dmom2,
	     pos+ds/2*pdir+ds*ds/8.*dmom1,
	     dmom3,dpos3);
  CalcDirMom(ds,dEdx,q,mom+ds*dmom3,
	     pos+ds*pdir+ds*ds/2.*dmom3,
	     dmom4,dpos4);*/
  double VECT[7], VOUT[7];
  VECT[0] = pos.x();
  VECT[1] = pos.y();
  VECT[2] = pos.z();
  VECT[6] = mom.Mag();
  VECT[3] = fabs(1./S(state_q_over_pt,0))*cos(S(state_phi,0))/VECT[6];
  VECT[4] =  fabs(1./S(state_q_over_pt,0))*sin(S(state_phi,0))/VECT[6];
  VECT[5] = fabs(1./S(state_q_over_pt,0))*S(state_tanl,0)/VECT[6];
  
  grkuta_(&q, &ds, VECT, VOUT, bfield);

  pos.SetXYZ(VOUT[0],VOUT[1],VOUT[2]);
  mom.SetXYZ(VOUT[3]*VECT[6],VOUT[4]*VECT[6],VOUT[5]*VECT[6]);
  
  // New state vector
  pt=mom.Perp();
  S(state_q_over_pt,0)=q/pt;
  S(state_tanl,0)=mom.z()/pt;
  S(state_phi,0)=mom.Phi();
  S(state_z,0)=pos.z();
  
  return ds;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DKalmanFilter::StepJacobian(DVector3 pos, DVector3 wire_orig,
				     DVector3 wiredir,
				     double ds,DMatrix S,
				     double dEdx,DMatrix &J){
  // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix J1(5,5),J2(5,5),J3(5,5),J4(5,5);  
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  DMatrix S1(5,1),S2(5,1),S3(5,1),S4(5,1);
  DVector3 dpos1,dpos2,dpos3,dpos4;

   // charge
  double q=(S(state_q_over_pt,0)>0)?1.:-1.;
  // Magnetic field
  DVector3 B;
  bfield->GetField(pos.x(),pos.y(),pos.z(), B(0), B(1), B(2));
  double Bz_=fabs(B.z());
  double qpt=1./S(state_q_over_pt,0);
  double sinphi=sin(S(state_phi,0));
  double cosphi=cos(S(state_phi,0));
  double D=S(state_D,0);

  CalcDerivAndJacobian(0.,pos,dpos1,B,S,dEdx,J1,D1);

  DVector3 mypos=pos+(ds/2.)*dpos1;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(), B(0), B(1), B(2));
  S1=S+(0.5*ds)*D1; 

  CalcDerivAndJacobian(ds/2.,mypos,dpos2,B,S1,dEdx,J2,D2);
  J2=J2+0.5*(J2*J1);

  mypos=pos+(ds/2.)*dpos2;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(), B(0), B(1), B(2));
  S2=S+(0.5*ds)*D2;

  CalcDerivAndJacobian(ds/2.,mypos,dpos3,B,S2,dEdx,J3,D3);
  J3=J3+0.5*(J3*J2);  

  mypos=pos+ds*dpos3;
  bfield->GetField(mypos.x(),mypos.y(),mypos.z(), B(0), B(1), B(2));
  S3=S+ds*D3;

  CalcDerivAndJacobian(ds,mypos,dpos4,B,S3,dEdx,J4,D4);
  J4=J4+J4*J3;

  // New Jacobian matrix
  J+=ds*(ONE_SIXTH*J1+ONE_THIRD*J2+ONE_THIRD*J3+ONE_SIXTH*J4);

  // change in position
  DVector3 dpos
    =ds*(ONE_SIXTH*dpos1+ONE_THIRD*dpos2+ONE_THIRD*dpos3+ONE_SIXTH*dpos4);
  double rc=sqrt(dpos.Perp2()
		 +2.*(D+qpt/qBr2p/Bz_)*(dpos.x()*sinphi-dpos.y()*cosphi)
		 +(D+qpt/qBr2p/Bz_)*(D+qpt/qBr2p/Bz_));
    
  J(state_D,state_D)=q*(dpos.x()*sinphi-dpos.y()*cosphi+D+qpt/qBr2p/Bz_)/rc;
  J(state_D,state_q_over_pt)=qpt*qpt/qBr2p/Bz_*(J(state_D,state_D)-1.);
  J(state_D,state_phi)=q*(D+qpt/qBr2p/Bz_)*(dpos.x()*cosphi+dpos.y()*sinphi)/rc;
  
  return NOERROR;
}

// Multiple scattering covariance matrix for central track parameters
jerror_t DKalmanFilter::GetProcessNoiseCentral(double mass_hyp,double ds,
					       DVector3 pos,double X0,
					       DMatrix Sc,DMatrix &Q){
  Q.Zero();
  if (isfinite(X0) && X0<1e8 && X0>0.){
    DMatrix Q1(5,5);
    double tanl=Sc(state_tanl,0);
    double q_over_pt=Sc(state_q_over_pt,0); 
    //  double X0=material->GetRadLen(pos.x(),pos.y(),pos.z());
    double my_ds=fabs(ds);
    
    Q1(state_phi,state_phi)=1.+tanl*tanl;
    Q1(state_tanl,state_tanl)=(1.+tanl*tanl)*(1.+tanl*tanl);
    Q1(state_q_over_pt,state_q_over_pt)=q_over_pt*q_over_pt*tanl*tanl;
    Q1(state_q_over_pt,state_tanl)=Q1(state_tanl,state_q_over_pt)
      =q_over_pt*tanl*(1.+tanl*tanl);
    Q1(state_D,state_D)=ds*ds/3.;
    Q1(state_D,state_tanl)=Q1(state_tanl,state_D)
      //=-my_ds/2.*cos(Sc(state_phi,0))*(1.+tanl*tanl);
      =my_ds/2.*(1.+tanl*tanl);
    Q1(state_D,state_phi)=Q1(state_phi,state_D)
      //my_ds/2.*sin(Sc(state_phi,0))*sqrt(1.+tanl*tanl);
      =my_ds/2.*sqrt(1.+tanl*tanl);
    Q1(state_D,state_q_over_pt)=Q1(state_q_over_pt,state_D)
      //      -my_ds/2.*cos(Sc(state_phi,0))*tanl;
      =my_ds/2.*q_over_pt*tanl;

    double p2=(1.+tanl*tanl)/q_over_pt/q_over_pt;
    double sig2_ms=0.0136*0.0136*(1.+mass_hyp*mass_hyp/p2)*my_ds/X0/p2
      *(1.+0.038*log(my_ds/X0*(1.+mass_hyp*mass_hyp/p2)))
      *(1.+0.038*log(my_ds/X0*(1.+mass_hyp*mass_hyp/p2)));
    //sig2_ms=0.;

    Q=sig2_ms*Q1;
  }
  //Q.Print();

  return NOERROR;

}


// Compute contributions to the covariance matrix due to multiple scattering
jerror_t DKalmanFilter::GetProcessNoise(double mass_hyp,double ds,double z,
					double X0,DMatrix S,DMatrix &Q){
  Q.Zero();
  if (isfinite(X0) && X0<1e8 && X0>0.){
    DMatrix Q1(5,5);
    double tx=S(state_tx,0),ty=S(state_ty,0);
    double one_over_p_sq=S(state_q_over_p,0)*S(state_q_over_p,0);
    double my_ds=fabs(ds);
    // double X0=material->GetRadLen(S(state_x,0),S(state_y,0),z);
    
    Q1(state_tx,state_tx)=(1.+tx*tx)*(1.+tx*tx+ty*ty);
    Q1(state_ty,state_ty)=(1.+ty*ty)*(1.+tx*tx+ty*ty);
    Q1(state_tx,state_ty)=Q1(state_ty,state_tx)=tx*ty*(1.+tx*tx+ty*ty);
    Q1(state_x,state_x)=ds*ds/3.;
    Q1(state_y,state_y)=ds*ds/3.;
    Q1(state_y,state_ty)=Q1(state_ty,state_y)
      //      =my_ds/2.*tx*(1.+tx*tx+ty*ty)/sqrt(tx*tx+ty*ty);
      = my_ds/2.*sqrt((1.+tx*tx+ty*ty)*(1.+ty*ty));
    Q1(state_x,state_tx)=Q1(state_tx,state_x)
      // =my_ds/2.*ty*sqrt((1.+tx*tx+ty*ty)/(tx*tx+ty*ty));
      = my_ds/2.*sqrt((1.+tx*tx+ty*ty)*(1.+tx*tx));
    
    double sig2_ms= 0.0136*0.0136*(1.+one_over_p_sq*mass_hyp*mass_hyp)
      *one_over_p_sq*my_ds/X0
      *(1.+0.038*log(my_ds/X0*(1.+one_over_p_sq*mass_hyp*mass_hyp)))
      *(1.+0.038*log(my_ds/X0*(1.+one_over_p_sq*mass_hyp*mass_hyp)));
	//  sig2_ms=0.;   
    Q=sig2_ms*Q1;
  }
  return NOERROR;
}

// Calculate the energy loss per unit length given properties of the material
// through which a particle of momentum p is passing
double DKalmanFilter::GetdEdx(double M,double q_over_p,double Z,
			      double A, double rho){
  if (rho<=0.) return 0.;
  //return 0;
  double betagamma=1./M/fabs(q_over_p);
  double beta2=1./(1.+M*M*q_over_p*q_over_p);
  if (beta2<EPS) return 0.;

  double Me=0.000511; //GeV
  double m_ratio=Me/M;
  double Tmax
    =2.*Me*betagamma*betagamma/(1.+2.*sqrt(1.+betagamma*betagamma)*m_ratio
				+m_ratio*m_ratio);
  double I0=12.*Z+7.; //eV
  if (Z>12) I0=9.76*Z+58.8*pow(Z,-0.19);

  // The next block of code approximates the "density effect" using the 
  // prescription employed by GEANT 3.21.
  double delta=0.;  //density effect  
  double X=log10(betagamma);
  double X0,X1;
  double C=-2.*log(I0/28.816/sqrt(rho*Z/A))-1.;
  double Cbar=-C;
  if (rho>0.01){ // not a gas
    if (I0<100){
      if (Cbar<=3.681) X0=0.2;
      else X0=-0.326*C-1.;
      X1=2.;
    }
    else{
      if (Cbar<=5.215) X0=0.2;
      else X0=0.326*Cbar-1.5;
      X1=3.;
    }
  }
  else{ // gases
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
    delta=4.606*X+C+(Cbar-4.606*X0)*pow(X1-X,3.)/pow(X1-X0,3.);
  else if (X>=X1)
    delta=4.606*X+C;

  I0*=1e-9; // convert to GeV
  return -0.0001535*Z/A*rho/beta2*(log(2.*Me*betagamma*betagamma*Tmax/I0/I0)
				   -2.*beta2-delta);
}

// Calculate the variance in the energy loss in a Gaussian approximation
double DKalmanFilter::GetEnergyVariance(double ds,double q_over_p,double Z,
			      double A, double rho){
  if (rho<=0.) return 0.;
  //return 0;
  double betagamma=1./MASS/fabs(q_over_p);
  double beta2=1./(1.+MASS*MASS*q_over_p*q_over_p);
  if (beta2<EPS) return 0.;

  double Me=0.000511; //GeV
  double m_ratio=Me/MASS;
  double Tmax
    =2.*Me*betagamma*betagamma/(1.+2.*sqrt(1.+betagamma*betagamma)*m_ratio
				+m_ratio*m_ratio);

  return 0.0001535*fabs(ds)*(rho*Z/A)*Tmax/beta2;
}


// Swim the state vector through the field from the start of the reference
// trajectory to the end
jerror_t DKalmanFilter::SwimToPlane(DMatrix &S){
  int max=forward_traj_cdc.size();
  double z,newz=0.,dedx=0.;
  
  // If we have trajectory entries for the CDC, start there
  if (max>0){
    max--;
    z=forward_traj_cdc[max].pos.Z();
    for (unsigned int m=max-1;m>0;m--){
      newz=forward_traj_cdc[m].pos.Z();  
      if (do_energy_loss){
	dedx=GetdEdx(MASS,S(state_q_over_p,0),forward_traj_cdc[m].Z,
		     forward_traj_cdc[m].A,forward_traj_cdc[m].density);
      }
      Step(z,newz,dedx,S);  
      z=newz;
    }
    if (do_energy_loss){
      dedx=GetdEdx(MASS,S(state_q_over_p,0),forward_traj_cdc[1].Z,
		   forward_traj_cdc[1].A,forward_traj_cdc[1].density);
    }
    newz=forward_traj_cdc[0].pos.Z(); 
    Step(z,newz,dedx,S);  
  }
  // Follow track into FDC
  max=forward_traj.size()-1;
  if (max>0){
    z=forward_traj[max].pos.Z(); 
    for (unsigned int m=max-1;m>0;m--){
      // newz=z+STEP_SIZE;  
    newz=forward_traj[m].pos.z();
    if (do_energy_loss){
      dedx=GetdEdx(MASS,S(state_q_over_p,0),forward_traj[m].Z,
		   forward_traj[m].A,forward_traj[m].density);
    }
    Step(z,newz,dedx,S); 
    z=newz;
    } 
    if (do_energy_loss){
      dedx=GetdEdx(MASS,S(state_q_over_p,0),forward_traj[1].Z,
		   forward_traj[1].A,forward_traj[1].density);
    }
    newz=forward_traj[0].pos.Z(); 
    Step(z,newz,dedx,S);  
  }
  z_=newz;
  
  return NOERROR;
}

// Swim the state vector from z_start to z_end through the field
jerror_t DKalmanFilter::SwimToPlane(double z_start,double z_end, DMatrix &S){
  int num_inc=(int)((z_end-z_start)/STEP_SIZE);
  double z=z_start;
  double dedx=0.;
   
  for (int i=0;i<num_inc;i++){
    z_=z+STEP_SIZE;  
    Step(z,z_,dedx,S);  
    z=z_;
  } 

  // Final step
  Step(z_,z_end,dedx,S);  
  z_=z_end;

  return NOERROR;
}


// Swim the state vector and the covariance matrix from z_start to z_end 
// through the field
jerror_t DKalmanFilter::SwimToPlane(double z_start,double z_end, DMatrix &S,
				    DMatrix &C){
  DMatrix J(5,5),J_T(5,5);
  int num_inc=(int)((z_end-z_start)/0.5);
  double z=z_start;
  double dedx=0.;
  for (int i=0;i<num_inc;i++){
     z_=z+0.5;  
     StepJacobian(z,z_,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);     
    C=J*(C*J_T);
    z=z_;
  }
 
  // Final step
  StepJacobian(z_,z_end,S,dedx,J);  
  J_T=DMatrix(DMatrix::kTransposed,J);     
  C=J*(C*J_T);
  z_=z_end;

  return NOERROR;
}
// Swim the state vector and the covariance matrix from the current position 
// to the position corresponding to the radius R
jerror_t DKalmanFilter::SwimToRadius(DVector3 &pos,double Rf,DMatrix &Sc,
				     DMatrix &Cc){
  DMatrix Jc(5,5),Jc_T(5,5),Q(5,5);
  double R=pos.Perp();
  double ds=0.25;
  double dedx=0.;
  unsigned int m=cdchits.size()-1;
  while (R<Rf){ 
    // Get dEdx for this step
    //double tanl=Sc(state_tanl,0);
    //double cosl=cos(atan(tanl));
    //double q_over_p=Sc(state_q_over_pt,0)*cosl;
    //dedx=GetdEdx(M,q_over_p,18.,39.9,0.00166);

   // Step the position, state vector, and covariance matrix through the field
    StepJacobian(pos,cdchits[m]->origin,cdchits[m]->dir,ds,Sc,dedx,Jc);
    Cc=Jc*(Cc*DMatrix(DMatrix::kTransposed,Jc));

    FixedStep(pos,ds,Sc,dedx);

    // New radius
    R=pos.Perp();
  }

  if (pos.z()>0 && pos.z()<200.){
    double ds1 =0.5,s1=0.;
    // Save the current state vector and position
    DMatrix S0(5,1);
    DVector3 newpos=pos;
    S0=Sc;
    int num_iter=0;
    // We've gone too far, so we back track...
    while (fabs(R-Rf)>EPS && num_iter<10){
      ds1/=2.;
      if (R>Rf) ds=-ds1;
      else ds=ds1;
      FixedStep(newpos,ds,S0,dedx);
      R=newpos.Perp();
      s1+=ds;
      num_iter++;
    }
    // Final step
    StepJacobian(pos,cdchits[m]->origin,cdchits[m]->dir,s1,Sc,dedx,Jc);
    Cc=Jc*(Cc*DMatrix(DMatrix::kTransposed,Jc));

    FixedStep(pos,s1,Sc,dedx);
  }

  // Wire position for wire#0
  double x=cdchits[0]->origin.x()
    +(pos.z()-cdchits[0]->origin.z())*cdchits[0]->dir.x(); 
  double y=cdchits[0]->origin.y()
    +(pos.z()-cdchits[0]->origin.z())*cdchits[0]->dir.y();
  double dx=pos.x()-x;
  double dy=pos.y()-y;

  //B-field at (x,y,z)
  double Bx,By,Bz;
  bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  double Rc=1./qBr2p/B/Sc(state_q_over_pt,0);
  double cosphi=cos(Sc(state_phi,0));
  double sinphi=sin(Sc(state_phi,0));
  double xc=pos.x()-Rc*sinphi;
  double yc=pos.y()+Rc*cosphi;
  double rc2=(x-xc)*(x-xc)+(y-yc)*(y-yc);
  double sign=-Sc(state_q_over_pt,0)/fabs(Sc(state_q_over_pt,0));
  if (rc2<Rc*Rc) sign*=-1.;
  
  //sign=1.;

  // New doca
  Sc(state_D,0)=sign*sqrt(dx*dx+dy*dy);  

  // Update the internal variables
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  return NOERROR;
}

// Swim the state vector and the covariance matrix from the current position 
// to the position corresponding to the radius R
jerror_t DKalmanFilter::SwimToRadius(DVector3 &pos,double Rf,DMatrix &Sc){
  DMatrix Jc(5,5),Jc_T(5,5),Q(5,5);
  double R=pos.Perp();
  double ds=0.25;
  double dedx=0.;
  DVector3 origin=cdchits[0]->origin;
  DVector3 dir=cdchits[0]->dir;
  while (R<Rf){ 
    // Get dEdx for this step
    //double tanl=Sc(state_tanl,0);
    //double cosl=cos(atan(tanl));
    //double q_over_p=Sc(state_q_over_pt,0)*cosl;
    //dedx=GetdEdx(M,q_over_p,18.,39.9,0.00166);

    // Step the position, state vecto
    FixedStep(pos,ds,Sc,dedx);

    // New radius
    R=pos.Perp();
  }

  // Update the internal variables
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  return NOERROR;
}


// Interface routine for Kalman filter
jerror_t DKalmanFilter::KalmanLoop(double mass_hyp,int pass){
  if (z_<0) return VALUE_OUT_OF_RANGE;
  this->pass=pass;
  
  DMatrix S(5,1),C(5,5),Sc(5,1),Cc(5,5);
  DMatrix C0(5,5);
  double chisq=1.e16,chisq_forward=1.e16,chisq_central=1.e16;
  chisq_=1.e16;
  // position along track. 
  DVector3 pos(x_,y_,z_); 

  // energy loss flag
  //do_energy_loss=false;
  do_energy_loss=true;

  // multiple scattering flag
  //do_multiple_scattering=false;
  do_multiple_scattering=true;

  // Initialize path length variable
  //  len=0;

  // deal with hits in FDC
  if (fdchits.size()>0){   
    // Order the hits
    sort(fdchits.begin(),fdchits.end(),DKalmanFDCHit_cmp);
    if (cdchits.size()>0){
      // Order the CDC hits by ring number
      sort(cdchits.begin(),cdchits.end(),DKalmanCDCHit_cmp);
    }

    // Initialize the state vector and covariance matrix
    S(state_x,0)=x_;
    S(state_y,0)=y_;
    S(state_tx,0)=tx_;
    S(state_ty,0)=ty_;
    S(state_q_over_p,0)=q_over_p_; 

    // Initial guess for forward representation covariance matrix
    C0(state_x,state_x)=1.;
    C0(state_y,state_y)=1.;
    C0(state_tx,state_tx)=0.010*0.010;
    C0(state_ty,state_ty)=0.010*0.010;
    C0(state_q_over_p,state_q_over_p)=0.04*q_over_p_*q_over_p_;

    DMatrix Slast(S);
    DMatrix Clast(C0);
    
    double chisq_iter=chisq;
    double chisq_min=chisq;
    double zvertex=65.;
    double scale_factor=200.,anneal_factor=1.;
    // Iterate over reference trajectories
    for (int iter2=0;iter2<20;iter2++){   
      //for (int iter2=0;iter2<1;iter2++){   
      //if (iter2>0) do_energy_loss=true;
      //if (iter2>0) do_multiple_scattering=true;
      
      double f=2.5;
      if (pass==kTimeBased)
	anneal_factor=scale_factor/pow(f,iter2)+1.;

       // Initialize path length variable
      len=0;

      // If we have cdc hits, swim through the field past these measurements
      // first
      if (cdchits.size()>0){
	SetCDCForwardReferenceTrajectory(S);
      }
      
      // Swim once through the field out to the most upstream FDC hit
      SetReferenceTrajectory(S);
      //C0=C;

      //printf("forward iteration %d cdc size %d\n",iter2,forward_traj_cdc.size());
      
      if (forward_traj.size()> 0){
	unsigned int num_iter=NUM_ITER;
	//num_iter=1;
	//if (cdchits.size()==0) num_iter=3;
	chisq_forward=1.e16;
	for (unsigned int iter=0;iter<20;iter++) {      
	  // perform the kalman filter 
	  
	  if (iter>0){
	    // Swim back to the first (most downstream) plane and use the new 
	    // values of S and C as the seed data to the Kalman filter 
	    SwimToPlane(S);
	  } 
	  
	  C=C0;
	  KalmanForward(anneal_factor,S,C,chisq);
	  // KalmanForward(1.,S,C,chisq);
	  
	  //printf("forward chi2 %f p %f \n",chisq,1./S(state_q_over_p,0));
	  // include any hits from the CDC on the trajectory
	  if (cdchits.size()>0 && forward_traj_cdc.size()>0){
	    // Proceed into CDC
	    KalmanForwardCDC(mass_hyp,anneal_factor,S,C,chisq);
	    //KalmanForwardCDC(mass_hyp,1.,S,C,chisq);
	  }

	  //printf("iter %d chi2 %f %f\n",iter2,chisq,chisq_forward);
	  if (isnan(chisq) || 
	      (fabs(chisq-chisq_forward)<0.1 || chisq>chisq_forward))
	    break;
	  chisq_forward=chisq;
	  Slast=S;
	  Clast=C;
	} //iteration
      }  

      //printf("iter2: %d chi2 %f %f\n",iter2,chisq_forward,chisq_iter);
      
      // Abort loop if the chisq is not changing much or increasing too much
      if ( isnan(chisq_forward) || (//iter2>7 && 
	  (fabs(chisq_forward-chisq_iter)<0.1 
	   || chisq_forward-chisq_iter>10.))) 
	break;
      chisq_iter=chisq_forward;
      
      // Find the state at the so-called "vertex" position
      ExtrapolateToVertex(Slast,Clast);
      C=Clast;
      S=Slast;
      zvertex=z_;
    }

    // Convert from forward rep. to central rep.
    ConvertStateVector(zvertex,0.,0.,Slast,Clast,Sc,Cc);

    // Track Parameters at "vertex"
    phi_=Sc(state_phi,0);
    q_over_pt_=Sc(state_q_over_pt,0);
    tanl_=Sc(state_tanl,0);
    x_=Slast(state_x,0);
    y_=Slast(state_y,0);
    z_=zvertex;

    if (DEBUG_LEVEL>0)
      cout
	<< "Vertex:  p " 
	<<   1./q_over_pt_/cos(atan(tanl_))
	<< " theta "  << 90.0-180./M_PI*atan(tanl_)
	<< " vertex " << x_ << " " << y_ << " " << z_ <<endl;
    
    // Covariance matrices...
    // ... central parameterization
    vector<double>dummy;
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	dummy.push_back(Cc(i,j));
      }
      cov.push_back(dummy);
    }
    // ... forward parameterization
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
      dummy.push_back(Clast(i,j));
      }
      fcov.push_back(dummy);
    }

    // total chisq and ndf
    chisq_=chisq_iter;
    ndf=2*fdchits.size()+cdchits.size()-5;
        
    if (DEBUG_HISTS){
      TH2F *fdc_xresiduals=(TH2F*)gROOT->FindObject("fdc_xresiduals");
      if (fdc_xresiduals){
	for (unsigned int i=0;i<fdchits.size();i++){
	  fdc_xresiduals->Fill(fdchits[i]->z,fdchits[i]->xres);
	}
      }
      TH2F *fdc_yresiduals=(TH2F*)gROOT->FindObject("fdc_yresiduals");
      if (fdc_yresiduals){
	for (unsigned int i=0;i<fdchits.size();i++){
	  fdc_yresiduals->Fill(fdchits[i]->z,fdchits[i]->yres);
	}
      }
      /*
	TH2F *fdc_ypulls=(TH2F*)gROOT->FindObject("fdc_ypulls");
	if (fdc_ypulls) fdc_ypulls->Fill(fdchits[id]->z,R(state_y,0)/sqrt(RC(1,1)));
	
	TH2F *fdc_xpulls=(TH2F*)gROOT->FindObject("fdc_xpulls");
	if (fdc_xpulls) fdc_xpulls->Fill(fdchits[id]->z,R(state_y,0)/sqrt(RC(0,0)));
      */
    }
       
    return NOERROR;
  }

  if (cdchits.size()>0 && tanl_>0.84){
    sort(fdchits.begin(),fdchits.end(),DKalmanFDCHit_cmp);
    if (cdchits.size()>0){
      // Order the CDC hits by ring number
      sort(cdchits.begin(),cdchits.end(),DKalmanCDCHit_cmp);
    }

    // Initialize the state vector and covariance matrix
    S(state_x,0)=x_;
    S(state_y,0)=y_;
    S(state_tx,0)=tx_;
    S(state_ty,0)=ty_;
    S(state_q_over_p,0)=q_over_p_; 

    // Initial guess for forward representation covariance matrix
    C0(state_x,state_x)=1;
    C0(state_y,state_y)=1;
    C0(state_tx,state_tx)=0.010*0.010;
    C0(state_ty,state_ty)=0.010*0.010;
    C0(state_q_over_p,state_q_over_p)=0.04*q_over_p_*q_over_p_;
    
    DMatrix Slast(S);
    DMatrix Clast(C0);
    
    double chisq_iter=chisq;
    double chisq_min=chisq;
    double zvertex=65.;
    double scale_factor=200.,anneal_factor=1.;
    // Iterate over reference trajectories
    for (int iter2=0;iter2<20;iter2++){   
    //for (int iter2=0;iter2<1;iter2++){   
      //if (iter2>0) do_energy_loss=true;
      //if (iter2>0) do_multiple_scattering=true;
      
      if (pass==kTimeBased){
	double f=1.75;
	anneal_factor=scale_factor/pow(f,iter2)+1.;
      }
  
      // Initialize path length variable
      len=0;
            
      SetCDCForwardReferenceTrajectory(S);
      
      if (forward_traj_cdc.size()> 0){
	unsigned int num_iter=NUM_ITER;
	//num_iter=1;
	//if (cdchits.size()==0) num_iter=3;
	chisq_forward=1.e16;
	for (unsigned int iter=0;iter<20;iter++) {      
	  // perform the kalman filter 
	  
	  if (iter>0){
	    // Swim back to the first (most downstream) plane and use the new 
	    // values of S and C as the seed data to the Kalman filter 
	    SwimToPlane(S);
	  } 
	  
	  C=C0;
	  chisq=0.;
	  KalmanForwardCDC(mass_hyp,anneal_factor,S,C,chisq);
	  
	  //printf("iter %d chi2 %f %f\n",iter,chisq,chisq_forward);
	  if (isnan(chisq) || 
	      (fabs(chisq-chisq_forward)<0.1 || chisq>chisq_forward))
	    break;
	  chisq_forward=chisq;
	  Slast=S;
	  Clast=C;
	} //iteration
      }  
      
      //      printf("iter2: %d factor %f chi2 %f %f\n",iter2,anneal_factor,chisq_forward,chisq_iter);
      
      // Abort loop if the chisq is not changing much or increasing too much
      if ( isnan(chisq_forward) || (//iter2>12 && 
				    (fabs(chisq_forward-chisq_iter)<0.1 
	   || chisq_forward-chisq_iter>10.))) 
	break;
      chisq_iter=chisq_forward;
      
      // Find the state at the so-called "vertex" position
      ExtrapolateToVertex(Slast,Clast);
      C=Clast;
      S=Slast;
      zvertex=z_;
    }

    // Convert from forward rep. to central rep.
    ConvertStateVector(zvertex,0.,0.,Slast,Clast,Sc,Cc);
    
    // Track Parameters at "vertex"
    phi_=Sc(state_phi,0);
    q_over_pt_=Sc(state_q_over_pt,0);
    tanl_=Sc(state_tanl,0);
    x_=Slast(state_x,0);
    y_=Slast(state_y,0);
    z_=zvertex;

    if (DEBUG_LEVEL>0)
      cout
	<< "Vertex:  p " 
	<<   1./q_over_pt_/cos(atan(tanl_))
	<< " theta "  << 90.0-180./M_PI*atan(tanl_)
	<< " vertex " << x_ << " " << y_ << " " << z_ <<endl;
    
    // Covariance matrices...
    // ... central parameterization
    vector<double>dummy;
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	dummy.push_back(Cc(i,j));
      }
      cov.push_back(dummy);
    }

    // total chisq and ndf
    chisq_=chisq_iter;
    ndf=cdchits.size()-5;

    return NOERROR;
  }
  
  // Fit in Central region:  deal with hits in the CDC 
  if (cdchits.size()>0){  
    // Order the CDC hits by radius
    sort(cdchits.begin(),cdchits.end(),DKalmanCDCHit_cmp);

    // Initialize the state vector and covariance matrix
    Sc(state_q_over_pt,0)=q_over_pt_;
    Sc(state_phi,0)=phi_;
    Sc(state_tanl,0)=tanl_;
    Sc(state_z,0)=z_;  

    //C0(state_z,state_z)=1.;
    C0(state_z,state_z)=2.0;
    C0(state_q_over_pt,state_q_over_pt)=0.20*0.20*q_over_pt_*q_over_pt_;
    C0(state_phi,state_phi)=0.01*0.01;
    C0(state_D,state_D)=1.0;
    double dlambda=0.05;
    //dlambda=0.1;
    C0(state_tanl,state_tanl)=(1.+tanl_*tanl_)*(1.+tanl_*tanl_)
      *dlambda*dlambda;

    // Initialization
    Cc=C0;
    DMatrix Scbest(Sc),Sclast(Sc);
    DMatrix Ccbest(Cc),Cclast(Cc);
    DVector3 pos0=pos;
    DVector3 best_pos=pos;
  
    // iteration 
    double scale_factor=1., anneal_factor=1.;
    double chisq_iter=chisq;
    double chisq_min=chisq;
    for (int iter2=0;iter2<20;iter2++){  
    //for (int iter2=0;iter2<1;iter2++){  
      //if (iter2>0) do_energy_loss=true;
      //if (iter2>0) do_multiple_scattering=true;

      // Initialize path length variable
      len=0;
      
      // Abort if the chisq for the previous iteration is junk
      if (chisq_central==0.) break;

      //      printf("-------------central iteration %d cdc size %d\n",iter2,central_traj.size());

      // Calculate an annealing factor for the measurement errors that depends 
      // on the iteration,so that we approach the "true' measurement errors
      // by the last iteration.
      double f=3.5;
      if (pass==kTimeBased)
	anneal_factor=scale_factor/pow(f,iter2)+1.;

      // Initialize trajectory deque and position
      SetCDCReferenceTrajectory(pos0,Sc);
              
      if (central_traj.size()==0){
	cout << "No CDC reference trajectory????" <<endl;
	return RESOURCE_UNAVAILABLE;
      }

      // Iteration for given reference trajectory 
      chisq=1.e16;
      for (int iter=0;iter<20;iter++){
	Cc=C0;
	if (iter>0){
	  SwimCentral(pos,Sc);
	  //anneal_factor=scale_factor/pow(f,iter)+1.;
	}
	//anneal_factor=1.;

	jerror_t error=NOERROR;
	error=KalmanCentral(mass_hyp,anneal_factor,Sc,Cc,pos,chisq_central);
	if (error!=NOERROR) break;
	if (chisq_central==0.) break;
      	
	//fom=anneal_factor*chisq_central;
	if (chisq_central>=1e16 ){
	  cout 
	    << "-- central fit failed --" <<endl;
	  return VALUE_OUT_OF_RANGE;
	}
	
	if (DEBUG_LEVEL>0)
	  cout 
	    << "iteration " << iter+1  << " factor " << anneal_factor 
	    << " chi2 " 
	    << chisq_central << " p " 
	    <<   1./Sc(state_q_over_pt,0)/cos(atan(Sc(state_tanl,0)))
	    << " theta "  << 90.-180./M_PI*atan(Sc(state_tanl,0)) 
	    << " vertex " << x_ << " " << y_ << " " << z_ <<endl;
	
	if (isnan(chisq_central) || (fabs(chisq_central-chisq)<0.1 
				     || (chisq_central>chisq )))
	  break; 
	// Save the current "best" state vector and covariance matrix
	Cclast=Cc;
	Sclast=Sc;
	pos0=pos;
	chisq=chisq_central;
      }
          
      // Abort loop if the chisq is not changing much or increasing too much
      if (isnan(chisq) || (iter2>10 && 
			   (fabs(chisq-chisq_iter)<0.1 
			    || chisq-chisq_iter>10.)))
	break;
      chisq_iter=chisq;

      // Find track parameters where track crosses beam line
      ExtrapolateToVertex(pos0,Sclast,Cclast); 
      Cc=Cclast;
      Sc=Sclast;
      best_pos=pos0;
    }

    if (chisq_iter==1.e16) {
      _DBG_ << "Central fit failed!" <<endl;
      return VALUE_OUT_OF_RANGE;
    }
   
    // Track Parameters at "vertex"
    phi_=Sclast(state_phi,0);
    q_over_pt_=Sclast(state_q_over_pt,0);
    tanl_=Sclast(state_tanl,0);
    x_=best_pos.x();
    y_=best_pos.y();
    z_=best_pos.z();

    if (DEBUG_LEVEL>0)
      cout
      << "Vertex:  p " 
      <<   1./Sclast(state_q_over_pt,0)/cos(atan(Sclast(state_tanl,0)))
      << " theta "  << 90.-180./M_PI*atan(Sclast(state_tanl,0)) 
      << " vertex " << x_<< " " << y_<< " " << z_<<endl;

    // Covariance matrix at vertex
    vector<double>dummy;
    for (unsigned int i=0;i<5;i++){
      dummy.clear();
      for(unsigned int j=0;j<5;j++){
	dummy.push_back(Cclast(i,j));
      }
      cov.push_back(dummy);
    }
    
    // total chisq and ndf
    chisq_=chisq_iter;
    ndf=cdchits.size()-5;
  }

  
  if (DEBUG_HISTS){ 
    TH2F *cdc_residuals=(TH2F*)gROOT->FindObject("cdc_residuals");
    if (cdc_residuals){
      for (unsigned int i=0;i<cdchits.size();i++)
	   cdc_residuals->Fill(cdchits[i]->ring,cdchits[i]->residual);
    }
    /*
    TH2F *cdc_pulls=(TH2F*)gROOT->FindObject("cdc_pulls");
    if (cdc_pulls){
	 //cdc_pulls->Fill(cdchits[cdc_index]->ring,dm/sqrt(var));
    }
    */
  }

  return NOERROR;
}

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b) ((b)>=0.0?fabs(a):-fabs(a))
// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
double DKalmanFilter::BrentsAlgorithm(double ds1,double ds2,
				      double dedx,DVector3 &pos,DVector3 origin,
				      DVector3 dir,  
				      DMatrix &Sc){
  int iter;
  double a,b,d=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-ds1;
  double cx=-ds1-ds2;
  
  a=(ax<cx?ax:cx);
  b=(ax>cx?ax:cx);
  x=w=v=bx;

  // Save the starting position 
  // DVector3 pos0=pos;
  // DMatrix S0(Sc);
  
  // Step to intermediate point
  FixedStep(pos,x,Sc,dedx);
  DVector3 wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
  double u_old=x;

  // initialization
  fw=fv=fx=(pos-wirepos).Mag();

  // main loop
  for (iter=1;iter<=ITMAX;iter++){
    xm=0.5*(a+b);
    tol2=2.0*(tol1=EPS*fabs(x)+ZEPS);
    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      if (pos.z()>=endplate_z-endplate_dz){
	double my_endz=endplate_z-endplate_dz;
	// Check if the minimum doca would occur outside the straw and if so, bring the state 
	// vector back to the end of the straw
	int iter=0;
	while (fabs(pos.z()-my_endz)>EPS && iter<20){
	  u=x-(my_endz-pos.z())*sin(atan(Sc(state_tanl,0)));
	  x=u;
	  // Function evaluation
	  FixedStep(pos,u_old-u,Sc,dedx);
	  u_old=u;
	  iter++;
	}
	//printf("new z %f ds %f \n",pos.z(),x);
      }
      if (pos.z()<=cdc_origin[2]){
	int iter=0;
	while (fabs(pos.z()-cdc_origin[2])>EPS && iter<20){
	  u=x-(cdc_origin[2]-pos.z())*sin(atan(Sc(state_tanl,0)));
	  x=u;
	  // Function evaluation
	  FixedStep(pos,u_old-u,Sc,dedx);
	  u_old=u;
	  iter++;
	}
	//printf("new z %f ds %f \n",pos.z(),x);
      }	
     
      return cx-x;
    }
    // trial parabolic fit
    if (fabs(e)>tol1){
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(q);
      etemp=e;
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
    FixedStep(pos,u_old-u,Sc,dedx);
    u_old=u;
    
    wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
    fu=(pos-wirepos).Mag();

    //printf("Brent: z %f d %f\n",pos.z(),fu);
    
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
  
  return cx-x;
}

// Routine for finding the minimum of a function bracketed between two values
// (see Numerical Recipes in C, pp. 404-405).
double DKalmanFilter::BrentsAlgorithm(double z,double dz,
				      double dedx,DVector3 origin,
				      DVector3 dir,DMatrix S){
  int iter;
  double a,b,d=0.,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0; // will be distance moved on step before last 
  double ax=0.;
  double bx=-dz;
  double cx=-2.*dz;
  
  a=(ax<cx?ax:cx);
  b=(ax>cx?ax:cx);
  x=w=v=bx;
  
  // Save the state vector after the last step
  DMatrix S0(5,1);
  S0=S;

  // Step to intermediate point
  Step(z,z+x,dedx,S0);
  DVector3 wirepos=origin+((z+x-origin.z())/dir.z())*dir;
  DVector3 pos(S0(state_x,0),S0(state_y,0),z+x);

  // initialization
  fw=fv=fx=(pos-wirepos).Mag();

  // main loop
  for (iter=1;iter<=ITMAX;iter++){
    xm=0.5*(a+b);
    tol2=2.0*(tol1=EPS*fabs(x)+ZEPS);
    if (fabs(x-xm)<=(tol2-0.5*(b-a))){
      return x;
    }
    // trial parabolic fit
    if (fabs(e)>tol1){
      r=(x-w)*(fx-fv);
      q=(x-v)*(fx-fw);
      p=(x-v)*q-(x-w)*r;
      q=2.0*(q-r);
      if (q>0.0) p=-p;
      q=fabs(q);
      etemp=e;
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
    S0=S;
    Step(z,z+u,dedx,S0);
    wirepos=origin+((z+u-origin.z())/dir.z())*dir;
    pos.SetXYZ(S0(state_x,0),S0(state_y,0),z+u);
    fu=(pos-wirepos).Mag();

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
  return x;
}


// Routine for finding the minimum of a function bracketed between two values
double DKalmanFilter::GoldenSection(double ds1,double ds2,
				    double dedx,DVector3 pos,DVector3 origin,
				    DVector3 dir,  
				    DMatrix Sc){
  
  double ds=ds1+ds2;
  double s0=0.;
  double s3=-ds; // need to backtrack
  double golden_ratio1=0.61803399;
  double golden_ratio2=1.-golden_ratio1;
  unsigned int iter=0;
  double d1,d2;
  double s1,s2;

  if (fabs(s3+ds1)<fabs(s0+ds1)){
    s1=-ds1;
    s2=-ds1+golden_ratio2*(s3+ds1);
  }
  else{
    s2=-ds1;
    s1=-ds1+golden_ratio2*(s0+ds1);
  }

  DMatrix S0(5,1);

  // Save the state vector and position after the last step
  DVector3 oldpos=pos;
  S0=Sc;

  // Step to one of the intermediate points
  FixedStep(pos,s1,S0,dedx);
  DVector3 wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
  d1=(pos-wirepos).Mag();

  // Step to the other intermediate point
  S0=Sc; 
  pos=oldpos;
  FixedStep(pos,s2,S0,dedx);
  wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
  d2=(pos-wirepos).Mag(); 

  // loop to find the minimum 
  while (fabs(s3-s0)>EPS*(fabs(s1)+fabs(s2))&& iter<MAX_ITER){
    iter++;
    S0=Sc;
    pos=oldpos;
    double stemp;
    if(d2<d1){
      stemp=golden_ratio1*s1+golden_ratio2*s3;
      s0=s1; 
      s1=s2; 
      s2=stemp;
      d1=d2;
      FixedStep(pos,s2,S0,dedx);
      wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
      d2=(pos-wirepos).Mag();
    }
    else{
      stemp=golden_ratio1*s2+golden_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      FixedStep(pos,s1,S0,dedx);
      wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
      d1=(pos-wirepos).Mag();
    }
  }

  if (d1<d2) {
    return s1;
  }
  else{
    return s2;
  }
}

// Routine for finding the minimum of a function bracketed between two values
jerror_t DKalmanFilter::GoldenSection(double &z,double dz,double dEdx,
				      DVector3 origin, DVector3 dir, 
				      DMatrix &S){
  unsigned int iter=0;
  double s0=0.;
  double s3=-2.*dz;
  double golden_ratio1=0.61803399;
  double golden_ratio2=1.-golden_ratio1;
  double s1=-dz;
  double s2=-dz+golden_ratio2*(s3+dz);
  DMatrix S0(5,1);
  S0=S;

  // Step to first intermediate point
  Step(z,z+s1,dEdx,S0);
  double x0=origin.x()+(z+s1-origin.z())*dir.x();
  double y0=origin.y()+(z+s1-origin.z())*dir.y();
  double d1=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
    +(S0(state_y,0)-y0)*(S0(state_y,0)-y0);

  // Step to second intermediate point
  S0=S;
  Step(z,z+s2,dEdx,S0);
  x0=origin.x()+(z+s2-origin.z())*dir.x();
  y0=origin.y()+(z+s2-origin.z())*dir.y();
  double d2=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
    +(S0(state_y,0)-y0)*(S0(state_y,0)-y0);

  // main loop
  while (fabs(s3-s0)>EPS*(fabs(s1)+fabs(s2))&& iter<MAX_ITER
	 && fabs(d1-d2)>EPS2){
    iter++;
    S0=S;
    double stemp;
    if(d2<d1){
      stemp=golden_ratio1*s1+golden_ratio2*s3;
      s0=s1; 
      s1=s2; 
      s2=stemp;
      d1=d2;
      Step(z,z+s2,dEdx,S0);
      x0=origin.x()+(z+s2-origin.z())*dir.x();
      y0=origin.y()+(z+s2-origin.z())*dir.y();
      d2=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
	+(S0(state_y,0)-y0)*(S0(state_y,0)-y0);
    }
    else{
      stemp=golden_ratio1*s2+golden_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      Step(z,z+s1,dEdx,S0);
      x0=origin.x()+(z+s1-origin.z())*dir.x();
      y0=origin.y()+(z+s1-origin.z())*dir.y();
      d1=(S0(state_x,0)-x0)*(S0(state_x,0)-x0)
	+(S0(state_y,0)-y0)*(S0(state_y,0)-y0);    
    }
  }
  if (d1<d2) {
    Step(z,z+s1,dEdx,S);
    z=z+s1;
  }
  else{
    Step(z,z+s2,dEdx,S);
    z=z+s2;
  }
  return NOERROR;
}




// Routine for finding the minimum of a function bracketed between two values
jerror_t DKalmanFilter::GoldenSection(double &ds,double doca,double dedx,
				    DVector3 &pos,
				    DVector3 origin,DVector3 dir,  
				    DMatrix &Sc,DMatrix &Jc){
  double s0=0.;
  double s3=-2.*ds; // need to backtrack
  double magic_ratio1=0.61803399;
  double magic_ratio2=1.-magic_ratio1;
  double s2=-ds;
  double s1=-ds*(1.-magic_ratio2);
  double d2=doca;
  DMatrix S0(5,1);
  
  // Save the state vector and position after the last step
  DVector3 oldpos=pos;
  S0=Sc;

  // Step to one of the intermediate points
  FixedStep(pos,s1,S0,dedx);
  double d1=S0(state_D,0);
  while (fabs(s3-s0)>EPS*(fabs(s1)+fabs(s2))){
    // Check for exit condition near zero
    if (fabs(s3)<EPS && fabs(s0)<EPS && fabs(s1)<EPS && fabs(s2)<EPS) break;
    S0=Sc;
    pos=oldpos;
    double stemp;
    if(fabs(d2)<fabs(d1)){
      stemp=magic_ratio1*s1+magic_ratio2*s3;
      s0=s1; 
      s1=s2; 
      s2=stemp;
      d1=d2;
      FixedStep(pos,s2,S0,dedx);
      d2=S0(state_D,0);
    }
    else{
      stemp=magic_ratio1*s2+magic_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      FixedStep(pos,s1,S0,dedx);
      d1=S0(state_D,0);
    }
  }
  pos=oldpos;
  if (fabs(d1)<fabs(d2)) {
    StepJacobian(pos,origin,dir,s1,Sc,dedx,Jc);
    ds=s1;
  }  
  else{
    StepJacobian(pos,origin,dir,s2,Sc,dedx,Jc);
    ds=s2;
  }
  return NOERROR;
}



// Kalman engine for Central tracks; updates the position on the trajectory
// after the last hit (closest to the target) is added
jerror_t DKalmanFilter::KalmanCentral(double mass_hyp,double anneal_factor,
				      DMatrix &Sc,DMatrix &Cc,
				      DVector3 &pos,double &chisq){
  DMatrix H(1,5);  // Track projection matrix
  DMatrix H_T(5,1); // Transpose of track projection matrix
  DMatrix J(5,5);  // State vector Jacobian matrix
  DMatrix JT(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,1);  // Kalman gain matrix
  double V=1.6*1.6/12.;  // Measurement variance
  double InvV; // inverse of variance
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix S0(5,1),S0_(5,1); // state vector

  // Initialize the chi2 for this part of the track
  chisq=0.;

  // path length
  double s=central_traj[0].s;
  // path length increment
  double ds=-CDC_STEP_SIZE,ds2=0.;

  // beginning position
  pos(0)=central_traj[0].pos.x()-Sc(state_D,0)*sin(Sc(state_phi,0));
  pos(1)=central_traj[0].pos.y()+Sc(state_D,0)*cos(Sc(state_phi,0));
  pos(2)=Sc(state_z,0);

  // Wire origin and direction
  unsigned int cdc_index=cdchits.size()-1;
  DVector3 origin=cdchits[cdc_index]->origin;
  DVector3 dir=cdchits[cdc_index]->dir;
  DVector3 wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;

  // doca variables
  double doca,old_doca=(pos-wirepos).Mag();

  // energy loss
  double dedx=0.;

  // Boolean for flagging when we are done with measurements
  bool more_measurements=true;

  // Initialize S0_ and perform the loop over the trajectory
  S0_=DMatrix(*central_traj[0].S);

  for (unsigned int k=1;k<central_traj.size();k++){
    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=DMatrix(*central_traj[k].S);
    J=DMatrix(*central_traj[k].J);
    Q=DMatrix(*central_traj[k].Q);

    // State S is perturbation about a seed S0
    dS=Sc-S0_;
    //dS.Zero();

    // Update the actual state vector and covariance matrix
    Sc=S0+J*dS;  
    JT=DMatrix(DMatrix::kTransposed,J);
    Cc=J*(Cc*JT)+Q;   

    // Copy covariance to the trajectory vector
    for (unsigned int i=0;i<5;i++){
      for (unsigned int j=0;j<5;j++){
	central_traj[k].C->operator()(i,j)=Cc(i,j);
      }
    }

    // update position based on new doca to reference trajectory
    pos(0)=central_traj[k].pos.x()-Sc(state_D,0)*sin(Sc(state_phi,0));
    pos(1)=central_traj[k].pos.y()+Sc(state_D,0)*cos(Sc(state_phi,0));
    pos(2)=Sc(state_z,0);

    // Save the current state of the reference trajectory
    S0_=S0;

    // new wire position
    wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;

    // new doca
    doca=(pos-wirepos).Mag();

    // Check if the doca is no longer decreasing
    if ((doca>old_doca) && (pos.z()<endplate_z && pos.z()>cdc_origin[2])
	&& more_measurements){
      //printf("----------------------------------------\n");
      //if (pos.z()<endplate_z && pos.z()>cdc_origin[2])
      {
	// Save values at end of current step
	DVector3 pos0=central_traj[k].pos;
	
	// dEdx for current position along trajectory
	double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	if (do_energy_loss){
	  dedx=GetdEdx(MASS,q_over_p,central_traj[k].Z,
		       central_traj[k].A,central_traj[k].density);
	}
	
	// Variables for the computation of D at the doca to the wire
	double D=Sc(state_D,0);
	double q=(Sc(state_q_over_pt,0)>0)?1.:-1.;
	double qpt=1./Sc(state_q_over_pt,0);
	double sinphi=sin(Sc(state_phi,0));
	double cosphi=cos(Sc(state_phi,0));
	double Bx,By,Bz;
	bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
	double Bz_=fabs(Bz);
	
	// We've passed the true minimum; now use Brent's algorithm
	// to find the best doca.  See Numerical Recipes in C, pp 404-405
	ds2=BrentsAlgorithm(ds,ds,dedx,pos,origin,dir,Sc);
	
	int numstep=(int)(ds2/CDC_STEP_SIZE);
	double myds=CDC_STEP_SIZE;
	if (ds2<0) myds*=-1.;
	double ds3=ds2-CDC_STEP_SIZE*numstep;
	// propagate covariance matrix along the reference trajectory.
	// We ignore the tiny amount of multiple scattering for these small 
	// steps.
	for (int j=0;j<abs(numstep);j++){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,myds,S0,dedx,J);

	  // Update covariance matrix
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT);
	  
	  // Step along reference trajectory 
	  FixedStep(pos0,myds,S0,dedx);	
	}
	if (fabs(ds3)>EPS2){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,ds3,S0,dedx,J);
	  
	  // Step along reference trajectory 
	  FixedStep(pos0,ds3,S0,dedx);
	  
	  // Update covariance matrix
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT);
	}

	// Compute the value of D (signed distance to the reference trajectory)
	// at the doca to the wire
	DVector3 dpos1=pos0-central_traj[k].pos;
	double rc=sqrt(dpos1.Perp2()
		       +2.*(D+qpt/qBr2p/Bz_)*(dpos1.x()*sinphi-dpos1.y()*cosphi)
		       +(D+qpt/qBr2p/Bz_)*(D+qpt/qBr2p/Bz_));
	Sc(state_D,0)=q*rc-qpt/qBr2p/Bz_;
		
	// wire position
	wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
	
	//doca 
	doca=(pos-wirepos).Perp();
      
	// Measurement
	double lambda=atan(Sc(state_tanl,0));
	q_over_p=Sc(state_q_over_pt,0)*cos(lambda);
	double sinl=sin(lambda);
	double cosl=cos(lambda);
	double one_over_beta=sqrt(1.+mass_hyp*mass_hyp*q_over_p*q_over_p);
	double s=central_traj[k].s+ds2;
	double measurement=0.;
	if (pass==kTimeBased){
	  measurement=CDC_DRIFT_SPEED*(cdchits[cdc_index]->t
				       -s*one_over_beta/SPEED_OF_LIGHT);

	  //printf("t %f \n",s*one_over_beta/SPEED_OF_LIGHT);
	  
	  // Measurement error
	  V=anneal_factor*0.000225;
	  /*
	    +DRIFT_SPEED*DRIFT_SPEED
	    *(VAR_S*one_over_beta*one_over_beta/SPEED_OF_LIGHT/SPEED_OF_LIGHT 
	      +mass_hyp*mass_hyp*mass_hyp*mass_hyp*s*s
	      /one_over_beta/one_over_beta/SPEED_OF_LIGHT/SPEED_OF_LIGHT
	      *q_over_p*q_over_p*cosl*cosl
	      *(Cc(state_q_over_pt,state_q_over_pt)
		+q_over_p*q_over_p*sinl*sinl*Cc(state_tanl,state_tanl)
		-2.*q_over_p*sinl*Cc(state_tanl,state_q_over_pt)
		));
	  */
	}
	
	// prediction for measurement  
	DVector3 diff=pos-wirepos;
	DVector3 mdir=pos-wirepos-(diff.Dot(dir))*dir;
	mdir.SetMag(1.0);
	double prediction=diff.Dot(mdir);
	
	// Projection matrix        
	sinphi=sin(Sc(state_phi,0));
	cosphi=cos(Sc(state_phi,0));
	double ux=dir.x();
	double uy=dir.y();
	double uxuy=ux*uy;
	double ux2=ux*ux;
	double uy2=uy*uy;
	double dx=diff.x();
	double dy=diff.y();
	if (prediction>0.){
	  H(0,state_D)=H_T(state_D,0)
	    =(dy*(uxuy*sinphi+(1.-uy2)*cosphi)-dx*((1.-ux2)*sinphi+uxuy*cosphi))/prediction;
	  H(0,state_phi)=H_T(state_phi,0)
	    =-Sc(state_D,0)*(dx*((1.-ux2)*cosphi-uxuy*sinphi)+dy*((1.-uy2)*sinphi-uxuy*cosphi))/prediction;
	  H(0,state_z)=H_T(state_z,0)
	    =-dir.z()*(dx*ux+dy*uy)/prediction;
	}
	
	// Difference and variance
	double var=V,var_pred=0.;
	if (prediction>0.){
	  var_pred=(H*(Cc*H_T))(0,0);
	  var+=var_pred;
	}
	double dm=measurement-prediction;

	if (var_pred<0.){
	  /*
	  Cc.Print();
	  cout << "Negative variance???" << var_pred << endl;
	  H.Print();
	  */
	  return VALUE_OUT_OF_RANGE;
	}
	
	// probability
	double p=exp(-0.5*dm*dm/var)/sqrt(2.*M_PI*var);
	p=1.;
	//if (fabs(dm)/sqrt(var)>3.) dm=3.*sqrt(var)*(dm>0?1.:-1.);

	if (DEBUG_LEVEL>0) 
	  cout 
	    << "ring " << cdchits[cdc_index]->ring << 
	    " Dm " << measurement << 
	    " Dm-Dpred " << dm << " sigma " 
	    << sqrt(var) << " " << sqrt(V)<< " p " 
	    << 1./(Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)))) 
	    << " theta " << 90.-180./M_PI*atan(Sc(state_tanl,0)) 
	    << " t " << s*one_over_beta/SPEED_OF_LIGHT
	    << " x " << pos.x() << " y " << pos.y() << " z " << pos.z()
	    << endl;
	
	// Inverse of variance
	InvV=1./(p*V+var_pred);
	
	// Compute Kalman gain matrix
	K=InvV*(Cc*H_T);
	
	// Update the state vector 
	dS=p*dm*K;
	//dS.Zero();
	Sc=Sc+dS;
	
	// Path length in active volume
	//path_length+=?

	// Update state vector covariance matrix
	Cc=Cc-(K*(H*Cc));  
	
	/* This is more accurate, but is it worth it??
	// update position on current trajectory based on corrected doca to 
	// reference trajectory
	pos=pos0;
	pos(0)+=-Sc(state_D,0)*sin(Sc(state_phi,0));
	pos(1)+= Sc(state_D,0)*cos(Sc(state_phi,0));
	pos(2)=Sc(state_z,0);   
      
	// wire position
	wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;

	// revised prediction
	diff=pos-wirepos;
	mdir=pos-wirepos-(diff.Dot(dir))*dir;
	mdir.SetMag(1.0);
	prediction=diff.Dot(mdir);
	*/

	// calculate the residual
	dm*=(1.-(H*K)(0,0));
	//dm=measurement-prediction;
	cdchits[cdc_index]->residual=dm;
	
	// Update chi2 for this hit
	var=V*(1.-(H*K)(0,0));
	chisq+=dm*dm/var;      

	// propagate the covariance matrix to the next point on the trajectory
	for (int j=0;j<abs(numstep);j++){
	  DMatrix Stemp(S0);
	
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,-myds,S0,dedx,J);
	  
	  // Update covariance matrix, ignoring multiple scattering
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT);

	  // Step along reference trajectory 
	  FixedStep(pos0,-myds,S0,dedx);	
	  
	  // Step to the next point on the trajectory
	  Sc=S0+J*(Sc-Stemp); 
	}
	if (fabs(ds3)>EPS){
	  // Compute the Jacobian matrix
	  StepJacobian(pos0,origin,dir,-ds3,S0,dedx,J);
	  
	  // Update covariance matrix
	  JT=DMatrix(DMatrix::kTransposed,J);
	  Cc=J*(Cc*JT);
	}

	for (unsigned int i=0;i<5;i++){
	  for (unsigned int j=0;j<5;j++){
	    central_traj[k].C->operator()(i,j)=Cc(i,j);
	  }
	}
	
	// Step to the next point on the trajectory
	Sc=S0_+J*(Sc-S0); 
	
	// update position on current trajectory based on corrected doca to 
	// reference trajectory
	pos=central_traj[k].pos;
	pos(0)+=-Sc(state_D,0)*sin(Sc(state_phi,0));
	pos(1)+= Sc(state_D,0)*cos(Sc(state_phi,0));
	pos(2)=Sc(state_z,0);   
      }
      
      // new wire origin and direction
      if (cdc_index>0){
	cdc_index--;
	origin=cdchits[cdc_index]->origin;
	dir=cdchits[cdc_index]->dir;
      }
      else{
	origin.SetXYZ(0.,0.,65.);
	dir.SetXYZ(0,0,1.);
	more_measurements=false;
      }
      
      // Update the wire position
      wirepos=origin+((pos.z()-origin.z())/dir.z())*dir;
      
      //s+=ds2;
      // new doca
      doca=(pos-wirepos).Mag();
    }

    old_doca=doca;
  } 

  if (DEBUG_LEVEL>0)
    cout 
      << " p " 
      << 1./(Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)))) 
      << " theta " << 90.-180./M_PI*atan(Sc(state_tanl,0)) 
      << " vertex " << pos.x() << " " << pos.y() <<" " << pos.z() <<endl;
  
  // update internal variables
  phi_=Sc(state_phi,0);
  q_over_pt_=Sc(state_q_over_pt,0);
  tanl_=Sc(state_tanl,0);

  x_=pos.x();
  y_=pos.y();
  z_=Sc(state_z,0);
  
  chisq*=anneal_factor;
  
  return NOERROR;
}


// Kalman engine for forward tracks
jerror_t DKalmanFilter::KalmanForward(double anneal_factor, DMatrix &S, 
				      DMatrix &C,
				      double &chisq){
  DMatrix M(2,1);  // measurement vector
  DMatrix Mpred(2,1); // prediction 
  DMatrix H(2,5);  // Track projection matrix
  DMatrix H_T(5,2); // Transpose of track projection matrix
  DMatrix J(5,5);  // State vector Jacobian matrix
  DMatrix J_T(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,2);  // Kalman gain matrix
  DMatrix V(2,2);  // Measurement covariance matrix
  DMatrix V1(2,2);  // same, with the contribution from C
  DMatrix R(2,1);  // Filtered residual
  DMatrix R_T(1,2);  // ...and its transpose
  DMatrix RC(2,2);  // Covariance of filtered residual
  DMatrix InvRC(2,2); // and its inverse
  DMatrix S0(5,1),S0_(5,1); //State vector
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix InvV(2,2); // Inverse of error matrix

  // Path length in active volume
  path_length=0;

  // overall path length
  double s=forward_traj[0].s;

  // Initialize chi squared
  chisq=0;

  // Initialize error matrix 
  V(0,0)=1.0*1.0/12;
  V(1,1)=0.32*0.32/12.;

  S0_=DMatrix(*forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size()-1;k++){
    //path length
    s-=STEP_SIZE*sqrt(1.+S(state_tx,0)*S(state_tx,0)
		      +S(state_ty,0)*S(state_ty,0));

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=DMatrix(*forward_traj[k].S);
    J=DMatrix(*forward_traj[k].J);
    Q=DMatrix(*forward_traj[k].Q);

    // State S is perturbation about a seed S0
    dS=S-S0_;

    // Update the actual state vector and covariance matrix
    S=S0+J*dS;
    C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;   

    // Save the current state of the reference trajectory
    S0_=S0;

    // Add the hit
    if (forward_traj[k].h_id>0){
      unsigned int id=forward_traj[k].h_id-1;
      
      double cosa=fdchits[id]->cosa;
      double sina=fdchits[id]->sina;
      double u=fdchits[id]->uwire;
      double v=fdchits[id]->vstrip;
      double x=S(state_x,0);
      double y=S(state_y,0);
      double tx=S(state_tx,0);
      double ty=S(state_ty,0);
      double du=x*cosa-y*sina-u;
      double tu=tx*cosa-ty*sina;
      double alpha=atan(tu);
      double cosalpha=cos(alpha);
      double sinalpha=sin(alpha);

      // The next measurement 
      M(0,0)=0.;
      M(1,0)=v;

      if (pass==kTimeBased){
	// Compute drift distance
	double one_over_beta
	  =sqrt(1.+S(state_q_over_p,0)*S(state_q_over_p,0)*MASS*MASS);
	double tflight=forward_traj[k].s*one_over_beta/SPEED_OF_LIGHT;
	double drift=DRIFT_SPEED*(fdchits[id]->t-tflight);  
	drift*=(du>0?1.:-1.);
	
	// Angles of incidence to the measurement plane
	double phi=atan2(S(state_y,0),S(state_x,0));
	double cosphi=cos(phi);
	double sinphi=sin(phi);

	// Drift distance
	M(0,0)=drift;
	
	// Correction for lorentz effect
	double nz=fdchits[id]->nz;
	double nr=fdchits[id]->nr;
	double dv=nz*drift*sinalpha*cosphi-nr*drift*cosalpha;
	M(1,0)=v+dv;// with correction for Lorentz effect
	
	// ... and its covariance matrix 
	V(0,0)=anneal_factor*0.000225;
	V(1,1)=fdc_y_variance(alpha,drift);

	// Variance due to Lorentz correction
	double var_alpha
	  =(C(state_tx,state_tx)*cosa*cosa+C(state_ty,state_ty)*sina*sina
	    -2.*sina*cosa*C(state_tx,state_ty))/(1.+tu*tu)/(1.+tu*tu);
	double var_phi=(y*y*C(state_x,state_x)+x*x*C(state_y,state_y)
			-2*x*y*C(state_x,state_y))/(x*x+y*y)/(x*x+y*y);
	
	V(1,1)+=dv*dv*V(0,0)/drift/drift
	  +drift*drift*(nz*cosalpha*cosphi+nr*sinalpha)
	  *(nz*cosalpha*cosphi+nr*sinalpha)*var_alpha
	  +drift*drift*nz*nz*sinalpha*sinalpha*sinphi*sinphi*var_phi;
	
	V(1,1)*=anneal_factor;
      }

      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      H(0,state_x)=H_T(state_x,0)=cosa*cosalpha;
      H(1,state_x)=H_T(state_x,1)=sina;
      H(0,state_y)=H_T(state_y,0)=-sina*cosalpha;
      H(1,state_y)=H_T(state_y,1)=cosa;
      double temp=tx*cosa-ty*sina;
      H(0,state_ty)=H_T(state_ty,0)
	=du*sina*temp/sqrt(1.+temp*temp)/(1.+temp*temp);
      H(0,state_tx)=H_T(state_tx,0)
	=-du*cosa*temp/sqrt(1.+temp*temp)/(1.+temp*temp);

      // Updated error matrix
      V1=V+H*(C*H_T);

      // Calculate the inverse of V
      double det=V1(0,0)*V1(1,1)-V1(0,1)*V1(1,0);
      if (det!=0){
	InvV(0,0)=V1(1,1)/det;
	InvV(1,0)=-V1(1,0)/det;
	InvV(0,1)=-V1(0,1)/det;
	InvV(1,1)=V1(0,0)/det;
      }
      else{
	cout 
	  << "Kalman filter:  Singular matrix..." << endl;
	return UNRECOVERABLE_ERROR;
      }
      
      // Compute Kalman gain matrix
      K=C*(H_T*InvV);
     
      // Update the state vector 
      Mpred(0,0)=du*cos(alpha);
      Mpred(1,0)= y*cosa+x*sina;
      S=S+K*(M-Mpred); 


      //.      printf("z %f Diff\n",forward_traj[k].pos.z());
      //(M-Mpred).Print();
      // Path length in active volume
      path_length+=STEP_SIZE*sqrt(1.+S(state_tx,0)*S(state_tx,0)
			    +S(state_ty,0)*S(state_ty,0));
      
      // Update state vector covariance matrix
      C=C-K*(H*C);    
      
      // Residuals
      x=S(state_x,0);
      y=S(state_y,0);
      tx=S(state_tx,0);
      ty=S(state_ty,0);
      du=x*cosa-y*sina-u;
      R(0,0)=M(0,0)-du*cos(atan(tx*cosa-ty*sina));
      R(1,0)=M(1,0)-(y*cosa+x*sina);
      R_T=DMatrix(DMatrix::kTransposed,R);
      RC=V-H*(C*H_T);

      fdchits[id]->xres=R(0,0);
      fdchits[id]->yres=R(1,0);

      // Calculate the inverse of RC
      det=RC(0,0)*RC(1,1)-RC(0,1)*RC(1,0);
      if (det!=0){
	InvRC(0,0)=RC(1,1)/det;
	InvRC(1,0)=-RC(1,0)/det;
	InvRC(0,1)=-RC(0,1)/det;
	InvRC(1,1)=RC(0,0)/det;
      }
      else{ 
	cout 
	  << "Kalman filter:  Singular matrix RC..." << endl;
	return UNRECOVERABLE_ERROR;
      }
      
      // Update chi2 for this segment
      chisq+=(R_T*(InvRC*R))(0,0);
    }

  }
  chisq*=anneal_factor;

  // Final position for this leg
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=forward_traj[forward_traj.size()-1].pos.Z();

  if (DEBUG_LEVEL>0)
    cout << "Position after forward filter: " << x_ << ", " << y_ << ", " << z_ <<endl;

  return NOERROR;
}
// Kalman engine for forward tracks -- this routine adds CDC hits
jerror_t DKalmanFilter::KalmanForwardCDC(double mass_hyp,double anneal,
					 DMatrix &S, 
					 DMatrix &C,double &chisq){
  DMatrix H(1,5);  // Track projection matrix
  DMatrix H_T(5,1); // Transpose of track projection matrix
  DMatrix J(5,5);  // State vector Jacobian matrix
  DMatrix J_T(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,1);  // Kalman gain matrix
  DMatrix S0(5,1),S0_(5,1); //State vector
  DMatrix dS(5,1);  // perturbation in state vector
  double V=0.000225*anneal; // measurement variance
  double InvV;  // inverse of variance

  // Path length in active volume
  path_length=0;
  double s=forward_traj_cdc[0].s,ds=0.;
  // position variables
  double x=S(state_x,0);
  double y=S(state_y,0);
  double z=forward_traj_cdc[0].pos.z();

  // wire information  
  unsigned int cdc_index=cdchits.size()-1;
  DVector3 origin=cdchits[cdc_index]->origin;
  DVector3 dir=cdchits[cdc_index]->dir;
  DVector3 wirepos=origin+((z-origin.z())/dir.z())*dir;
  bool more_measurements=true;

  // doca variables
  double doca,old_doca=sqrt((x-wirepos.x())*(x-wirepos.x())
			    +(y-wirepos.y())*(y-wirepos.y()));
  /*
  printf("p %f theta %f\n",1./S(state_q_over_p,0),
	 90-180/M_PI*atan(1./sqrt(S(state_tx,0)*S(state_tx,0)
				  +S(state_ty,0)*S(state_ty,0))));
  */
  //C.Print();
  
  // loop over entries in the trajectory
  S0_=DMatrix(*forward_traj_cdc[0].S);
  for (unsigned int k=1;k<forward_traj_cdc.size()-1;k++){
    z=forward_traj_cdc[k].pos.z();

    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=DMatrix(*forward_traj_cdc[k].S);
    J=DMatrix(*forward_traj_cdc[k].J);
    Q=DMatrix(*forward_traj_cdc[k].Q);

    // State S is perturbation about a seed S0
    dS=S-S0_;

    // Update the actual state vector and covariance matrix
    S=S0+J*dS;
    C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;   

    // path length
    s-=STEP_SIZE*sqrt(1.+S(state_tx,0)*S(state_tx,0)
		      +S(state_ty,0)*S(state_ty,0));
   
    // Save the current state of the reference trajectory
    S0_=S0;

    // new wire position
    wirepos=origin+((z-origin.z())/dir.z())*dir;

    // new doca
    x=S(state_x,0);
    y=S(state_y,0);
    doca=sqrt((x-wirepos.x())*(x-wirepos.x())
	      +(y-wirepos.y())*(y-wirepos.y()));

    // Check if the doca is no longer decreasing
    //    printf("x %f xw %f y %f yw %f doca %f old_doca %f\n",x,wirepos.x(),
    //	   y,wirepos.y(),doca,old_doca);
    if ((doca>old_doca) && z<endplate_z && more_measurements){
      // Get energy loss 
      double dedx=0.;
      if (do_energy_loss){
	dedx=GetdEdx(MASS,S(state_q_over_p,0),forward_traj_cdc[k].Z,
		     forward_traj_cdc[k].A,forward_traj_cdc[k].density);
      }

      // We have bracketed the minimum doca 
      double dz=BrentsAlgorithm(z,-CDC_FORWARD_STEP_SIZE,dedx,origin,dir,S);
      double newz=z+dz;
      ds=Step(z,newz,dedx,S);
      s+=ds;

      // Step reference trajectory by dz
      Step(z,newz,dedx,S0);

      // propagate error matrix to z-position of hit
      StepJacobian(z,newz,S0,dedx,J);
      C=J*(C*DMatrix(DMatrix::kTransposed,J));

      // Wire position at current z
      wirepos=origin+((newz-origin.z())/dir.z())*dir;
      double xw=wirepos.x();
      double yw=wirepos.y();

      // predicted doca taking into account the orientation of the wire
      double ux=dir.x();
      double uy=dir.y();
      double uxuy=ux*uy;
      double ux2=ux*ux;
      double uy2=uy*uy;
      double dy=S(state_y,0)-yw;
      double dx=S(state_x,0)-xw;      
      double d=sqrt(dx*dx*(1.-ux2)+dy*dy*(1.-uy2)-2.*dx*dy*uxuy);
     
      // Track projection
      if (d>0.){
	H(0,state_x)=H_T(state_x,0)=(dx*(1.-ux2)-dy*uxuy)/d;
	H(0,state_y)=H_T(state_y,0)=(dy*(1.-uy2)-dx*uxuy)/d;
      }

      //H.Print();
	
      // The next measurement
      double dm=0.;
      double V=1.6*1.6/12.;
      if (pass==kTimeBased){
	double one_over_beta
	  =sqrt(1.+mass_hyp*mass_hyp*S(state_q_over_p,0)*S(state_q_over_p,0));
	//s=forward_traj_cdc[k].s;
	dm=CDC_DRIFT_SPEED*(cdchits[cdc_index]->t-s*one_over_beta
				   /SPEED_OF_LIGHT);
	// variance
	V=0.000225*anneal;
      }
      // variance including prediction
      double var=V,var_pred=0.;
      if (d>0.){
	var_pred=(H*(C*H_T))(0,0);
	var+=var_pred;
      }
      if (var<0.){
	cout << "Negative variance???" << endl;
	return VALUE_OUT_OF_RANGE;
      }
      
      //printf("prediction %f measurement %f V %f %f sig %f\n",d,dm,V,var,sqrt(var));

       // probability
      double p=exp(-0.5*(dm-d)*(dm-d)/var);   
      p=1.;
      //if (fabs(dm-d)/sqrt(var)>3.) dm=d+3.*sqrt(var)*(dm-d>0?1.:-1.);
      
      // Inverse of covariance matrix 
      InvV=p/(p*V+var_pred);
      
      // Compute Kalman gain matrix
      K=InvV*(C*H_T);
      
      //printf("invV %f\n",InvV);
      //C.Print();
      
      //K.Print();
      
      // Update the state vector 
      S=S+(dm-d)*K;
      
      //printf("State\n");
      //S.Print();
      
      // Path length in active volume
      // path_length+=?
      
      //printf("correction to C\n");
      //(K*(H*C)).Print();
      
      // Update state vector covariance matrix
      C=C-K*(H*C);    
      
      // doca after correction
      dy=S(state_y,0)-yw;
      dx=S(state_x,0)-xw;      
      d=sqrt(dx*dx*(1.-ux2)+dy*dy*(1.-uy2)-2.*dx*dy*uxuy);
	
      // Residual
      double res=dm-d;
	
	// Update chi2 for this segment
      chisq+=anneal*res*res/(V-(H*(C*H_T))(0,0));
	/*
	  printf("res %f chisq contrib %f varpred %f\n",res,
	  anneal*res*res/(V-(H*(C*H_T))(0,0)),
	  ( H*(C*H_T))(0,0)
	  );
	*/
	     
      // multiple scattering
      if (do_multiple_scattering){
	GetProcessNoise(MASS,ds,newz,forward_traj_cdc[k].X0,S0,Q);
      }
      // Step C back to the z-position on the reference trajectory
      StepJacobian(newz,z,S0,dedx,J);
      C=J*(C*DMatrix(DMatrix::kTransposed,J))+Q;
      
      // new wire origin and direction
      if (cdc_index>0){
	cdc_index--;
	origin=cdchits[cdc_index]->origin;
	dir=cdchits[cdc_index]->dir;
      }
      else{
	origin.SetXYZ(0.,0.,65.);
	dir.SetXYZ(0,0,1.);
	more_measurements=false;
      }
      
      // Step S to current position on the reference trajectory
      s+=Step(newz,z,dedx,S);

      // Update the wire position
      wirepos=origin+((z-origin.z())/dir.z())*dir;
      
      //s+=ds2;
      // new doca
      x=S(state_x,0);
      y=S(state_y,0);
      doca=sqrt((x-wirepos.x())*(x-wirepos.x())
		+(y-wirepos.y())*(y-wirepos.y()));
    }
    old_doca=doca;
  }

  // Final position for this leg
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=forward_traj_cdc[forward_traj_cdc.size()-1].pos.Z();
  
  if (DEBUG_LEVEL>0)
    cout << "Position after forward cdc filter: " << x_ << ", " << y_ << ", " << z_ <<endl;
  
  return NOERROR;
}

// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.  Converts to the central
// track representation at the end.
jerror_t DKalmanFilter::ExtrapolateToVertex(DMatrix &S,DMatrix &C){
  DMatrix J(5,5);  //.Jacobian matrix
  DMatrix JT(5,5); // and its transpose
  DMatrix Q(5,5); // multiple scattering matrix
  DMatrix Sc(5,1),Cc(5,5);  // central representation

  // position variables
  double z=z_,newz=z_;
  double dz=-STEP_SIZE;
  double ds=sqrt(1.+S(state_tx,0)*S(state_tx,0)+S(state_ty,0)*S(state_ty,0))
    *dz;
  double r2_old=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);

  double dEdx=0.;

  // Check the direction of propagation
  DMatrix S0(5,1);
  S0=S;
  Step(z,z+dz,dEdx,S0);
  double r2=S0(state_x,0)*S0(state_x,0)+S0(state_y,0)*S0(state_y,0);
  if (r2>r2_old) dz*=-1.;

  // material properties
  double Z,A,density,X0;
  DVector3 pos;  // current position along trajectory

  while (z>Z_MIN && sqrt(r2_old)<65. && z<Z_MAX){
    // get material properties from the Root Geometry
    pos.SetXYZ(S(state_x,0),S(state_y,0),z);
    if (RootGeom->FindMat(pos,density,A,Z,X0)!=NOERROR){
      _DBG_ << "Material error in ExtrapolateToVertex! " << endl;
      break;
    }

    // Get dEdx for the upcoming step
    if (do_energy_loss){
      dEdx=GetdEdx(MASS,S(state_q_over_p,0),Z,A,density); 
    }
 
    // Get the contribution to the covariance matrix due to multiple 
    // scattering
    if (do_multiple_scattering)
      GetProcessNoise(MASS,ds,z,X0,S,Q);

    newz=z+dz;
    // Compute the Jacobian matrix
    StepJacobian(z,newz,S,dEdx,J);

    // Propagate the covariance matrix
    JT=DMatrix(DMatrix::kTransposed,J);
    C=J*(C*JT)+Q;

    // Step through field
    ds=Step(z,newz,dEdx,S);
    r2=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);
    if (r2>r2_old){
      DVector3 origin(0,0,65.);
      DVector3 dir(0,0,1.);

      // Find position of best doca to beam line
      //GoldenSection(newz,dz,dEdx,origin,dir,S);
      dz=BrentsAlgorithm(newz,dz,dEdx,origin,dir,S);
      Step(newz,newz+dz,dEdx,S);
      newz+=dz;
      break;
    }
    r2_old=r2;
    z=newz;
  }
  // update internal variables
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=newz;

  return NOERROR;
}


// Propagate track to point of distance of closest approach to origin
jerror_t DKalmanFilter::ExtrapolateToVertex(DVector3 &pos,
					    DMatrix &Sc,DMatrix &Cc){
  DMatrix Jc(5,5);  //.Jacobian matrix
  DMatrix JcT(5,5); // and its transpose
  DMatrix Q(5,5); // multiple scattering matrix

  // Initialize the beam position = center of target, and the direction
  DVector3 origin(0,0,65.);  
  DVector3 dir(0,0,1.);

  // Position and step variables
  double r=pos.Perp();

  // Check if we are outside the nominal beam radius 
  if (r>BEAM_RADIUS){
    double ds=-CDC_STEP_SIZE; // step along path in cm
    double r_old=r;
    Sc(state_D,0)=r;
    
    // Energy loss
    double dedx=0.;
    
    // Check direction of propagation
    DMatrix S0(5,1);
    S0=Sc; 
    DVector3 pos0=pos;
    FixedStep(pos0,ds,S0,dedx);
    r=pos0.Perp();
    if (r>r_old) ds*=-1.;
    
    // Track propagation loop
    while (Sc(state_z,0)>Z_MIN && Sc(state_z,0)<Z_MAX  
	   && r<R_MAX){ 
      // get material properties from the Root Geometry
      double density,A,Z,X0;
      if (RootGeom->FindMat(pos,density,A,Z,X0)!=NOERROR){
	_DBG_ << "Material error in ExtrapolateToVertex! " << endl;
	break;
      }
      
      // Get dEdx for the upcoming step
      double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
      if (do_energy_loss){
	dedx=GetdEdx(MASS,q_over_p,Z,A,density); 
      }
      
      // Compute the Jacobian matrix
      StepJacobian(pos,origin,dir,ds,Sc,dedx,Jc);
      
      // Multiple scattering
      if (do_multiple_scattering)
	GetProcessNoiseCentral(MASS,ds,pos,X0,Sc,Q);
      // Propagate the covariance matrix
      JcT=DMatrix(DMatrix::kTransposed,Jc);
      Cc=Jc*(Cc*JcT)+Q;
      
      // Propagate the state through the field
      FixedStep(pos,ds,Sc,dedx);
      
      r=pos.Perp();
      if (r>r_old) {
	// We've passed the true minimum; now use the Brent's algorithm
	// to find the best doca.  See Numerical Recipes in C, pp 404-405
	ds=BrentsAlgorithm(ds,ds,dedx,pos,origin,dir,Sc);
		
	break;
      }
      r_old=r;
    }   
  } // if (r>BEAM_RADIUS)
  
  return NOERROR;
}

