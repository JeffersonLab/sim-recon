//************************************************************************
// DKalmanFilter.cc
//************************************************************************

#include "DKalmanFilter.h"
#include "CDC/DCDCTrackHit.h"
#include "DANA/DApplication.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DGeometry.h"
#include <TDecompLU.h>

#include <TH2F.h>
#include <TROOT.h>

#include <math.h>

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define EPS 3.0e-8
#define EPS2 1.e-3
#define DEDX_ENDPLATE (3.95e-3) // Carbon
#define DEDX_AIR  (2.19e-6) // GeV/cm
#define DEDX_LH2 ( 2.856e-4) 
#define DEDX_ROHACELL (0.4e-3) // low density carbon for now!!
#define BEAM_RADIUS  0.1 
#define MAX_ITER 25
#define STEP_SIZE 0.5
#define NUM_ITER 5
#define Z_MIN 45.
#define Z_MAX 85.
#define SPEED_OF_LIGHT 29.98
#define CDC_DRIFT_SPEED 55e-4

// Local boolean routines for sorting
bool DKalmanHit_cmp(DKalmanHit_t *a, DKalmanHit_t *b){
  return a->z<b->z;
}
bool DKalmanFDCHit_cmp(DKalmanFDCHit_t *a, DKalmanFDCHit_t *b){
  return a->z<b->z;
}
bool DKalmanCDCHit_cmp(DKalmanCDCHit_t *a, DKalmanCDCHit_t *b){
  double ra=a->origin.Perp();
  double rb=b->origin.Perp();
  return (rb>ra);
}

DKalmanFilter::DKalmanFilter(const DMagneticFieldMap *bfield,
			     const DGeometry *dgeom){
  this->bfield=bfield;
  this->geom=dgeom;
  hits.clear();

  // Get the position of the exit of the CDC endplate from DGeometry
  geom->GetCDCEndplate(endplate_z,endplate_dz,endplate_rmin,endplate_rmax);

  // Start counter material
 

   // Get dimensions of target wall
  vector<double>target_center;
  geom->Get("//posXYZ[@volume='Target']/@X_Y_Z",target_center);
  geom->Get("//tubs[@name='CYLW']/@Rio_Z",targ_wall);
  targ_wall[2]+=target_center[2];
  // Target material
  geom->Get("//tubs[@name='LIH2']/@Rio_Z",target);
  target[2]+=target_center[2];

  ndf=0;

  DEBUG_HISTS=true;
  // DEBUG_HISTS=false;
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
  /*
    double p=fabs(1./q_over_p_);
    double factor=sqrt(1.+tx_*tx_+ty_*ty_);
    mom.SetXYZ(p*tx_/factor,p*ty_/factor,p/factor);
  */
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
  hit->covu=hit->covv=0.000225;

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
  hit->dir=(1./cdchit->wire->udir.z())*cdchit->wire->udir;

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
jerror_t DKalmanFilter::CalcDeriv(double z,DMatrix S, double dEdx, DMatrix &D){
  double x=S(state_x,0), y=S(state_y,0),tx=S(state_tx,0),ty=S(state_ty,0);
  double q_over_p=S(state_q_over_p,0);
  double factor=sqrt(1.+tx*tx+ty*ty);
  double mass=0.14; // pion for now
  double E=sqrt(1./q_over_p/q_over_p+mass*mass); 

  //B-field at (x,y,z)
  double Bx,By,Bz;
  bfield->GetField(x,y,z, Bx, By, Bz);

  D(state_x,0)=tx;
  D(state_y,0)=ty;
  D(state_tx,0)=qBr2p*q_over_p*factor*(tx*ty*Bx-(1.+tx*tx)*By+ty*Bz);
  D(state_ty,0)=qBr2p*q_over_p*factor*((1.+ty*ty)*Bx-tx*ty*By-tx*Bz);
  D(state_q_over_p,0)=-q_over_p*q_over_p*q_over_p*E*dEdx*factor;

  return NOERROR;
}


// Calculate the derivative of the state vector with respect to z and the 
// Jacobian matrix relating the state vector at z to the state vector at z+dz.
jerror_t DKalmanFilter::CalcDerivAndJacobian(double z,DMatrix S,double dEdx,
					     DMatrix &J,DMatrix &D){
  double x=S(state_x,0), y=S(state_y,0),tx=S(state_tx,0),ty=S(state_ty,0);
  double q_over_p=S(state_q_over_p,0);
  double factor=sqrt(1.+tx*tx+ty*ty);
  double mass=0.14; // pion for now
  double E=sqrt(1./q_over_p/q_over_p+mass*mass); 

  //B-field and field gradient at (x,y,z)
  double Bx,By,Bz,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bfield->GetField(x,y,z, Bx, By, Bz);
  bfield->GetFieldGradient(x,y,z,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,
			   dBzdy,dBzdz);

  // Derivative of S with respect to z
  D(state_x,0)=tx;
  D(state_y,0)=ty;
  D(state_tx,0)=qBr2p*q_over_p*factor*(tx*ty*Bx-(1.+tx*tx)*By+ty*Bz);
  D(state_ty,0)=qBr2p*q_over_p*factor*((1.+ty*ty)*Bx-tx*ty*By-tx*Bz);
  D(state_q_over_p,0)=-q_over_p*q_over_p*q_over_p*E*dEdx*factor;

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
  J(state_q_over_p,state_q_over_p)=-dEdx*factor/E
    *(2.+3.*mass*mass*q_over_p*q_over_p);
  
  return NOERROR;
}

// Reference trajectory for forward tracks in CDC region
jerror_t DKalmanFilter::SetCDCForwardReferenceTrajectory(DMatrix &S){
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  DMatrix S0(5,1);
  bool create_new_point=true;

  // first point
  temp.h_id=0;
  temp.num_hits=0;
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z_);
  temp.S= new DMatrix(S);
  temp.Q= new DMatrix(Q);
  temp.J= new DMatrix(J);
  temp.s= len;
  forward_traj_cdc.push_front(temp);

  double dEdx=0.;
  double z=z_,newz,ds=0.;

  // doca variables
  double doca,old_doca;

  // Loop over CDC hits
  for (unsigned int m=0;m<cdchits.size();m++){
    // Position along wire
    DVector3 wirepos=cdchits[m]->origin
      +(z-cdchits[m]->origin.z())*cdchits[m]->dir;
    
    // doca
    old_doca=(temp.pos-wirepos).Mag();

    while(z<endplate_z){
      newz=z+STEP_SIZE;
  
      // Get dEdx for the upcoming step

      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
      
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      GetProcessNoise(0.14,ds,10029.,S,Q);
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory
      temp.s=len;
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),newz);
      if (create_new_point){
	temp.S= new DMatrix(S);	
	temp.Q= new DMatrix(Q);
	temp.J= new DMatrix(J);
	
	forward_traj_cdc.push_front(temp);    
      }
      else{	  
	// In order to prevent unnecessary deleting and creating entries 
	// in the list, we overwrite the values for the trajectory point 
	// corresponding to the true doca and the point immediately after 
	// that.
	for (unsigned int i=0;i<5;i++){
	  forward_traj_cdc[0].S->operator()(i,0)=S(i,0);
	  for (unsigned int j=0;j<5;j++){	      	      
	    forward_traj_cdc[0].Q->operator()(i,j)=Q(i,j);
	    forward_traj_cdc[0].J->operator()(i,j)=J(i,j);
	  }
	}
	forward_traj_cdc[0].pos=temp.pos;
	forward_traj_cdc[0].s=len;
	forward_traj_cdc[0].h_id=temp.h_id;
	create_new_point=true;
      }
      z=newz;
      
      // update position along wire
      wirepos=cdchits[m]->origin
	+(z-cdchits[m]->origin.z())*cdchits[m]->dir;
    
      // new doca 
      doca=(temp.pos-wirepos).Mag();
      
      // Check if we've passed the true minimum doca...
      if (doca>old_doca){	
	if (forward_traj_cdc.size()<3){
	  _DBG_ << "Too few steps before first CDC measurement!" << endl;
	  continue;
	}

	// We've bracketed the minimum; now use the golden section algorithm
	// to find the best doca.  See Numerical Recipes in C, pp 401-402.     
	GoldenSection(newz,STEP_SIZE,dEdx,cdchits[m]->origin,cdchits[m]->dir,S);
	// Trajectory data at point just before minimum doca
	z=forward_traj_cdc[2].pos.z();
	len=forward_traj_cdc[2].s;
	for (unsigned int i=0;i<5;i++)
	  S(i,0)=forward_traj_cdc[2].S->operator()(i,0);
	
	// Step to new position
	ds=Step(z,newz,dEdx,S);
	len+=ds;
	
	// Get the contribution to the covariance matrix due to multiple 
	// scattering
	GetProcessNoise(0.14,ds,10029.,S,Q);

	// Compute the Jacobian matrix
	StepJacobian(newz,z,S,dEdx,J);

	// update the trajectory
	z=newz;
	forward_traj_cdc[1].h_id=m+1;
	forward_traj_cdc[1].s=len;
	forward_traj_cdc[1].pos.SetXYZ(S(state_x,0),S(state_y,0),z);
	for (unsigned int i=0;i<5;i++){
	  forward_traj_cdc[1].S->operator()(i,0)=S(i,0);
	  for (unsigned int j=0;j<5;j++){	      
	    forward_traj_cdc[1].J->operator()(i,j)=J(i,j);
	    forward_traj_cdc[1].Q->operator()(i,j)=Q(i,j);
	  }
	}

	create_new_point=false;
	break;
      }
      old_doca=doca;
    }
  }
  /*
  for (unsigned int m=0;m<forward_traj_cdc.size();m++){
    printf("id %d x %f y %f z %f s %f \n",
	   forward_traj_cdc[m].h_id,forward_traj_cdc[m].pos.x(),
	   forward_traj_cdc[m].pos.y(),forward_traj_cdc[m].pos.z(),
	   forward_traj_cdc[m].s);
  }    
  */

  // position at the end of the swim
  z_=z;
  x_=S(state_x,0);
  y_=S(state_y,0);

  return NOERROR;
}

// Reference trajectory for central tracks
jerror_t DKalmanFilter::SetCDCReferenceTrajectory(DMatrix &Sc){
  DMatrix Jc(5,5),Q(5,5);    
  DKalmanState_t temp;
  DMatrix S0(5,1);

  // Position, step, radius, etc. variables
  DVector3 pos(x_,y_,z_),origin(0,0,65.),dir(0,0,1.),wirepos; 
  DVector3 old_origin=origin,old_dir=dir;
  double r,q_over_p;
  double dedx=0,ds=STEP_SIZE;
  bool create_new_point=true;
  double doca=0.,old_doca=0.;
  double sign=1.;

  // First point at "vertex"
  temp.h_id=temp.num_hits=0;
  temp.pos=pos;
  temp.S=new DMatrix(Sc);
  temp.Q=new DMatrix(Q);
  temp.J=new DMatrix(Jc);
  central_traj.push_front(temp);
 
  // Loop over CDC hits
  for (unsigned int m=0;m<cdchits.size();m++){
    // Reset ds
    ds=STEP_SIZE;
    
    //B-field at (x,y,z)
    double Bx,By,Bz;
    bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
    double B=sqrt(Bx*Bx+By*By+Bz*Bz);
    double Rc=1./qBr2p/B/Sc(state_q_over_pt,0);
    double cosphi=cos(Sc(state_phi,0));
    double sinphi=sin(Sc(state_phi,0));
    double xc=pos.x()-Rc*sinphi;
    double yc=pos.y()+Rc*cosphi;
    //double rc2=(wirepos.x()-xc)*(wirepos.x()-xc)
    //  +(wirepos.y()-yc)*(wirepos.y()-yc);
    double rc2=xc*xc+yc*yc;
    sign=Sc(state_q_over_pt,0)/fabs(Sc(state_q_over_pt,0));
    if (rc2<Rc*Rc) sign*=-1.;

    Sc(state_D,0)=sign*pos.Perp();

    // position on wire
    wirepos=cdchits[m]->origin
      +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir;
    
    // doca
    old_doca=(pos-wirepos).Mag();

    dedx=0.; //temp

    while(true){       
      // update path length
      len+=ds;
      
      // Propagate the state through the field
      Step(pos,origin,dir,ds,Sc,dedx);
      
      // Compute the Jacobian matrix
      StepJacobian(pos,origin,dir,-ds,Sc,dedx,Jc);
      
      // Process noise
      q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
      if (fabs(Sc(state_D,0))<0.8){
	GetProcessNoiseCentral(0.14,ds,10029.,Sc,Q); // use Ar for now
      }
      else{
	GetProcessNoiseCentral(0.14,ds,30420.,Sc,Q); // use Air for now
      }
   
      // Position along wire
      wirepos=cdchits[m]->origin
	+(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir;
    
      // new doca 
      doca=(pos-wirepos).Mag(); 

      Sc(state_D,0)=sign*pos.Perp();

      temp.pos=pos;	
      temp.s=len;
      if (create_new_point){
	temp.S= new DMatrix(Sc);	
	temp.Q= new DMatrix(Q);
	temp.J= new DMatrix(Jc);
	
	central_traj.push_front(temp);    
      }
      else{	  
	// In order to prevent unnecessary deleting and creating entries 
	// in the list, we overwrite the values for the trajectory point 
	// corresponding to the true doca and the point immediately after 
	// that.
	for (unsigned int i=0;i<5;i++){
	  central_traj[0].S->operator()(i,0)=Sc(i,0);
	  for (unsigned int j=0;j<5;j++){	      	      
	    central_traj[0].Q->operator()(i,j)=Q(i,j);
	    central_traj[0].J->operator()(i,j)=Jc(i,j);
	  }
	}
	central_traj[0].pos=pos;
	central_traj[0].s=len;
	central_traj[0].h_id=temp.h_id;
	create_new_point=true;
      }

      // Check if we've passed the true minimum doca...
      // if (fabs(doca)<fabs(Sc(state_D,0) )
      if (doca>old_doca)
	  {
	double ds1=central_traj[0].s-central_traj[1].s;

	// We've bracketed the minimum; now use the golden section algorithm
	// to find the best doca.  See Numerical Recipes in C, pp 401-402.     
	ds=GoldenSection(ds1,ds,dedx,pos,cdchits[m]->origin,cdchits[m]->dir,Sc);
	
	// Trajectory data at point just before minimum doca
	for (unsigned int i=0;i<5;i++)
	  Sc(i,0)=central_traj[2].S->operator()(i,0);
	pos=central_traj[2].pos;
	  
	// Energy loss and process noise (multiple scattering)
	q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
	if (fabs(doca)<0.8){
	  GetProcessNoiseCentral(0.14,ds,10029.,Sc,Q); // use Ar for now
	  dedx=GetdEdx(0.14,q_over_p,18.,39.9,0.00166);
	}
	else{
	  GetProcessNoiseCentral(0.14,ds,30420.,Sc,Q); // use Air for now
	  dedx=GetdEdx(0.14,q_over_p,7.,14.,1.205e-3);
	}
	dedx=0.;
	
	//pos.Print();
	// Step to the minimum doca 
	Step(pos,origin,dir,ds,Sc,dedx);
	//pos.Print();

	Sc(state_D,0)=sign*pos.Perp();

	
	// Compute the Jacobian matrix
	StepJacobian(pos,origin,dir,-ds,Sc,dedx,Jc);
	
	// Update the trajectory
	for (unsigned int i=0;i<5;i++){
	  central_traj[1].S->operator()(i,0)=Sc(i,0);
	  for (unsigned int j=0;j<5;j++){	      
	    central_traj[1].J->operator()(i,j)=Jc(i,j);
	    central_traj[1].Q->operator()(i,j)=Q(i,j);
	  }
	}
	len=central_traj[1].s=central_traj[2].s+ds;
	central_traj[1].pos=pos;
	central_traj[1].h_id=m+1;

	create_new_point=false;
	
	break;
      } 
	
      // Save the doca
      old_doca=doca;
      //doca=Sc(state_D,0);
    }
    old_origin=origin;
    old_dir=dir;
  }

  // Swim out
  r=pos.Perp();
  ds=STEP_SIZE;
  while (r<endplate_rmax && pos.z()<endplate_z && pos.z()>0.){
    // update path length
    len+=ds;

    // Propagate the state through the field
    Step(pos,origin,dir,ds,Sc,dedx);

    Sc(state_D,0)=sign*pos.Perp();

    // Compute the Jacobian matrix
    StepJacobian(pos,origin,dir,-ds,Sc,dedx,Jc);
      
    // Process noise
    q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
    if (fabs(Sc(state_D,0))<0.8){
      GetProcessNoiseCentral(0.14,ds,10029.,Sc,Q); // use Ar for now
    }
    else{
      GetProcessNoiseCentral(0.14,ds,30420.,Sc,Q); // use Air for now
    }
    
    // Update the trajectory
    temp.pos=pos;	
    temp.s=len;
    temp.S= new DMatrix(Sc);	
    temp.Q= new DMatrix(Q);
    temp.J= new DMatrix(Jc);
    central_traj.push_front(temp);    

    r=pos.Perp();
  }

  /*
    for (unsigned int m=0;m<central_traj.size();m++){
    printf("id %d x %f y %f z %f s %f d %f\n",
    central_traj[m].h_id,central_traj[m].pos.x(),
    central_traj[m].pos.y(),central_traj[m].pos.z(),central_traj[m].s,
    central_traj[m].S->operator()(state_D,0));
    }
  */

  // Position at the end of the swim
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  return NOERROR;
}

// Reference trajectory for trajectories with hits in the forward direction
jerror_t DKalmanFilter::SetReferenceTrajectory(DMatrix &S){ 
  // Order the hits
  sort(fdchits.begin(),fdchits.end(),DKalmanFDCHit_cmp);
  
  DMatrix J(5,5),Q(5,5);    
  DKalmanState_t temp;
  double ds=0.; // path length increment

  // first point
  temp.h_id=0;
  temp.num_hits=0;
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),z_);

  if (forward_traj_cdc.size()>0){
    J=*(forward_traj_cdc[0].J);
    Q=*(forward_traj_cdc[0].Q);
  }
  temp.S= new DMatrix(S);
  temp.Q= new DMatrix(Q);
  temp.J= new DMatrix(J);
  temp.s= len;
  forward_traj.push_front(temp);

  double dEdx=0.;

   // progress in z from hit to hit
  double z=z_;
  double newz;  
  for (unsigned int m=0;m<fdchits.size();m++){
    int num=int((fdchits[m]->z-z)/STEP_SIZE);
    newz=fdchits[m]->z-STEP_SIZE*double(num);
    
    if (newz-z>0.){
      // Get dEdx for the upcoming step

      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
      
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      GetProcessNoise(0.14,ds,10029.,S,Q);
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory
      temp.s=len;
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),newz);
      temp.S=new DMatrix(S);
      temp.Q=new DMatrix(Q);
      temp.J=new DMatrix(J);
      forward_traj.push_front(temp);

      // update z
      z=newz;
    }

    for (int k=0;k<num;k++){
      newz=z+STEP_SIZE;

      // Step through field
      ds=Step(z,newz,dEdx,S);
      len+=ds;
      
      // Get the contribution to the covariance matrix due to multiple 
      // scattering
      GetProcessNoise(0.14,ds,10029.,S,Q);
      
      // Compute the Jacobian matrix
      StepJacobian(newz,z,S,dEdx,J);
      
      // update the trajectory
      temp.s=len;
      temp.pos.SetXYZ(S(state_x,0),S(state_y,0),newz);
      temp.S=new DMatrix(S);
      temp.Q=new DMatrix(Q);
      temp.J=new DMatrix(J);
      forward_traj.push_front(temp);
      
      //update z
      z=newz;
    }
    // flag this as an active layer
    forward_traj[0].h_id=m+1;
  }
  // One more step after last hit point
  newz=z+STEP_SIZE;
  ds=Step(z,newz,dEdx,S);
  len+=ds;
      
  // Get the contribution to the covariance matrix due to multiple 
  // scattering
  GetProcessNoise(0.14,ds,10029.,S,Q);
      
  // Compute the Jacobian matrix
  StepJacobian(newz,z,S,dEdx,J);
      
  // update the trajectory
  temp.s=len;
  temp.pos.SetXYZ(S(state_x,0),S(state_y,0),newz);
  temp.S=new DMatrix(S);
  temp.Q=new DMatrix(Q);
  temp.J=new DMatrix(J);
  forward_traj.push_front(temp);

  /*
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



// Step the state vector through the field from oldz to newz.
// Uses the 4th-order Runga-Kutte algorithm.
double DKalmanFilter::Step(double oldz,double newz, double dEdx,DMatrix &S){
  double delta_z=newz-oldz;
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);

  double s=sqrt(1.+S(state_tx,0)*S(state_tx,0)+S(state_ty,0)*S(state_ty,0))
    *delta_z;
  CalcDeriv(oldz,S,dEdx,D1);
  CalcDeriv(oldz+delta_z/2.,S+0.5*delta_z*D1,dEdx,D2);
  CalcDeriv(oldz+delta_z/2.,S+0.5*delta_z*D2,dEdx,D3);
  CalcDeriv(oldz+delta_z,S+delta_z*D3,dEdx,D4);
	
  S+=delta_z*((1./6.)*D1+(1./3.)*D2+(1./3.)*D3+(1./6.)*D4);

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
  CalcDerivAndJacobian(oldz,S,dEdx,J1,D1);
  CalcDerivAndJacobian(oldz+delta_z/2.,S+0.5*delta_z*D1,dEdx,J2,D2);
  J2=J2+0.5*(J2*J1);
  CalcDerivAndJacobian(oldz+delta_z/2.,S+0.5*delta_z*D2,dEdx,J3,D3);
  J3=J3+0.5*(J3*J2);
  CalcDerivAndJacobian(oldz+delta_z,S+delta_z*D3,dEdx,J4,D4);
  J4=J4+J4*J3;

  J+=delta_z*((1./6.)*J1+(1./3.)*J2+(1./3.)*J3+(1./6.)*J4);
  
  return NOERROR;
}

// Calculate the derivative for the alternate set of parameters {q/pT, phi, 
// tan(lambda),D,z}
jerror_t DKalmanFilter::CalcDeriv(double ds,DVector3 pos,DVector3 &dpos,
				  DVector3 wire_orig,
				  DVector3 wiredir,
				  DMatrix S,double dEdx,
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
  double D=S(state_D,0);
  double sign=D/fabs(D);
  double z=S(state_z,0);
  double p=pt/cosl;
  double mass=0.14; //pion for now
  DVector3 wire_pos=wire_orig+(z-wire_orig.z())*wiredir;
  double E=sqrt(p*p+mass*mass);
  double dx=pos.x()-wire_pos.x();
  double dy=pos.y()-wire_pos.y();

  //B-field and field gradient at (x,y,z)
  double Bx,By,Bz,B;
  bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
  B=sqrt(Bx*Bx+By*By+Bz*Bz);

  // Derivative of S with respect to s
  D1(state_q_over_pt,0)=qBr2p*q_over_pt*q_over_pt*sinl*(By*cosphi-Bx*sinphi)
    -q_over_pt*E/p/p*dEdx;
  D1(state_phi,0)=qBr2p*q_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_tanl,0)=qBr2p*q_over_pt*(By*cosphi-Bx*sinphi)/cosl;

  if (fabs(D)>EPS2)
    D1(state_D,0)=(dx*(cosphi*cosl-wiredir.x()*sinl)
		   +dy*(sinphi*cosl-wiredir.y()*sinl))/D;

  if (sign*(D+D1(state_D,0)*ds)<0)
    D1(state_D,0)*=-1.;

  D1(state_z,0)=sinl;

  // New direction
  dpos.SetXYZ(cosl*cosphi,cosl*sinphi,sinl);
  // Include second order correction 
  /*
  dpos.SetXYZ(cosl*(cosphi+sinphi*ds*qBr2p*q_over_pt*B),
	      cosl*(sinphi-cosphi*ds*qBr2p*q_over_pt*B),sinl);
  */

  return NOERROR;
}

// Calculate the derivative and Jacobian matrices for the alternate set of 
// parameters {q/pT, phi, tan(lambda),D,z}
jerror_t DKalmanFilter::CalcDerivAndJacobian(double ds,DVector3 pos,
					     DVector3 &dpos,
					     DVector3 wire_orig,
					     DVector3 wiredir,
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
  double D=S(state_D,0);
  double z=S(state_z,0);
  double p=pt/cosl;
  double q=pt*q_over_pt;
  double mass=0.14; //pion for now
  DVector3 wire_pos=wire_orig+(z-wire_orig.z())*wiredir;
  double E=sqrt(p*p+mass*mass);
  double dx=pos.x()-wire_pos.x();
  double dy=pos.y()-wire_pos.y();

  //B-field and field gradient at (x,y,z)
  double Bx,By,Bz,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
  bfield->GetFieldGradient(pos.x(),pos.y(),pos.z(),dBxdx,dBxdy,dBxdz,dBydx,
			   dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  
  // New direction
  dpos.SetXYZ(cosl*cosphi,cosl*sinphi,sinl);
  // Include second order correction 
  //dpos.SetXYZ(cosl*(cosphi+sinphi*ds*qBr2p*q_over_pt*B),
  //cosl*(sinphi-cosphi*ds*qBr2p*q_over_pt*B),sinl);


  // Derivative of S with respect to s
  D1(state_q_over_pt,0)=qBr2p*q_over_pt*q_over_pt*sinl*(By*cosphi-Bx*sinphi)
    -q_over_pt*E/p/p*dEdx;

  D1(state_phi,0)=qBr2p*q_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_tanl,0)=qBr2p*q_over_pt*(By*cosphi-Bx*sinphi)/cosl;

  if (fabs(D)>EPS2){
    D1(state_D,0)=(dx*(dpos.X()-wiredir.x()*sinl)
		   +dy*(dpos.Y()-wiredir.y()*sinl))/D;
  }
  D1(state_z,0)=sinl;

  // Jacobian matrix elements
  J1(state_phi,state_phi)=qBr2p*q_over_pt*sinl*(By*cosphi-Bx*sinphi);
  J1(state_phi,state_q_over_pt)=qBr2p*(Bx*cosphi*sinl+By*sinphi*sinl
				       -Bz*cosl);
  J1(state_phi,state_tanl)=qBr2p*q_over_pt*(Bx*cosphi*cosl+By*sinphi*cosl
					    +Bz*sinl)/(1.+tanl*tanl);
  J1(state_phi,state_z)
    =qBr2p*q_over_pt*(dBxdz*cosphi*sinl+dBydz*sinphi*sinl-dBzdz*cosl);

  J1(state_tanl,state_phi)=-qBr2p*q_over_pt*(By*sinphi+Bx*cosphi)/cosl;
  J1(state_tanl,state_q_over_pt)=qBr2p*(By*cosphi-Bx*sinphi)/cosl;
  J1(state_tanl,state_tanl)=qBr2p*q_over_pt*sinl*(By*cosphi-Bx*sinphi);
  J1(state_tanl,state_z)=qBr2p*q_over_pt*(dBydz*cosphi-dBxdz*sinphi)/cosl;
  
  J1(state_q_over_pt,state_phi)
    =-qBr2p*q_over_pt*q_over_pt*sinl*(By*sinphi+Bx*cosphi);

  double temp=sqrt(1.+cosl*cosl*q_over_pt*q_over_pt*mass*mass);
  J1(state_q_over_pt,state_q_over_pt)
    =2.*qBr2p*q_over_pt*sinl*(By*cosphi-Bx*sinphi)
    -q*dEdx*cosl*q_over_pt*(2.+3.*cosl*cosl*q_over_pt*q_over_pt*mass*mass)
    /temp;
  J1(state_q_over_pt,state_tanl)
    =qBr2p*q_over_pt*q_over_pt*cosl*cosl*cosl*(By*cosphi-Bx*sinphi)
    +q*dEdx*sinl*cosl*cosl/temp*q_over_pt*q_over_pt
    *(1.+2.*cosl*cosl*mass*mass*q_over_pt*q_over_pt);
  J1(state_q_over_pt,state_z)
    =qBr2p*q_over_pt*q_over_pt*sinl*(dBydz*cosphi-dBxdz*sinphi);

  J1(state_z,state_tanl)=cosl/(1.+tanl*tanl);

  if (fabs(D)>EPS2){
    J1(state_D,state_q_over_pt)=(dx*cosl*sinphi*ds*qBr2p*B-dy*cosl*cosphi*ds*qBr2p*B)/D;

    J1(state_D,state_phi)=cosl*(dy*(cosphi+sinphi*ds*qBr2p*q_over_pt*B)
				-dx*(sinphi-cosphi*ds*qBr2p*q_over_pt*B))/D;
    
    J1(state_D,state_tanl)
      =-1./D/(1+tanl*tanl)*(dx*((cosphi+sinphi*ds*qBr2p*q_over_pt*B)*sinl+wiredir.x()*cosl)  
			    +dy*((sinphi-cosphi*ds*qBr2p*q_over_pt*B)*sinl+wiredir.y()*cosl));

    J1(state_D,state_D)=-D1(state_D,0)/D;

    // printf("dx %f dy %f D %f J %f D1 %f \n",dx,dy,D,J1(state_D,state_D),D1(state_D,0));

    if (wiredir.x()!=0. || wiredir.y()!=0.)
      J1(state_D,state_z)=(-1./D)*(wiredir.x()*(dpos.X()-wiredir.x()*sinl)
				   +wiredir.y()*(dpos.Y()-wiredir.y()*sinl));
   
  }

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
jerror_t DKalmanFilter::Step(DVector3 &pos, DVector3 wire_orig,DVector3 wiredir,
				     double ds,DMatrix &S,double dEdx){
  // Matrices for intermediate steps
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  DMatrix S1(5,1),S2(5,1),S3(5,1),S4(5,1);
  DVector3 dpos1,dpos2,dpos3,dpos4;

  CalcDeriv(0.,pos,dpos1,wire_orig,wiredir,S,dEdx,D1);

  DVector3 mypos=pos+(ds/2.)*dpos1;
  S1=S+(0.5*ds)*D1; 

  CalcDeriv(ds/2.,mypos,dpos2,wire_orig,wiredir,S1,dEdx,D2);

  mypos=pos+(ds/2.)*dpos2;
  S2=S+(0.5*ds)*D2; 

  CalcDeriv(ds/2.,mypos,dpos3,wire_orig,wiredir,S2,dEdx,D3);

  mypos=pos+(ds/2.)*dpos3;
  S3=S+ds*D3;

  CalcDeriv(ds,mypos,dpos4,wire_orig,wiredir,S3,dEdx,D4);
   
  // New state vector
  S+=ds*((1./6.)*D1+(1./3.)*D2+(1./3.)*D3+(1./6.)*D4);
  
  // New position
  pos+=ds*((1./6.)*dpos1+(1./3.)*dpos2+(1./3.)*dpos3+(1./6.)*dpos4);

  return NOERROR;
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

  CalcDerivAndJacobian(0.,pos,dpos1,wire_orig,wiredir,S,dEdx,J1,D1);

  DVector3 mypos=pos+(ds/2.)*dpos1;
  S1=S+(0.5*ds)*D1; 

  CalcDerivAndJacobian(ds/2.,mypos,dpos2,wire_orig,wiredir,S1,dEdx,J2,D2);
  J2=J2+0.5*(J2*J1);

  mypos=pos+(ds/2.)*dpos2;
  S2=S+(0.5*ds)*D2;

  CalcDerivAndJacobian(ds/2.,mypos,dpos3,wire_orig,wiredir,S2,dEdx,J3,D3);
  J3=J3+0.5*(J3*J2);  

  mypos=pos+(ds/2.)*dpos3;
  S3=S+(ds/2.)*D3;

  CalcDerivAndJacobian(ds,mypos,dpos4,wire_orig,wiredir,S3,dEdx,J4,D4);
  J4=J4+J4*J3;

  // New Jacobian matrix
  J+=ds*((1./6.)*J1+(1./3.)*J2+(1./3.)*J3+(1./6.)*J4);

  return NOERROR;
}

// Multiple scattering covariance matrix for central track parameters
jerror_t DKalmanFilter::GetProcessNoiseCentral(double mass_hyp,double ds,
					       double X0,DMatrix Sc,
					       DMatrix &Q){
  DMatrix Q1(5,5);
  double tanl=Sc(state_tanl,0);
  double q_over_pt=Sc(state_q_over_pt,0);
 
  Q1(state_phi,state_phi)=1.+tanl*tanl;
  Q1(state_tanl,state_tanl)=(1.+tanl*tanl)*(1.+tanl*tanl);
  Q1(state_q_over_pt,state_q_over_pt)=q_over_pt*q_over_pt*tanl*tanl;
  Q1(state_q_over_pt,state_tanl)=Q1(state_tanl,state_q_over_pt)
    =q_over_pt*tanl*(1.+tanl*tanl);

  double my_ds=fabs(ds);
  double p2=(1.+tanl*tanl)/q_over_pt/q_over_pt;
  double sig2_ms= 0.0136*0.0136*(1.+mass_hyp*mass_hyp/p2)
    *my_ds/X0/p2*(1.+0.038*log(my_ds/X0))*(1.+0.038*log(my_ds/X0));

  sig2_ms=0.;

  Q=sig2_ms*Q1;

  return NOERROR;

}
					       

// Compute contributions to the covariance matrix due to multiple scattering
jerror_t DKalmanFilter::GetProcessNoise(double mass_hyp,double ds,
					double X0,DMatrix S,DMatrix &Q){
  DMatrix Q1(5,5);
  double tx=S(state_tx,0),ty=S(state_ty,0);
  double one_over_p_sq=S(state_q_over_p,0)*S(state_q_over_p,0);
  double my_ds=fabs(ds);

  Q1(state_tx,state_tx)=(1.+tx*tx)*(1.+tx*tx+ty*ty);
  Q1(state_ty,state_ty)=(1.+ty*ty)*(1.+tx*tx+ty*ty);
  Q1(state_tx,state_ty)=Q1(state_ty,state_tx)=tx*ty*(1.+tx*tx+ty*ty);

  double sig2_ms= 0.0136*0.0136*(1.+one_over_p_sq*mass_hyp*mass_hyp)
    *one_over_p_sq*my_ds/X0*(1.+0.038*log(my_ds/X0))*(1.+0.038*log(my_ds/X0));

  sig2_ms=0.;

  Q=sig2_ms*Q1;

  return NOERROR;
}

// Calculate the energy loss per unit length given properties of the material
// through which a particle of momentum p is passing
double DKalmanFilter::GetdEdx(double M,double q_over_p,double Z,
			      double A, double rho){

  double betagamma=1./M/fabs(q_over_p);
  double beta2=1./(1.+M*M*q_over_p*q_over_p);
  double Me=511.; //keV
  double m_ratio=Me*1.e-6/M;
  double Tmax
    =2.*Me*betagamma*betagamma/(1.+2.*sqrt(1.+betagamma*betagamma)*m_ratio
				+m_ratio*m_ratio);
  double I0=(12.*Z+7.)*1e-3; //keV

  return 0.;

  if (Z>12) I0=(9.76*Z+58.8*pow(Z,-0.19))*1e-3;
  return -0.0001535*Z/A*rho/beta2*(log(2.*Me*betagamma*betagamma*Tmax/I0/I0)
			       -2.*beta2);
}

// Swim the state vector from z_start to z_end through the field
jerror_t DKalmanFilter::SwimToPlane(double z_start,double z_end, DMatrix &S){
  int num_inc=(int)((z_end-z_start)/STEP_SIZE);
  double z=z_start;
  double dedx=0.;
  for (int i=0;i<num_inc;i++){
    z_=z+STEP_SIZE;  
    dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
    dedx=0.;
    Step(z,z_,dedx,S);  
    z=z_;
  }
  // Final step
  dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
  dedx=0.;
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
    dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
    StepJacobian(z,z_,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);     
    C=J*(C*J_T);
    z=z_;
    printf("z %f\n",z);
  }
  // Final step
  dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
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
    double tanl=Sc(state_tanl,0);
    double cosl=cos(atan(tanl));
    double q_over_p=Sc(state_q_over_pt,0)*cosl;
    dedx=GetdEdx(0.14,q_over_p,18.,39.9,0.00166);

   // Step the position, state vector, and covariance matrix through the field
    StepJacobian(pos,cdchits[m]->origin,cdchits[m]->dir,ds,Sc,dedx,Jc);
    Cc=Jc*(Cc*DMatrix(DMatrix::kTransposed,Jc));

    Step(pos,cdchits[m]->origin,cdchits[m]->dir,ds,Sc,dedx);

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
      Step(newpos,cdchits[m]->origin,cdchits[m]->dir,ds,S0,dedx);
      R=newpos.Perp();
      s1+=ds;
      num_iter++;
    }
    // Final step
    StepJacobian(pos,cdchits[m]->origin,cdchits[m]->dir,s1,Sc,dedx,Jc);
    Cc=Jc*(Cc*DMatrix(DMatrix::kTransposed,Jc));

    Step(pos,cdchits[m]->origin,cdchits[m]->dir,s1,Sc,dedx);
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


// Interface routine for Kalman filter
jerror_t DKalmanFilter::KalmanLoop(double mass_hyp){
  if (z_<0) return VALUE_OUT_OF_RANGE;

  DMatrix S(5,1),C(5,5),Sc(5,1),Cc(5,5);
  DMatrix C0(5,5);
  DMatrix Sbest(5,1),Cbest(5,5),Scbest(5,1),Ccbest(5,5);
  double chisq=0.,chisq_forward=0.,chisq_central=0.;
  chisq_=1.e8;
  // position along track. 
  DVector3 pos(x_,y_,z_); 

  // Vector of cdc residuals  
  vector<double>best_cdc_resid(cdchits.size()); 
  vector<double>best_cdc_pulls(cdchits.size());

  // Initialize path length variable
  len=0;
 
  // Deal with hits in the CDC 
  if (cdchits.size()>0 && fdchits.size()==0){
     // Order the CDC hits by radius
    sort(cdchits.begin(),cdchits.end(),DKalmanCDCHit_cmp);

    // Initialize the state vector and covariance matrix
    Sc(state_q_over_pt,0)=q_over_pt_;
    Sc(state_phi,0)=phi_;
    Sc(state_tanl,0)=tanl_;
    Sc(state_z,0)=z_;  

    //B-field at (x,y,z)
    double Bx,By,Bz;
    bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
    double B=sqrt(Bx*Bx+By*By+Bz*Bz);
    double Rc=1./qBr2p/B/Sc(state_q_over_pt,0);
    double cosphi=cos(Sc(state_phi,0));
    double sinphi=sin(Sc(state_phi,0));
    double xc=pos.x()-Rc*sinphi;
    double yc=pos.y()+Rc*cosphi;
    DVector3 wirepos=cdchits[0]->origin+(pos.z()-cdchits[0]->origin.z())*
      cdchits[0]->dir;
    double rc2=(wirepos.x()-xc)*(wirepos.x()-xc)
      +(wirepos.y()-yc)*(wirepos.y()-yc);
    double sign=Sc(state_q_over_pt,0)/fabs(Sc(state_q_over_pt,0));
    if (rc2<Rc*Rc) sign*=-1.;    
    Sc(state_D,0)=sign*(pos-wirepos).Perp();   

    SetCDCReferenceTrajectory(Sc);
    
    // So far the candidate is our best result for the state vector...
    Scbest=Sc;
    Ccbest=Cc;
    
    // Compute track quantities in the forward representation if we have FDC
    // hits.
    if (hits.size()>0){
      tx_=cos(Sc(state_phi,0))/Sc(state_tanl,0);
      ty_=sin(Sc(state_phi,0))/Sc(state_tanl,0);
      q_over_p_=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
    }
  }

  // deal with hits in FDC
  if (fdchits.size()>0){   
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

    // If we have cdc hits, swim through the field past these measurements
    // first
    if (cdchits.size()>0){
      // Order the CDC hits by radius
      sort(cdchits.begin(),cdchits.end(),DKalmanCDCHit_cmp);

      SetCDCForwardReferenceTrajectory(S);
    }

    // Swim once through the field out to the most upstream FDC hit
    SetReferenceTrajectory(S);
 
    unsigned int num_iter=NUM_ITER;
    last_iter=false;
    num_iter=3;
    for (unsigned int iter=0;iter<num_iter;iter++) {
       // perform the kalman filter 

      if (iter>0){
	// Swim back to the first (most downstream) plane and use the new 
      // values of S and C as the seed data to the Kalman filter 
	SwimToPlane(z_,forward_traj[0].pos.Z(),S);
      } 
  
      if (iter==num_iter-1) last_iter=true;
      C=C0;
      KalmanForward(mass_hyp,S,C,chisq_forward);
    
      if (cdchits.size()>0){
	double f=2.5;
	double anneal=1.;
	anneal=99./pow(f,iter)+1.;

	// Proceed into CDC
	KalmanForwardCDC(mass_hyp,anneal,S,C,chisq_forward);
      }
    } // iteration  
    
    // Find the state at the so-called "vertex" position
    ExtrapolateToVertex(mass_hyp,S,C);
    
    // total chisq and ndf
    chisq_=chisq_forward;
    ndf=fdchits.size()+cdchits.size();

    return NOERROR;
  }

  // Fit in Central region
  if (false){
  //if (cdchits.size()>0){   
    // Guess for covariance matrix in central representation
    Cc(state_z,state_z)=0.25;
    Cc(state_q_over_pt,state_q_over_pt)=0.04*q_over_pt_*q_over_pt_;
    Cc(state_phi,state_phi)=0.025*0.025;
    Cc(state_D,state_D)=0.25;
    double theta=90.-180./M_PI*atan(tanl_);
    double dlambda=2.1e-3-3.4e-4*theta+3.5e-5*theta*theta;
    dlambda=0.125;
    Cc(state_tanl,state_tanl)=(1.+tanl_*tanl_)*(1.+tanl_*tanl_)
      *dlambda*dlambda;
    
    // position vector
    DVector3 pos=central_traj[0].pos;
      
    C0=Cc;

    // Calculate a scale factor for the measurement errors that depends 
    // on the iteration,so that we approach the "true' measurement errors
    // by the last iteration.
    double f=2.5;
    double anneal_factor=1.;
    //anneal_factor=499./pow(f,iter)+1.;
    //if (iter>0) anneal_factor=1.;
    
    //Cc=C0;
    jerror_t error=NOERROR;
    error=KalmanCentral(mass_hyp,anneal_factor,Sc,Cc,pos,chisq_central);
    //if (error!=NOERROR) break;
    
    //printf("chi2 %f p %f theta %f\n",chisq,1./Sc(state_q_over_pt,0)/cos(atan(Sc(state_tanl,0))),M_PI/2.-atan(Sc(state_tanl,0)));
    //printf("iter %d chi2*anneal %f\n",iter,chisq*anneal_factor);
    
    if (chisq_central>=1e8 ){
      _DBG_ << "-- central fit failed --" <<endl;
      return VALUE_OUT_OF_RANGE;
    }
    
    /*
      if (DEBUG_HISTS){
      TH2F *cdc_residuals=(TH2F*)gROOT->FindObject("cdc_residuals");
      TH2F *cdc_pulls_histo=(TH2F*)gROOT->FindObject("cdc_pulls");
      if (cdc_residuals && cdc_pulls_histo)
      for (unsigned int m=0;m<cdchits.size();m++){
      cdc_residuals->Fill(cdchits[m]->origin.Perp(),best_cdc_resid[m]);
      cdc_pulls_histo->Fill(cdchits[m]->origin.Perp(),best_cdc_pulls[m]);
      }
      }
    */ 
  
    // Find track parameters where track crosses beam line
    ExtrapolateToVertex(mass_hyp,pos,Sc,Cc);

    // total chisq and ndf
    chisq_=chisq_central;
    ndf=cdchits.size();
  }

  return NOERROR;
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
  Step(pos,origin,dir,s1,S0,dedx);
  DVector3 wirepos=origin+(pos.z()-origin.z())*dir;
  d1=(pos-wirepos).Mag();

  // Step to the other intermediate point
  S0=Sc; 
  pos=oldpos;
  Step(pos,origin,dir,s2,S0,dedx);
  wirepos=origin+(pos.z()-origin.z())*dir;
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
      Step(pos,origin,dir,s2,S0,dedx);
      wirepos=origin+(pos.z()-origin.z())*dir;
      d2=(pos-wirepos).Mag();
    }
    else{
      stemp=golden_ratio1*s2+golden_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      Step(pos,origin,dir,s1,S0,dedx);
      wirepos=origin+(pos.z()-origin.z())*dir;
      d1=(pos-wirepos).Mag();
    }
  }
  if (d1<d2) {
    return s1+ds;
  }
  else{
    return s2+ds;
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
  Step(pos,origin,dir,s1,S0,dedx);
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
      Step(pos,origin,dir,s2,S0,dedx);
      d2=S0(state_D,0);
    }
    else{
      stemp=magic_ratio1*s2+magic_ratio2*s0;
      s3=s2;
      s2=s1;
      s1=stemp;
      d2=d1;
      Step(pos,origin,dir,s1,S0,dedx);
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
  //double V=0.8*0.8/12.;  // Measurement variance
  double V=anneal_factor*0.000225;
  double InvV; // inverse of variance
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix S0(5,1),S0_(5,1); // state vector
  DMatrix M(1,1);
  
  // Initialize the chi2 for this part of the track
  chisq=0.;

  //Track projection matrix
  H(0,state_D)=H_T(state_D,0)=1.;

  S0_=DMatrix(*central_traj[0].S);
  for (unsigned int k=1;k<central_traj.size()-1;k++){
    // Get the state vector, jacobian matrix, and multiple scattering matrix 
    // from reference trajectory
    S0=DMatrix(*central_traj[k].S);
    J=DMatrix(*central_traj[k].J);
    Q=DMatrix(*central_traj[k].Q);

    double sign=S0(state_D,0)/fabs(S0(state_D,0));
    double sign2=S0_(state_D,0)/fabs(S0_(state_D,0));
    if (sign!=sign2) S0_(state_D,0)*=-1.;
    Sc(state_D,0)=sign*fabs(Sc(state_D,0));

    // State S is perturbation about a seed S0
    dS=Sc-S0_;

    //dS.Print();
    //Cc.Print();
    
    // Update the actual state vector and covariance matrix
    Sc=S0+J*dS;
    Cc=J*(Cc*DMatrix(DMatrix::kTransposed,J))+Q;   

    // update position
    double dz=Sc(state_z,0)-S0(state_z,0);
    pos(0)=central_traj[k].pos.x()+dz*cos(S0(state_phi,0))/S0(state_tanl,0);
    pos(1)=central_traj[k].pos.y()+dz*sin(S0(state_phi,0))/S0(state_tanl,0);
    pos(2)=Sc(state_z,0);

    // Save the current state of the reference trajectory
    S0_=S0;

    // Add the hit
    if (central_traj[k].h_id>0){
      unsigned int cdc_index=central_traj[k].h_id-1;
      double sign=Sc(state_D,0)/fabs(Sc(state_D,0));
      if (cdchits[cdc_index]->stereo!=0.)
	H(0,state_D)=H_T(state_D,0)=sign*cos(cdchits[cdc_index]->stereo);
      else
	H(0,state_D)=H_T(state_D,0)=sign;

      // Measurement
      double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
      double one_over_beta=sqrt(1.+mass_hyp*mass_hyp*q_over_p*q_over_p);
      M(0,0)=CDC_DRIFT_SPEED*(cdchits[cdc_index]->t
		 -central_traj[k].s*one_over_beta/SPEED_OF_LIGHT);
      /*
      printf("index %d d %f %f s %f stereo %f\n",cdc_index,M(0,0),cdchits[cdc_index]->d,
	     central_traj[k].s,cdchits[cdc_index]->stereo);
      */
      double cosa=cos(cdchits[cdc_index]->stereo);
      DVector3 mdir(-sin(Sc(state_phi,0))/cosa,
		    cos(Sc(state_phi,0))/cosa,
		    0.);
      DVector3 wirepos=cdchits[cdc_index]->origin
	+(pos.z()-cdchits[cdc_index]->origin.z())*cdchits[cdc_index]->dir;

      double psi=atan2(pos.y(),pos.x());
      H(0,state_D)=H_T(state_D,0)=(-cos(psi)*sin(Sc(state_phi,0))
				   +sin(psi)*cos(Sc(state_phi,0)))/cosa;
      H(0,state_phi)=H_T(state_phi,0)
	=(-(pos.x()-wirepos.x())*cos(Sc(state_phi,0))
	  -(pos.y()-wirepos.y())*sin(Sc(state_phi,0)))/cosa;
      if(cdchits[cdc_index]->stereo!=0.){
	H(0,state_z)=H_T(state_z,0)
	  =(sin(Sc(state_phi,0))*cdchits[cdc_index]->dir.x()
	    -cos(Sc(state_phi,0))*cdchits[cdc_index]->dir.y())/cosa;
      }
      else H(0,state_z)=H_T(state_z,0)=0.;
						

      // Difference and variance
      double dm=(M-H*Sc)(0,0);
      double var=(V+H*(Cc*H_T))(0,0);
      // probability

      printf("dm1 %f %f\n",dm,(H*Sc)(0,0));
      dm=M(0,0)-fabs((pos-wirepos).Dot(mdir));
      double p=exp(-0.5*dm*dm/var);

      p=1.;

      // Inverse of variance
      InvV=p/(p*V+(H*(Cc*H_T))(0,0));
       
      // Compute Kalman gain matrix
      K=InvV*(Cc*H_T);
      /*
      Cc.Print();
      H_T.Print();
      K.Print();
      */
      // Update the state vector 
      //dS=K*(M-H*Sc);

      dS=dm*K;
      /*
      printf("dS\n");
      dS.Print();
      */

      Sc=Sc+dS;

      // update position
      double dz=Sc(state_z,0)-S0(state_z,0);
      pos(0)=central_traj[k].pos.x()+dz*cos(S0(state_phi,0))/S0(state_tanl,0);
      pos(1)=central_traj[k].pos.y()+dz*sin(S0(state_phi,0))/S0(state_tanl,0);
      pos(2)=Sc(state_z,0);
      
      // Path length in active volume
      //path_length+=?
      
      // Update state vector covariance matrix
      Cc=Cc-(K*(H*Cc));  

      /*
      if (DEBUG_HISTS){
	cdc_resid[cdc_index]= M(0,0)-(H*Sc)(0,0);
	if (V-(H*(Cc*H_T))(0,0)>0.){
	  cdc_pulls[cdc_index]= (M(0,0)-(H*Sc)(0,0))/sqrt((V-(H*(Cc*H_T))(0,0)));
	}
	else
	  cdc_pulls[cdc_index]=1000.;
      }
      */
      wirepos=
	cdchits[cdc_index]->origin+(pos.z()
				    -cdchits[cdc_index]->origin.z())*

	cdchits[cdc_index]->dir;
      mdir.SetXYZ(-sin(Sc(state_phi,0)),cos(Sc(state_phi,0)),0.);
      dm=M(0,0)-(pos-wirepos).Dot(mdir);

      // Update chi2 for this hit
      chisq+=dm*dm/(V-(H*(Cc*H_T))(0,0));
      
    }
  }
  
  //Sc.Print();

  // update internal variables
  phi_=Sc(1,0);
  q_over_pt_=Sc(0,0);
  tanl_=Sc(2,0);

  //pos.Print();
  x_=pos.x();
  y_=pos.y();
  z_=Sc(state_z,0);
  
  return NOERROR;
}


// Kalman engine for forward tracks
jerror_t DKalmanFilter::KalmanForward(double mass_hyp, DMatrix &S, DMatrix &C,
				      double &chisq){
  DMatrix M(2,1);  // measurement vector
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

  // Initialize chi squared
  chisq=0;

  S0_=DMatrix(*forward_traj[0].S);
  for (unsigned int k=1;k<forward_traj.size()-1;k++){
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

      // Compute drift distance
      double one_over_beta
	=sqrt(1.+S(state_q_over_p,0)*S(state_q_over_p,0)*mass_hyp*mass_hyp);
      double drift=DRIFT_SPEED*(fdchits[id]->t
				-forward_traj[k].s*one_over_beta
				/SPEED_OF_LIGHT);
      double cosa=fdchits[id]->cosa;
      double sina=fdchits[id]->sina;
      double u=fdchits[id]->uwire;
      double v=fdchits[id]->vstrip;
      drift*=(S(state_x,0)*cosa-S(state_y,0)*sina-u)>0?1.:-1.;

      // The next measurement 
      M(0,0)=u+drift;
      M(1,0)=v;
      
      // ... and its covariance matrix 
      V(0,0)=fdchits[id]->covu;
      V(1,1)=fdchits[id]->covv;

      // To transform from (x,y) to (u,v), need to do a rotation:
      //   u = x*cosa-y*sina
      //   v = y*cosa+x*sina
      H(0,state_x)=H_T(state_x,0)=cosa;
      H(1,state_x)=H_T(state_x,1)=sina;
      H(0,state_y)=H_T(state_y,0)=-sina;
      H(1,state_y)=H_T(state_y,1)=cosa;

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
	_DBG_ << "Kalman filter:  Singular matrix..." << endl;
	return UNRECOVERABLE_ERROR;
      }
      
      // Compute Kalman gain matrix
      K=C*(H_T*InvV);

      // Update the state vector 
      S=S+K*(M-H*S); 
      
      // Path length in active volume
      path_length+=1.0*sqrt(1.+S(state_tx,0)*S(state_tx,0)
			    +S(state_ty,0)*S(state_ty,0));
      
      // Update state vector covariance matrix
      C=C-K*(H*C);    
      
      // Residuals
      R=M-H*S;
      R_T=DMatrix(DMatrix::kTransposed,R);
      RC=V-H*(C*H_T);

      if (DEBUG_HISTS && last_iter){
	TH2F *fdc_xresiduals=(TH2F*)gROOT->FindObject("fdc_xresiduals");
	if (fdc_xresiduals) fdc_xresiduals->Fill(fdchits[id]->z,R(state_x,0));
	
	TH2F *fdc_yresiduals=(TH2F*)gROOT->FindObject("fdc_yresiduals");
	if (fdc_yresiduals) fdc_yresiduals->Fill(fdchits[id]->z,R(state_y,0));  
	
	TH2F *fdc_ypulls=(TH2F*)gROOT->FindObject("fdc_ypulls");
	if (fdc_ypulls) fdc_ypulls->Fill(fdchits[id]->z,R(state_y,0)/sqrt(RC(1,1)));
	
	TH2F *fdc_xpulls=(TH2F*)gROOT->FindObject("fdc_xpulls");
	if (fdc_xpulls) fdc_xpulls->Fill(fdchits[id]->z,R(state_y,0)/sqrt(RC(0,0)));
      }
      
      // Calculate the inverse of RC
      det=RC(0,0)*RC(1,1)-RC(0,1)*RC(1,0);
      if (det!=0){
	InvRC(0,0)=RC(1,1)/det;
	InvRC(1,0)=-RC(1,0)/det;
	InvRC(0,1)=-RC(0,1)/det;
	InvRC(1,1)=RC(0,0)/det;
      }
      else{
	_DBG_ << "Kalman filter:  Singular matrix RC..." << endl;
	return UNRECOVERABLE_ERROR;
      }
      
      // Update chi2 for this segment
      chisq+=(R_T*(InvRC*R))(0,0);
    }

  }

  // Final position for this leg
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=forward_traj[forward_traj.size()-1].pos.Z();
  
  return NOERROR;
}
// Kalman engine for forward tracks -- this routine adds CDC hits
jerror_t DKalmanFilter::KalmanForwardCDC(double mass_hyp,double anneal,
					 DMatrix &S, 
					 DMatrix &C,double &chisq){
  DMatrix M(1,1);  // measurement vector
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

  S0_=DMatrix(*forward_traj_cdc[0].S);
  for (unsigned int k=1;k<forward_traj_cdc.size()-1;k++){
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

    // Save the current state of the reference trajectory
    S0_=S0;

    // Add the hit
    if (forward_traj_cdc[k].h_id>0){
      unsigned int id=forward_traj_cdc[k].h_id-1;
      // z position
      double z=forward_traj_cdc[k].pos.z();

      // The next measurement
      double one_over_beta
	=sqrt(1.+mass_hyp*mass_hyp*S(state_q_over_p,0)*S(state_q_over_p,0));
      double dm=CDC_DRIFT_SPEED*(cdchits[id]->t
				 -forward_traj_cdc[k].s*one_over_beta
				 /SPEED_OF_LIGHT);
      // wire position
      DVector3 wirepos=cdchits[id]->origin
	+(z-cdchits[id]->origin.z())*cdchits[id]->dir;
      double xw=wirepos.x();
      double yw=wirepos.y();

      // position on track
      double x=S(state_x,0),y=S(state_y,0);
      
      // direction variables
      double tx=S(state_tx,0),ty=S(state_ty,0);
      double cosa=cos(cdchits[id]->stereo);
      double fac=1./sqrt(tx*tx+ty*ty);

      // predicted doca
      double d=(-(x-xw)*ty+(y-yw)*tx)*fac/cosa;
      double sign=(d>0)?1.:-1.;
      dm*=sign;
    
      // Track projection
      H(0,state_x)=H_T(state_x,0)=-ty*fac/cosa;
      H(0,state_y)=H_T(state_y,0)=tx*fac/cosa;
      H(0,state_tx)=H_T(state_tx,0)
	=((x-xw)*tx*ty+(y-yw)*ty*ty)*fac*fac*fac/cosa;
      H(0,state_ty)=H_T(state_ty,0)
	=(-(x-xw)*tx*tx-(y-yw)*tx*ty)*fac*fac*fac/cosa;

      // Inverse of covariance matrix // In verse of variance
      InvV=1./(V+(H*(C*H_T))(0,0));

      // Compute Kalman gain matrix
      K=InvV*(C*H_T);
      
      // Update the state vector 
      S=S+(dm-d)*K;
      
      // Path length in active volume
      // path_length+=?
      
      // Update state vector covariance matrix
      C=C-K*(H*C);    
      
      // doca after correction
      x=S(state_x,0);
      y=S(state_y,0);
      tx=S(state_tx,0);
      ty=S(state_ty,0);
      fac=1./sqrt(tx*tx+ty*ty);
      d=(-(x-xw)*ty+(y-yw)*tx)*fac/cosa;
     
      // Residual
      double res=dm-d;

      // Update chi2 for this segment
      chisq+=res*res/(V-(H*(C*H_T))(0,0));
    }

  }

  // Final position for this leg
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=forward_traj_cdc[forward_traj_cdc.size()-1].pos.Z();
  
  return NOERROR;
}

// Extrapolate to the point along z of closest approach to the beam line using 
// the forward track state vector parameterization.  Converts to the central
// track representation at the end.
jerror_t DKalmanFilter::ExtrapolateToVertex(double mass_hyp,DMatrix S, 
					    DMatrix C){
  DMatrix J(5,5);  //.Jacobian matrix
  DMatrix JT(5,5); // and its transpose
  DMatrix Q(5,5); // multiple scattering matrix
  DMatrix Sc(5,1),Cc(5,5);  // central representation

  // position variables
  double ds=0.;
  double z=z_,newz;
  double dz=-STEP_SIZE;
  double r2_old=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);

  double dEdx=0.;

  // Check the direction of propagation
  DMatrix S0(5,1);
  S0=S;
  Step(z,z+dz,dEdx,S0);
  double r2=S0(state_x,0)*S0(state_x,0)+S0(state_y,0)*S0(state_y,0);
  if (r2>r2_old) dz*=-1.;
  
  while (z>Z_MIN && z<Z_MAX){
    newz=z+dz;
    // Step through field
    ds=Step(z,newz,dEdx,S);
    r2=S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0);
    if (r2>r2_old){
      DVector3 origin(0,0,65.);
      DVector3 dir(0,0,1.);

      // Find position of best doca to beam line
      GoldenSection(newz,dz,dEdx,origin,dir,S);
      
      break;
    }
    r2_old=r2;
    z=newz;
  }

  // Convert from forward rep. to central rep.
  ConvertStateVector(newz,0.,0.,S,C,Sc,Cc);

  // Track Parameters at "vertex"
  phi_=Sc(state_phi,0);
  q_over_pt_=Sc(state_q_over_pt,0);
  tanl_=Sc(state_tanl,0);
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=newz;

  return NOERROR;
}


// Propagate track to point of distance of closest approach to origin
jerror_t DKalmanFilter::ExtrapolateToVertex(double mass_hyp,DVector3 pos,
					    DMatrix Sc,DMatrix Cc){
  DMatrix Jc(5,5);  //.Jacobian matrix
  DMatrix JcT(5,5); // and its transpose
  DMatrix Q(5,5); // multiple scattering matrix

  // Initialize the beam position = center of target, and the direction
  DVector3 origin(0,0,65.);  
  DVector3 dir(0,0,1.);

  // Position and step variables
  double r=pos.Perp();
  double ds=-0.25; // step along path in cm
  double r_old=r;
  Sc(state_D,0)=r;

  // Energy loss
  double dedx=0.;

  // Check direction of propagation
  DMatrix S0(5,1);
  S0=Sc; 
  DVector3 pos0=pos;
  Step(pos0,origin,dir,ds,S0,dedx);
  r=pos0.Perp();
  if (r>r_old) ds*=-1.;
 
  // Track propagation loop
  while (Sc(state_z,0)>0. && r>BEAM_RADIUS && r<65.){ 
    dedx=0.;

    // Compute the Jacobian matrix
    StepJacobian(pos,origin,dir,ds,Sc,dedx,Jc);

    // Process noise
    double q_over_p=Sc(state_q_over_pt,0)*cos(atan(Sc(state_tanl,0)));
    GetProcessNoiseCentral(0.14,ds,30420.,Sc,Q); // use Air for now

    // Propagate the covariance matrix
    JcT=DMatrix(DMatrix::kTransposed,Jc);
    Cc=Jc*(Cc*JcT)+Q;
   
    // Propagate the state through the field
    Step(pos,origin,dir,ds,Sc,dedx);

    r=pos.Perp();
    if (r>r_old) {
   
    
      break;
    }
    r_old=r;
  }

  // Track Parameters at "vertex"
  phi_=Sc(1,0);
  q_over_pt_=Sc(0,0);
  tanl_=Sc(2,0);
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();

  printf("p %f Best vertex \n",1./q_over_pt_/cos(atan(tanl_)));
  pos.Print();

  return NOERROR;
}


// ------------------------ UNUSED CODE -------------------------------------
#if 0
// Kalman engine for Central tracks; updates the position on the trajectory
// after the last hit (closest to the target) is added
jerror_t DKalmanFilter::KalmanCentral(double mass_hyp,double anneal_factor,
				      DMatrix &Sc,DMatrix &Cc,
				      DVector3 &pos,double &chisq){
  DMatrix H(1,5);  // Track projection matrix
  DMatrix H_T(5,1); // Transpose of track projection matrix
  DMatrix Jc(5,5);  // State vector Jacobian matrix
  DMatrix Jc_T(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,1);  // Kalman gain matrix
  double V=0.8*0.8/12.;  // Measurement variance
  // double V=0.000225;
  double InvV; // inverse of variance
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix S0(5,1); // state vector
  DMatrix M(1,1);
  
  // Initialize the chi2 for this part of the track
  chisq=0.;

  //Track projection matrix
  H(0,state_D)=H_T(state_D,0)=1.;
  
  // path length increment, dedx, and position variables.
  double ds=-0.5,ds2=-0.1;
  double dedx=0.;
  double R,r,x,y,dx,dy;
  double q=Sc(state_q_over_pt,0)/fabs(Sc(state_q_over_pt,0));

  // vector of state information
  vector<DKalmanCentral_t>cdc_state;
  DKalmanCentral_t temp;
  
  // First position 
  temp.pos=pos;
  temp.q_over_pt=Sc(state_q_over_pt,0);
  temp.phi=Sc(state_phi,0);
  temp.tanl=Sc(state_tanl,0);
  temp.D=Sc(state_D,0);
  temp.z=Sc(state_z,0);
  cdc_state.push_back(temp);

  // Loop over CDC hits
  DVector3 wirepos; 
  double doca=0.;
  for (unsigned int m=0;m<cdchits.size();m++){
    // Reset ds
    ds=-0.2;

    // Find the doca at the current position to set Sc(state_D,0)
    x=cdchits[m]->origin.x()
      +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.x(); 
    y=cdchits[m]->origin.y()
      +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.y();
    dx=pos.x()-x;
    dy=pos.y()-y;
    
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
    double sign=q;
    if (rc2<Rc*Rc) sign*=-1.;
      
    // doca
    Sc(state_D,0)=sign*sqrt(dx*dx+dy*dy);   

    // Current positions of wire and track
    R=sqrt(x*x+y*y);
    r=pos.Perp();
    doca=Sc(state_D,0);

     // Check that we are not already at the minimum doca
    S0=Sc;
    DVector3 pos0=pos;
    Step(pos0,cdchits[m]->origin,cdchits[m]->dir,ds,S0,0.);
    if (fabs(Sc(state_D,0))>fabs(S0(state_D,0))){
      while (r>targ_wall[1] && fabs(doca)>EPS){
	// Check that z value makes sense 
	if (pos.z()<0 || pos.z()>400){
	  chisq=1e8;
	  _DBG_<<" -- z value out of range!! -- " <<endl;
	  return VALUE_OUT_OF_RANGE;
	}
     
	// some track parameters for this step
	double tanl=Sc(state_tanl,0);
	double cosl=cos(atan(tanl));
	double q_over_p=Sc(state_q_over_pt,0)*cosl;
	
	// Energy loss and process noise (multiple scattering)
	if (fabs(doca)<0.8){
	  ds=-0.257*doca+0.0056;
	  if (fabs(doca)<0.1) ds=-0.02;
	  GetProcessNoiseCentral(mass_hyp,ds,10029.,Sc,Q); // use Ar for now
	  dedx=GetdEdx(0.14,q_over_p,18.,39.9,0.00166);
	}
	else{
	  GetProcessNoiseCentral(mass_hyp,ds,30420.,Sc,Q); // use Air for now
	  dedx=GetdEdx(0.14,q_over_p,7.,14.,1.205e-3);
	}
	dedx=0.;

	// Propagate the state and the covariance matrix through the field
	StepJacobian(pos,cdchits[m]->origin,cdchits[m]->dir,ds,Sc,dedx,Jc);
	
	// Transpose of Jacobian matrix
	Jc_T=DMatrix(DMatrix::kTransposed,Jc);
	
	// Propagate the covariance matrix
	Cc=Jc*(Cc*Jc_T)+Q;
	
	// Add the current state to the vector
	temp.pos=pos;
	temp.q_over_pt=Sc(state_q_over_pt,0);
	temp.phi=Sc(state_phi,0);
	temp.tanl=Sc(state_tanl,0);
	temp.D=Sc(state_D,0);
	temp.z=Sc(state_z,0);
	cdc_state.push_back(temp);
	
	if ((fabs(doca)<fabs(Sc(state_D,0)))
	    || (doca/fabs(doca)!=Sc(state_D,0)/fabs(Sc(state_D,0)))
	    ){
	  // We've bracketed the minimum; now use the golden section algorithm
	  // to find the best doca.  See Numerical Recipes in C, pp 401-402.	
	  ds=GoldenSection(ds,ds,doca,dedx,pos,cdchits[m]->origin,cdchits[m]->dir,Sc);
	  
	  // Back up to the position just before the doca
	  unsigned int c_id=cdc_state.size()-3;
	  Sc(state_q_over_pt,0)=cdc_state[c_id].q_over_pt;
	  Sc(state_phi,0)=cdc_state[c_id].phi;
	  Sc(state_tanl,0)=cdc_state[c_id].tanl;
	  Sc(state_D,0)=cdc_state[c_id].D;
	  Sc(state_z,0)=cdc_state[c_id].z;
	  pos=cdc_state[c_id].pos;
	  	  
	  // Energy loss and process noise (multiple scattering)
	  if (fabs(doca)<0.8){
	    GetProcessNoiseCentral(mass_hyp,ds,10029.,Sc,Q); // use Ar for now
	    dedx=GetdEdx(0.14,q_over_p,18.,39.9,0.00166);
	  }
	  else{
	    GetProcessNoiseCentral(mass_hyp,ds,30420.,Sc,Q); // use Air for now
	    dedx=GetdEdx(0.14,q_over_p,7.,14.,1.205e-3);
	  }
	  dedx=0.;

	  // Step to the doca 
	  StepJacobian(pos,cdchits[m]->origin,cdchits[m]->dir,ds,Sc,dedx,Jc);

	  // Transpose of Jacobian matrix
	  Jc_T=DMatrix(DMatrix::kTransposed,Jc);
	  
	  // Propagate the covariance matrix
	  Cc=Jc*(Cc*Jc_T)+Q;
	  
	// Update the appropriate entry in the list with the new doca info
	  cdc_state[c_id+1].pos=pos;
	  cdc_state[c_id+1].q_over_pt=Sc(state_q_over_pt,0);
	  cdc_state[c_id+1].phi=Sc(state_phi,0);
	  cdc_state[c_id+1].tanl=Sc(state_tanl,0);
	  cdc_state[c_id+1].D=Sc(state_D,0);
	  cdc_state[c_id+1].z=Sc(state_z,0);
	  
	  // Remove the point on the trajectory after the doca
	  cdc_state.pop_back();
	  
	  break;
	}
	
	// Find the radius at the end of this step
	x=cdchits[m]->origin.x()
	  +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.x(); 
	y=cdchits[m]->origin.y()
	  +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.y();
	R=sqrt(x*x+y*y);
	r=pos.Perp();
	
	// Save the doca
	doca=Sc(state_D,0);
      }    
    }
    // The next part of the code applies a correction to the state vector using the hit information

    H(0,state_D)=H_T(state_D,0)=Sc(state_D,0)/fabs(Sc(state_D,0));

    // Inverse of variance
    InvV=1./(V+(H*(Cc*H_T))(0,0));

    double chi2_temp=(M(0,0)-(H*Sc)(0,0))*(M(0,0)-(H*Sc)(0,0))*InvV;
    double chi2_cut = 9.0;
       
    // Compute Kalman gain matrix
    K=InvV*(Cc*H_T);

    // Update the state vector 
    M(0,0)=cdchits[m]->d;
    M(0,0)=0.;
    
    //if (sqrt(InvV)*fabs(M(0,0)-(H*Sc)(0,0))<2.5)
       {  
      dS=K*(M-H*Sc);
      
      Sc=Sc+dS;
           
      // Update the position
      double lambda=atan(Sc(state_tanl,0));
      double sinl=sin(lambda);
      double cosl=cos(lambda);
      sinphi=sin(Sc(state_phi,0));
      cosphi=cos(Sc(state_phi,0));
      ds=(Sc(state_z,0)-pos.Z())/sinl;
      pos(0)+=cosl*cosphi*ds;
      pos(1)+=cosl*sinphi*ds;
      pos(2)=Sc(state_z,0);
      
      // Path length in active volume
      //path_length+=?
      
      // Update state vector covariance matrix
      Cc=Cc-(K*(H*Cc));  
    
      if (DEBUG_HISTS){
	cdc_resid[m]= M(0,0)-(H*Sc)(0,0);
	if (V-(H*(Cc*H_T))(0,0)>0.){
	  cdc_pulls[m]= (M(0,0)-(H*Sc)(0,0))/sqrt((V-(H*(Cc*H_T))(0,0)));
	}
	else
	  cdc_pulls[m]=1000.;
      }
      
      // Update chi2 for this hit
      chisq+=(M(0,0)-(H*Sc)(0,0))*(M(0,0)-(H*Sc)(0,0))/(V-(H*(Cc*H_T))(0,0));
    } 
    //printf("p %f theta %f\n",1./fabs(Sc(state_q_over_pt,0))/cos(atan(Sc(state_tanl,0))),M_PI/2.-atan(Sc(state_tanl,0)));

    // Save the doca
    if (m<cdchits.size()-1){
      double x=cdchits[m+1]->origin.x()
	+(pos.z()-cdchits[m+1]->origin.z())*cdchits[m+1]->dir.x(); 
      double y=cdchits[m+1]->origin.y()
	+(pos.z()-cdchits[m+1]->origin.z())*cdchits[m+1]->dir.y();
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
      double sign=q;
      if (rc2<Rc*Rc) sign*=-1.;
      
      // New doca
      Sc(state_D,0)=sign*sqrt(dx*dx+dy*dy); //D for current state vector

      doca=Sc(state_D,0);
    }
    wirepos=cdchits[m]->origin +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir; 
  }

  // update internal variables
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();
  
  return NOERROR;
}
#endif
