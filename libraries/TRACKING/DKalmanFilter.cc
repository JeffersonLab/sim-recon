//************************************************************************
// DKalmanFilter.cc
//************************************************************************

#include "DKalmanFilter.h"
#include "CDC/DCDCTrackHit.h"
#include "DANA/DApplication.h"
#include "HDGEOMETRY/DMagneticFieldMap.h"
#include "HDGEOMETRY/DGeometry.h"
#include <TDecompLU.h>

#include <math.h>

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define EPS 1.0e-8
#define DEDX_ENDPLATE (3.95e-3) // Carbon
#define DEDX_AIR  (2.19e-6) // GeV/cm
#define DEDX_LH2 ( 2.856e-4) 
#define DEDX_ROHACELL (0.4e-3) // low density carbon for now!!
#define BEAM_RADIUS  0.1 
#define NUM_ITER 10

// Local boolean routines for sorting
bool DKalmanHit_cmp(DKalmanHit_t *a, DKalmanHit_t *b){
  return a->z>b->z;
}
bool DKalmanCDCHit_cmp(DKalmanCDCHit_t *a, DKalmanCDCHit_t *b){
  double ra=a->origin.Perp();
  double rb=b->origin.Perp();
  return (rb<ra);
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
  hit->dir=cdchit->wire->udir;

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
double DKalmanFilter::StepJacobian(double oldz,double newz,DMatrix &S,
				     double dEdx,DMatrix &J){
   // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix J1(5,5),J2(5,5),J3(5,5),J4(5,5);  
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  
  double delta_z=newz-oldz;
  double s=sqrt(1.+S(state_tx,0)*S(state_tx,0)+S(state_ty,0)*S(state_ty,0))
    *delta_z;
  CalcDerivAndJacobian(oldz,S,dEdx,J1,D1);
  CalcDerivAndJacobian(oldz+delta_z/2.,S+0.5*delta_z*D1,dEdx,J2,D2);
  J2=J2+0.5*J2*J1;
  CalcDerivAndJacobian(oldz+delta_z/2.,S+0.5*delta_z*D2,dEdx,J3,D3);
  J3=J3+0.5*J3*J2;
  CalcDerivAndJacobian(oldz+delta_z,S+delta_z*D3,dEdx,J4,D4);
  J4=J4+J4*J3;

  S+=delta_z*((1./6.)*D1+(1./3.)*D2+(1./3.)*D3+(1./6.)*D4);
  J+=delta_z*((1./6.)*J1+(1./3.)*J2+(1./3.)*J3+(1./6.)*J4);
  
  return s;
}

// Calculate the derivative and Jacobian matrices for the alternate set of 
// parameters {q/pT, phi, tan(lambda),D,z}
jerror_t DKalmanFilter::CalcDerivAndJacobian(DVector3 &pos,
					     DVector3 wire_pos,double ds,
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
  double D=S(state_D,0);
  double z=S(state_z,0);

  //B-field and field gradient at (x,y,z)
  double Bx,By,Bz,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
  bfield->GetFieldGradient(pos.x(),pos.y(),pos.z(),dBxdx,dBxdy,dBxdz,dBydx,
			   dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);
  double R=1./qBr2p/fabs(q_over_pt)/B;
  double xc=pos.x()-R*sinphi;
  double yc=pos.y()+R*cosphi;
  double rc=sqrt((wire_pos.x()-xc)*(wire_pos.x()-xc)
		 +(wire_pos.y()-yc)*(wire_pos.y()-yc));
  double a=D*D-(R-rc)*(R-rc);
  if (a>0.) a=sqrt(R*rc*a);
  else a=0.;

  // Derivative of S with respect to s
  D1(state_q_over_pt,0)=qBr2p*q_over_pt*q_over_pt*sinl*(By*cosphi-Bx*sinphi);
  D1(state_phi,0)=qBr2p*q_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_tanl,0)=qBr2p*q_over_pt*(By*cosphi-Bx*sinphi)/cosl;
  D1(state_D,0)=-D1(state_phi,0)/D*(rc*ds*cosl/R-a);
  D1(state_z,0)=sinl;

  // Jacobian matrix elements
  J1(state_phi,state_phi)=qBr2p*q_over_pt*sinl*(By*cosphi-Bx*sinphi);
  J1(state_phi,state_q_over_pt)=qBr2p*(Bx*cosphi*sinl+By*sinphi*sinl
				       -Bz*cosl);
  J1(state_phi,state_tanl)=qBr2p*q_over_pt*(Bx*cosphi*cosl+By*sinphi*cosl
					    -Bz*sinl)/(1.+tanl*tanl);
  J1(state_tanl,state_phi)
    =-qBr2p*q_over_pt*cosl*(1.+tanl*tanl)*(By*sinphi+Bx*cosphi);
  J1(state_tanl,state_q_over_pt)=
    qBr2p*cosl*(1.+tanl*tanl)*(By*cosphi-Bx*sinphi);
  J1(state_tanl,state_tanl)=qBr2p*q_over_pt*sinl*(By*cosphi-Bx*sinphi);
  
  J1(state_q_over_pt,state_phi)
    =-qBr2p*q_over_pt*q_over_pt*sinl*(By*sinphi+Bx*cosphi);
  J1(state_q_over_pt,state_q_over_pt)
    =2.*qBr2p*q_over_pt*sinl*(By*cosphi-Bx*sinphi);
  J1(state_q_over_pt,state_tanl)
    =qBr2p*q_over_pt*q_over_pt*cosl*cosl*cosl*(By*cosphi-Bx*sinphi);

  J1(state_z,state_tanl)=cosl/(1.+tanl*tanl);

  // New position 
  pos.SetXYZ(pos.x()+cosl*cosphi*ds,pos.y()+cosl*sinphi*ds,pos.z()+sinl*ds);

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
  Sc(state_D,0)=sqrt((x-wire_x)*(x-wire_x)+(y-wire_y)*(y-wire_y));
  Sc(state_z,0)=z;

  DMatrix J(5,5);
  J(state_tanl,state_tx)=-tx*tanl*tanl*tanl;
  J(state_tanl,state_ty)=-ty*tanl*tanl*tanl;
  J(state_z,state_x)=1./tx;
  J(state_z,state_y)=1./ty;
  J(state_q_over_pt,state_q_over_p)=1./cosl;
  J(state_q_over_pt,state_tx)=-tx*q_over_p*tanl*tanl*tanl*factor;
  J(state_q_over_pt,state_ty)=-ty*q_over_p*tanl*tanl*tanl*factor;
  J(state_phi,state_tx)=-ty*factor*factor;
  J(state_phi,state_ty)=tx*factor*factor;
  J(state_D,state_x)=(x-wire_x)/S(state_D,0);
  J(state_D,state_y)=(y-wire_y)/S(state_D,0);

  Cc=J*(C*DMatrix(DMatrix::kTransposed,J));

  return NOERROR;
}

// Runga-Kutte for alternate parameter set {q/pT,phi,tanl(lambda),D,z}
jerror_t DKalmanFilter::StepJacobian(DVector3 &pos,DVector3 wire_pos,
				     double ds,DMatrix &S,
				     double dEdx,DMatrix &J){
  // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix J1(5,5),J2(5,5),J3(5,5),J4(5,5);  
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  
  CalcDerivAndJacobian(pos,wire_pos,ds/2.,S,dEdx,J1,D1);
  CalcDerivAndJacobian(pos,wire_pos,0.,S+0.5*ds*D1,dEdx,J2,D2);
  J2=J2+0.5*J2*J1;
  CalcDerivAndJacobian(pos,wire_pos,ds/2.,S+0.5*ds*D2,dEdx,J3,D3);
  J3=J3+0.5*J3*J2;
  CalcDerivAndJacobian(pos,wire_pos,0.,S+ds*D3,dEdx,J4,D4);
  J4=J4+J4*J3;

  // New state vector and Jacobian matrix
  S+=ds*((1./6.)*D1+(1./3.)*D2+(1./3.)*D3+(1./6.)*D4);
  J+=ds*((1./6.)*J1+(1./3.)*J2+(1./3.)*J3+(1./6.)*J4);
  
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

  if (Z>12) I0=(9.76*Z+58.8*pow(Z,-0.19))*1e-3;
  return -0.0001535*Z/A*rho/beta2*(log(2.*Me*betagamma*betagamma*Tmax/I0/I0)
			       -2.*beta2);
}

// Interface routine for Kalman filter
jerror_t DKalmanFilter::KalmanLoop(double mass_hyp){
  DMatrix S(5,1),C(5,5),Sc(5,1),Cc(5,5);

  if (cdchits.size()>0){
    // Order the CDC hits by radius
    sort(cdchits.begin(),cdchits.end(),DKalmanCDCHit_cmp);

    Sc(state_q_over_pt,0)=q_over_pt_;
    Sc(state_phi,0)=phi_;
    Sc(state_tanl,0)=tanl_;
    Sc(state_z,0)=z_;

  }
  
  // Initialize the covariance matrix
  for (unsigned int i=0;i<5;i++) C(i,i)=Cc(i,i)=1;

  // Fit in Forward region
  if (hits.size()>0){   
    // Initialize the state vector 
    S(state_x,0)=x_;
    S(state_y,0)=y_;
    S(state_tx,0)=tx_;
    S(state_ty,0)=ty_;
    S(state_q_over_p,0)=q_over_p_; 
    
    KalmanForward(mass_hyp,S,C); 

    // Convert to Central Track parameter set
    double wire_x=0,wire_y=0;
    if (cdchits.size()>0){
      wire_x=cdchits[0]->origin.x()
	+(z_-cdchits[0]->origin.z())*cdchits[0]->dir.x(); 
      wire_y=cdchits[0]->origin.y()
	+(z_-cdchits[0]->origin.z())*cdchits[0]->dir.y();
    }
    ConvertStateVector(z_,wire_x,wire_y,S,C,Sc,Cc);
  }

  // Fit in Central region
  if (cdchits.size()>0){
    KalmanCentral(mass_hyp,Sc,Cc);
  }
  
  // Find track parameters where track crosses beam line
  ExtrapolateToVertex(mass_hyp,Sc,Cc);

  return NOERROR;
}

// Kalman engine for Central tracks
jerror_t DKalmanFilter::KalmanCentral(double mass_hyp,DMatrix &Sc,DMatrix &Cc){
  DVector3 pos(x_,y_,z_);
  DMatrix H(1,5);  // Track projection matrix
  DMatrix H_T(5,1); // Transpose of track projection matrix
  DMatrix Jc(5,5);  // State vector Jacobian matrix
  DMatrix Jc_T(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,1);  // Kalman gain matrix
  double V=0.15*0.15;  // Measurement variance
  double InvV; // inverse of variance
  DMatrix S0(5,1); //State vector
  DMatrix dS(5,1);  // perturbation in state vector

  // Track projection matrix
  H(0,3)=1.;
  H_T=DMatrix(DMatrix::kTransposed,H);
    
  // Set the seed
  S0=Sc;
  // Find the doca at the first position to set S0(state_D,0)
  double x=cdchits[0]->origin.x()
    +(pos.z()-cdchits[0]->origin.z())*cdchits[0]->dir.x(); 
  double y=cdchits[0]->origin.y()
    +(pos.z()-cdchits[0]->origin.z())*cdchits[0]->dir.y();
  double dx=pos.x()-x;
  double dy=pos.y()-y;
  S0(state_D,0)=Sc(state_D,0)=sqrt(dx*dx+dy*dy);
 
  // Loop over CDC hits
  double ds=-0.5;
  double dedx=0.;
  double R,r;
  DVector3 wire_pos;
 
  //for (int iter=0;iter<NUM_ITER;iter++)
    {
    double chisq=0;
    unsigned int index=0;
    
    central_traj.clear();

    for (unsigned int m=0;m<cdchits.size();m++){
      double doca=0.,doca_min=1000.;
      // Find the doca at the current position to set S0(state_D,0)
      x=cdchits[m]->origin.x()
	+(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.x(); 
      y=cdchits[m]->origin.y()
	+(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.y();
      dx=pos.x()-x;
      dy=pos.y()-y;
      S0(state_D,0)=Sc(state_D,0)=sqrt(dx*dx+dy*dy);   
      
      ds=-0.3;    
      do{
	dS=Sc-S0;
	
	// Propagate the state and the covariance matrix through the field
	wire_pos.SetXYZ(x,y,0.);
	StepJacobian(pos,wire_pos,ds,S0,dedx,Jc);
	Jc_T=DMatrix(DMatrix::kTransposed,Jc);
	Cc=Jc*(Cc*Jc_T);
	
	// Update the "improved" state vector
	Sc=S0+Jc*dS;
	
	DKalmanState_t temp;
	for (unsigned int i=0;i<5;i++){
	  temp.S[i]=Sc(i,0);
	    for (unsigned int j=0;j<5;j++){
	      temp.C[i][j]=Cc(i,j);
	      temp.JT[i][j]=Jc_T(i,j);
	    }
	}
	temp.measurement=false;
	central_traj.push_front(temp);
   	index++;

	// Find the radius at the end of this step
	x=cdchits[m]->origin.x()
	  +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.x(); 
	y=cdchits[m]->origin.y()
	  +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.y();
	R=sqrt(x*x+y*y);
	r=pos.Perp(); 
      }while(r>R+0.8 && r<65.);
      
      ds=-0.1;
      while (r>R-1.6){
	doca=Sc(state_D,0);
	if (doca>doca_min){
	  break;
	}
	doca_min=doca;
	
	// Perturbation to state vector
	dS=Sc-S0;      
	
	// Propagate the state and the covariance matrix through the field
	wire_pos.SetXYZ(x,y,0.);
	StepJacobian(pos,wire_pos,ds,S0,dedx,Jc);
	Jc_T=DMatrix(DMatrix::kTransposed,Jc);
	Cc=Jc*(Cc*Jc_T);
	
	// Update the "improved" state vector
	Sc=S0+Jc*dS;
	
	DKalmanState_t temp;
	for (unsigned int i=0;i<5;i++){
	  temp.S[i]=Sc(i,0);
	  for (unsigned int j=0;j<5;j++){
	    temp.C[i][j]=Cc(i,j);
	    temp.JT[i][j]=Jc_T(i,j);
	  }
	}
	temp.measurement=false;
	central_traj.push_front(temp);
	index++;

	// Find the radius at the end of this step
	x=cdchits[m]->origin.x()
	  +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.x(); 
	y=cdchits[m]->origin.y()
	  +(pos.z()-cdchits[m]->origin.z())*cdchits[m]->dir.y();
	R=sqrt(x*x+y*y);
	r=pos.Perp(); 
      }
      // Inverse of variance
      InvV=1./(V+(H*(Cc*H_T))(0,0));
      
      // Compute Kalman gain matrix
      K=InvV*(Cc*H_T);
      
      // Update the state vector 
      Sc=Sc+K*(cdchits[m]->d-doca);
            
      // Path length in active volume
      //path_length+=?
      
      // Update state vector covariance matrix
      Cc=Cc-K*(H*Cc);  

      // Update the list of state vectors and covariances
      DKalmanState_t temp;
      for (unsigned int i=0;i<5;i++){
	temp.S[i]=Sc(i,0);
	for (unsigned int j=0;j<5;j++){
	  temp.C[i][j]=Cc(i,j);
	  temp.JT[i][j]=Jc_T(i,j);
	}
      } 
      temp.measurement=true;
      central_traj.push_front(temp);
      index++;
      
      // Update chi2 for this hit
      chisq+=(cdchits[m]->d-Sc(state_D,0))*(cdchits[m]->d-Sc(state_D,0))
	/(V-(H*(Cc*H_T))(0,0));
    }
    
    /*
      DMatrix Ck(5,5),Jk(5,5),Sk(5,1),Ak(5,5),AkT(5,5);
      DMatrix Ck_old(5,5),Sk_old(5,5);
    for (unsigned int i=1; i<central_traj.size();i++){
    for (unsigned int j=0;j<5;j++){
    Sk(j,0)=central_traj[i].S[j];
    for (unsigned int k=0;k<5;k++){
    Jk(j,k)=central_traj[i].JT[j][k];
    Ck(j,k)=central_traj[i].C[j][k];
	}
	}
	printf("Measurement? %d\n",central_traj[i].measurement);
	Sk.Print();
	Ck.Print();
      if (central_traj[i].measurement==false){
      Ak=Ck*(Jk*DMatrix(DMatrix::kInverted,Ck));
      AkT=DMatrix(DMatrix::kTransposed,Ak);
      Sc=Sk+Ak*(Sc-Sk);
      Cc=Ck+Ak*((Cc-Ck)*AkT);
	Sk_old=Sk;
	Ck_old=Ck;
	}
	else{
	Ak=Ck*(Jk*DMatrix(DMatrix::kInverted,Ck_old));
	AkT=DMatrix(DMatrix::kTransposed,Ak);
	Sc=Sk+Ak*(Sc-Sk_old);
	Cc=Ck+Ak*((Cc-Ck_old)*AkT);
	}
	printf("Smoothed\n");
	Sc.Print();
      Cc.Print();
      
      }
    */
    }
  x_=pos.x();
  y_=pos.y();
  z_=pos.z();
  

  return NOERROR;
}


// Kalman engine for forward tracks
jerror_t DKalmanFilter::KalmanForward(double mass_hyp, DMatrix &S, DMatrix &C){
  DMatrix M(2,1);  // measurement vector
  DMatrix M_pred(2,1); // prediction for hit position
  DMatrix H(2,5);  // Track projection matrix
  DMatrix H_T(5,2); // Transpose of track projection matrix
  DMatrix J(5,5);  // State vector Jacobian matrix
  DMatrix J_T(5,5); // transpose of this matrix
  DMatrix Q(5,5);  // Process noise covariance matrix
  DMatrix K(5,2);  // Kalman gain matrix
  DMatrix V(2,2);  // Measurement covariance matrix
  DMatrix V1(2,2);  // same, with the contribution from C
  DMatrix X(2,1);  // Position on track
  DMatrix R(2,1);  // Filtered residual
  DMatrix R_T(1,2);  // ...and its transpose
  DMatrix RC(2,2);  // Covariance of filtered residual
  DMatrix InvRC(2,2); // and its inverse
  DMatrix S0(5,1); //State vector
  DMatrix SBest(5,1);
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix CBest(5,5);   // Covariance matrix for state vector
  DMatrix InvV(2,2); // Inverse of error matrix

  // path length increment
  double ds=0;
  // Path length in active volume
  path_length=0;
  // Energy loss
  double dedx=0.;

  //Set the seed
  S0=S;
  
  // Initialize the chisq
  chisq_=1e8;

  // Order the hits from the most downstream to the most upstream
  sort(hits.begin(),hits.end(),DKalmanHit_cmp);

  // Track projection matrix
  H(0,0)=H(1,1)=1.;
  H_T=DMatrix(DMatrix::kTransposed,H);

  // Loop over hits, updating the state vector at each step 
  double endz=hits[0]->z,newz,oldz; // we increment in z
  z_=endz;
  for (int iter=0;iter<NUM_ITER;iter++){
    double chisq=0.;

    for (unsigned int k=0;k<hits.size();k++){      
      endz=hits[k]->z;
    
      // The next measurement 
      M(0,0)=hits[k]->x;
      M(1,0)=hits[k]->y; 
      
      // ... and its covariance matrix  
      V(0,0)=hits[k]->covx;
      V(1,0)=V(0,1)=hits[k]->covxy;
      V(1,1)=hits[k]->covy;

      // Propagate state vector and covariance matrices to next measurement   
      int num_inc=int((z_-endz)/0.25);
      oldz=z_;
      newz=oldz;
      for (int j=0;j<num_inc-1;j++){
	dS=S-S0;
	newz=oldz-0.25;
	
	// Get dEdx for this step
	dedx=GetdEdx(0.14,S0(state_q_over_p,0),7.,14.,1.205e-3);
	
	// Step the seed state vector through the field and get the Jacobian 
	ds=StepJacobian(oldz,newz,S0,dedx,J);
	
	// Get the covariance matrix due to the multiple scattering
	GetProcessNoise(mass_hyp,ds,30420.,S0,Q); // use air for now
	
	// Rotate the covariance matrix and add the multiple scattering elements
	J_T=DMatrix(DMatrix::kTransposed,J);
	C=J*(C*J_T)+Q;
	
	// Update the "improved" state vector
	S=S0+J*dS;
	
	oldz=newz;			  
      }
      if (fabs(endz-oldz)>EPS){ 
	dS=S-S0;
	
	// Get dEdx for this step
	dedx=GetdEdx(0.14,S0(state_q_over_p,0),7.,14.,1.205e-3);
	
	// Step the seed state vector through the field and get the Jacobian 
	ds=StepJacobian(oldz,endz,S0,dedx,J);
	
	// Get the covariance matrix due to the multiple scattering
	GetProcessNoise(mass_hyp,ds,30420.,S0,Q); // use air for now
	
	// Rotate the covariance matrix and add the multiple scattering elements
	J_T=DMatrix(DMatrix::kTransposed,J);
	C=J*(C*J_T)+Q;
	
	// Update the "improved" state vector
	S=S0+J*dS;
	
	oldz=newz;	
      }  
      
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
      R(state_x,0)=M(0,0)-S(state_x,0);
      R(state_y,0)=M(1,0)-S(state_y,0);
      R_T=DMatrix(DMatrix::kTransposed,R);
      RC=V-H*(C*H_T);
 
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
      
      // increment z position
      z_=endz;
    }
    // Save the "best" values for S and C
    if (chisq<chisq_){
      SBest=S;
      CBest=C;
      chisq_=chisq;
    }

    // Swim back to first point to try again with "improved" S and C
    int num_inc=(int)((hits[0]->z-z_)/0.5);
    oldz=endz;
    for (int i=0;i<num_inc;i++){
      newz=oldz+0.5;  
      dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
      StepJacobian(oldz,newz,S,dedx,J);  
      J_T=DMatrix(DMatrix::kTransposed,J);     
      C=J*(C*J_T);
      oldz=newz;
    }

    // Scale the latest iteration of the covariance matrix to try to minimize
    // bias
    for (int i=0;i<5;i++)C(i,i)*=10000.;

  }

  // Propagate track to entrance to CDC endplate
  int  num_inc=(int)(z_-endplate_z);
  oldz=endz;
  endz=endplate_z;
  S=SBest; // Start with the state vector at the most upstream point
  C=CBest;
  for (int i=0;i<num_inc;i++){
    newz=oldz-0.5;  
    dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
    StepJacobian(oldz,newz,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);     
    C=J*(C*J_T);
    oldz=newz;
  }
  if (oldz!=endz){
    dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
    StepJacobian(oldz,endz,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);
    C=J*(C*J_T);    
  }

  // Next treat the CDC endplate
  oldz=endz;
  double r=0.;
  for (int i=0;i<4;i++){
    newz=oldz-endplate_dz/4.;
    r=sqrt(S(state_x,0)*S(state_x,0)+S(state_y,0)*S(state_y,0));
    if (r>endplate_rmin)
      dedx=GetdEdx(0.14,S(state_q_over_p,0),6.,12.,2.265);
    else
      dedx=GetdEdx(0.14,S(state_q_over_p,0),7.,14.,1.205e-3);
    StepJacobian(oldz,newz,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);  
    C=J*(C*J_T);
    oldz=newz;
  }

  // Final position for this leg
  x_=S(state_x,0);
  y_=S(state_y,0);
  z_=newz;

  return NOERROR;
}

// Propagate track to point of distance of closest approach to origin
jerror_t DKalmanFilter::ExtrapolateToVertex(double mass_hyp,DMatrix Sc,
					    DMatrix Cc){
  DMatrix Jc(5,5);  //.Jacobian matrix
  DMatrix JcT(5,5); // and its transpose
 
  // Initialize the position vector
  DVector3 pos(x_,y_,z_);
  DVector3 wire_pos(0,0,0);

  // Track propagation loop
  double r=pos.Perp();
  double ds=-0.25; // step along path in cm
  double r_old=r;
  while (Sc(state_z,0)>0. && r>BEAM_RADIUS){
    double cosl=cos(atan(Sc(state_tanl,0)));
    double dedx=0.;
    double q_over_p=Sc(state_q_over_pt,0)*cosl;
    if (r<target[1] && Sc(state_z,0)<target[2]){
      dedx=GetdEdx(mass_hyp,q_over_p,1.,1.,0.0708);
      ds=-0.1;
    }
    else
      dedx=GetdEdx(mass_hyp,q_over_p,7.,14.,1.205e-3);
    dedx=0;

    StepJacobian(pos,wire_pos,ds,Sc,dedx,Jc);
    JcT=DMatrix(DMatrix::kTransposed,Jc);
    Cc=Jc*(Cc*JcT);

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

  return NOERROR;
}
