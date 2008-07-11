//************************************************************************
// DKalmanFilter.cc
//************************************************************************

#include "DKalmanFilter.h"
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

bool DKalmanHit_cmp(DKalmanHit_t *a, DKalmanHit_t *b){
  return a->z>b->z;
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
  x_=pos(0);
  y_=pos(1);
  z_=pos(2);
  tx_= mom(0)/mom(2);
  ty_= mom(1)/mom(2);
  q_over_p_=q/mom.Mag();
  
  return NOERROR;
}

// Return the momentum at the distance of closest approach to the origin.
void DKalmanFilter::GetMomentum(DVector3 &mom){
  /*
    double p=fabs(1./q_over_p_);
    double factor=sqrt(1.+tx_*tx_+ty_*ty_);
    mom.SetXYZ(p*tx_/factor,p*ty_/factor,p/factor);
  */

  mom.SetXYZ(cos(phi_)/q_over_pt_,sin(phi_)/q_over_pt_,tanl_/q_over_pt_);
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
jerror_t DKalmanFilter::CalcDerivAndJacobian(DVector3 &pos,double ds,
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
  double d=S(state_D,0);
  double z=S(state_z,0);

  //B-field and field gradient at (x,y,z)
  double Bx,By,Bz,dBxdx,dBxdy,dBxdz,dBydx,dBydy,dBydz,dBzdx,dBzdy,dBzdz;
  bfield->GetField(pos.x(),pos.y(),pos.z(), Bx, By, Bz);
  bfield->GetFieldGradient(pos.x(),pos.y(),pos.z(),dBxdx,dBxdy,dBxdz,dBydx,
			   dBydy,dBydz,dBzdx,dBzdy,dBzdz);
  double B=sqrt(Bx*Bx+By*By+Bz*Bz);

  // Derivative of S with respect to s
  D1(state_q_over_pt,0)=qBr2p*q_over_pt*q_over_pt*sinl*(By*cosphi-Bx*sinphi);
  D1(state_phi,0)=qBr2p*q_over_pt*(Bx*cosphi*sinl+By*sinphi*sinl-Bz*cosl);
  D1(state_tanl,0)=qBr2p*q_over_pt*(By*cosphi-Bx*sinphi)/cosl;
  D1(state_D,0)=0.;
  D1(state_z,0)=sinl;

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
jerror_t DKalmanFilter::StepJacobian(DVector3 &pos,double ds,
				     DMatrix &S,
				     double dEdx,DMatrix &J){
   // Initialize the Jacobian matrix
  J.Zero();
  for (int i=0;i<5;i++) J(i,i)=1.;
  // Matrices for intermediate steps
  DMatrix J1(5,5),J2(5,5),J3(5,5),J4(5,5);  
  DMatrix D1(5,1),D2(5,1),D3(5,1),D4(5,1);
  
  CalcDerivAndJacobian(pos,ds/2.,S,dEdx,J1,D1);
  CalcDerivAndJacobian(pos,0.,S+0.5*ds*D1,dEdx,J2,D2);
  J2=J2+0.5*J2*J1;
  CalcDerivAndJacobian(pos,ds/2.,S+0.5*ds*D2,dEdx,J3,D3);
  J3=J3+0.5*J3*J2;
  CalcDerivAndJacobian(pos,0.,S+ds*D3,dEdx,J4,D4);
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


// Routine that performs the main loop of the Kalman engine
jerror_t DKalmanFilter::KalmanLoop(double mass_hyp){
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
  DMatrix S(5,1),S0(5,1); //State vector
  DMatrix SBest(5,1);
  DMatrix dS(5,1);  // perturbation in state vector
  DMatrix C(5,5),CBest(5,5);   // Covariance matrix for state vector
  DMatrix InvV(2,2); // Inverse of error matrix

  // path length increment
  double ds=0;
  // Path length in active volume
  path_length=0;
  // Energy loss
  double dedx=0.;

  // Initialize the state vector 
  S0(0,0)=x_;
  S0(1,0)=y_;
  S0(2,0)=tx_;
  S0(3,0)=ty_;
  S0(4,0)=q_over_p_;

  // Initialize the covariance matrix
  for (unsigned int i=0;i<5;i++) C(i,i)=1;
  
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
	dedx=GetdEdx(0.14,S0(4,0),7.,14.,1.205e-3);
	
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
	dedx=GetdEdx(0.14,S0(4,0),7.,14.,1.205e-3);
	
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
      R(0,0)=M(0,0)-S(0,0);
      R(1,0)=M(1,0)-S(1,0);
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
      dedx=GetdEdx(0.14,S(4,0),7.,14.,1.205e-3);
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
    dedx=GetdEdx(0.14,S(4,0),7.,14.,1.205e-3);
    StepJacobian(oldz,newz,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);     
    C=J*(C*J_T);
    oldz=newz;
  }
  if (oldz!=endz){
    dedx=GetdEdx(0.14,S(4,0),7.,14.,1.205e-3);
    StepJacobian(oldz,endz,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);
    C=J*(C*J_T);    
  }

  // Next treat the CDC endplate
  oldz=endz;
  double r=0.;
  for (int i=0;i<4;i++){
    newz=oldz-endplate_dz/4.;
    r=sqrt(S(0,0)*S(0,0)+S(1,0)*S(1,0));
    if (r>endplate_rmin)
      dedx=GetdEdx(0.14,S(4,0),6.,12.,2.265);
    else
      dedx=GetdEdx(0.14,S(4,0),7.,14.,1.205e-3);
    StepJacobian(oldz,newz,S,dedx,J);  
    J_T=DMatrix(DMatrix::kTransposed,J);  
    C=J*(C*J_T);
    oldz=newz;
  }
  endz=newz;

  // Convert state vector into a version more useful for central region
  DMatrix Sc(5,1);
  DMatrix Cc(5,5);
  DMatrix Jc(5,5);  
  ConvertStateVector(endz,0.,0.,S,C,Sc,Cc);
  DVector3 pos(S(0,0),S(1,0),endz);
 
  // Propagate track to point of distance of closest approach to origin
  r=pos.Perp();
  ds=-1.; // 1 cm step along path
  double r_old=r;
  while (Sc(4,0)>0. && r>BEAM_RADIUS){
    dedx=0.;
    StepJacobian(pos,ds,Sc,dedx,Jc);
    r=pos.Perp();
    if (r>r_old) { // Passed through minimum doca to beamline
      break;
    }
    // Store current parameters
    r_old=r;
    SBest=Sc;
    x_=pos.x();
    y_=pos.y();
  }
  
  z_=SBest(4,0);
  phi_=SBest(1,0);
  q_over_pt_=SBest(0,0);
  tanl_=SBest(2,0);

  return NOERROR;
}
