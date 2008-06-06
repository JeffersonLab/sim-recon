#include "DRiemannFit.h"
#include <TDecompLU.h>
#include <math.h>

#include <iostream>
#include <algorithm>
using std::cerr;
using std::endl;

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define Z_TARGET 65.0
#define EPS 1.0e-8
#define MIN_TANL 2.0

// Boolean function for sorting hits
bool DRiemannFit_hit_cmp(DRiemannHit_t *a,DRiemannHit_t *b){
  return (a->z>b->z);
}

/// Add a hit to the list of hits using cylindrical coordinates
jerror_t DRiemannFit::AddHit(double r, double phi, double z)
{
  return AddHitXYZ(r*cos(phi), r*sin(phi), z);
}


 /// Add a hit to the list of hits using Cartesian coordinates
jerror_t DRiemannFit::AddHitXYZ(double x,double y, double z){
  DRiemannHit_t *hit = new DRiemannHit_t;
  hit->x = x;
  hit->y = y;
  hit->z = z;
  hit->covx=1.;
  hit->covy=1.;
  hit->covxy=0.;
  
  hits.push_back(hit);

  return NOERROR;
  
} 

/// Add a hit to the list of hits using Cartesian coordinates
jerror_t DRiemannFit::AddHit(double x,double y, double z,double covx,
				double covy, double covxy){
  DRiemannHit_t *hit = new DRiemannHit_t;
  hit->x = x;
  hit->y = y;
  hit->z = z;
  hit->covx=covx;
  hit->covy=covy;
  hit->covxy=covxy;
  
  hits.push_back(hit);

  return NOERROR;
}

// Fit sequence
jerror_t DRiemannFit::DoFit(double rc_input){
  jerror_t error=NOERROR;

  if (CovR_!=NULL) delete CovR_;
  if (CovRPhi_!=NULL) delete CovRPhi_; 
  CovR_=NULL;
  CovRPhi_=NULL;

  error=FitCircle(rc_input);
  error=FitLine();
  q=GetCharge();;

  return error;
}


// Calculate the (normal) eigenvector corresponding to the eigenvalue lambda
jerror_t DRiemannFit::CalcNormal(DMatrix A,double lambda,DMatrix &N){
  double sum=0;

  N(0,0)=1.;
  N(1,0)=N(0,0)*(A(1,0)*A(0,2)-(A(0,0)-lambda)*A(1,2))
    /(A(0,1)*A(2,1)-(A(1,1)-lambda)*A(0,2));
  N(2,0)=N(0,0)*(A(2,0)*(A(1,1)-lambda)-A(1,0)*A(2,1))
    /(A(1,2)*A(2,1)-(A(2,2)-lambda)*(A(1,1)-lambda));
  
  // Normalize: n1^2+n2^2+n3^2=1
  for (int i=0;i<3;i++){
    sum+=N(i,0)*N(i,0);
  }
  for (int i=0;i<3;i++){
    N(i,0)/=sqrt(sum);
  }

  return NOERROR;
}

// Riemann Circle fit with correction for non-normal track incidence
jerror_t DRiemannFit::FitCircle(double rc){
  // Covariance matrices
  DMatrix CRPhi(hits.size(),hits.size());
  DMatrix CR(hits.size(),hits.size());
  // Auxiliary matrices for correcting for non-normal track incidence to FDC
  // The correction is 
  //    CRPhi'= C*CRPhi*C+S*CR*S, where S(i,i)=R_i*kappa/2
  //                                and C(i,i)=sqrt(1-S(i,i)^2)  
  DMatrix C(hits.size(),hits.size());
  DMatrix S(hits.size(),hits.size());
  for (unsigned int i=0;i<hits.size();i++){
    double rtemp=sqrt(hits[i]->x*hits[i]->x+hits[i]->y*hits[i]->y); 
    double stemp=rtemp/4./rc;
    double ctemp=1.-stemp*stemp;
    if (ctemp>0){
      S(i,i)=stemp;
      C(i,i)=sqrt(ctemp);
    }
    else{
      S(i,i)=0.;
      C(i,i)=1.;
    }
  }

  // Covariance matrices
  if (CovRPhi_==NULL){
    CovRPhi_ = new DMatrix(hits.size(),hits.size());
    for (unsigned int i=0;i<hits.size();i++){
      double Phi=atan2(hits[i]->y,hits[i]->x);
      CovRPhi_->operator()(i,i)
      =(Phi*cos(Phi)-sin(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covx
      +(Phi*sin(Phi)+cos(Phi))*(Phi*sin(Phi)+cos(Phi))*hits[i]->covy
      +2.*(Phi*sin(Phi)+cos(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covxy;
    }
  }
  if (CovR_==NULL){
    CovR_= new DMatrix(hits.size(),hits.size());
    for (unsigned int m=0;m<hits.size();m++){
      double Phi=atan2(hits[m]->y,hits[m]->x);
      CovR_->operator()(m,m)=cos(Phi)*cos(Phi)*hits[m]->covx
	+sin(Phi)*sin(Phi)*hits[m]->covy
	+2.*sin(Phi)*cos(Phi)*hits[m]->covxy;
    } 
  }
  for (unsigned int i=0;i<hits.size();i++){
    for (unsigned int j=0;j<hits.size();j++){
      CR(i,j)=CovR_->operator()(i, j);
      CRPhi(i,j)=CovRPhi_->operator()(i, j);
    }
  }

  // Correction for non-normal incidence of track on FDC 
  CRPhi=C*CRPhi*C+S*CR*S;
  for (unsigned int i=0;i<hits.size();i++)
    for (unsigned int j=0;j<hits.size();j++)
      CovRPhi_->operator()(i,j)=CRPhi(i,j);
  return FitCircle();
}



// Riemann Circle fit:  points on a circle in x,y project onto a plane cutting
// the circular paraboloid surface described by (x,y,x^2+y^2).  Therefore the
// task of fitting points in (x,y) to a circle is transormed to the task of
// fitting planes in (x,y, w=x^2+y^2) space
//
jerror_t DRiemannFit::FitCircle(){  
  if (hits.size()==0) return RESOURCE_UNAVAILABLE;
  DMatrix X(hits.size(),3);
  DMatrix Xavg(1,3);
  DMatrix A(3,3);
  double B0,B1,B2,Q,Q1,R,sum,diff;
  double theta,lambda_min=0.;
  // Column and row vectors of ones
  DMatrix Ones(hits.size(),1),OnesT(1,hits.size());
  DMatrix W_sum(1,1);
  DMatrix W(hits.size(),hits.size());
  // Eigenvector
  DMatrix N1(3,1);

  // Make sure hit list is ordered in z
  std::sort(hits.begin(),hits.end(),DRiemannFit_hit_cmp);
 
  // Covariance matrix
  DMatrix CRPhi(hits.size(),hits.size());
  if (CovRPhi_==NULL){
    CovRPhi_ = new DMatrix(hits.size(),hits.size());
    for (unsigned int i=0;i<hits.size();i++){
      double Phi=atan2(hits[i]->y,hits[i]->x);
      CovRPhi_->operator()(i,i)
      =(Phi*cos(Phi)-sin(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covx
      +(Phi*sin(Phi)+cos(Phi))*(Phi*sin(Phi)+cos(Phi))*hits[i]->covy
      +2.*(Phi*sin(Phi)+cos(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covxy;
    }
  }
  for (unsigned int i=0;i<hits.size();i++)
    for (unsigned int j=0;j<hits.size();j++)
      CRPhi(i,j)=CovRPhi_->operator()(i, j);
 
  // The goal is to find the eigenvector corresponding to the smallest 
  // eigenvalue of the equation
  //            lambda=n^T (X^T W X - W_sum Xavg^T Xavg)n
  // where n is the normal vector to the plane slicing the cylindrical 
  // paraboloid described by the parameterization (x,y,w=x^2+y^2),
  // and W is the weight matrix, assumed for now to be diagonal.
  // In the absence of multiple scattering, W_sum is the sum of all the 
  // diagonal elements in W.

  for (unsigned int i=0;i<hits.size();i++){
    X(i,0)=hits[i]->x;
    X(i,1)=hits[i]->y;
    X(i,2)=hits[i]->x*hits[i]->x+hits[i]->y*hits[i]->y;
    Ones(i,0)=OnesT(0,i)=1.;
  }

  // Check that CRPhi is invertible 
  TDecompLU lu(CRPhi);
  if (lu.Decompose()==false){
    return UNRECOVERABLE_ERROR; // error placeholder
  }
  W=DMatrix(DMatrix::kInverted,CRPhi);
  W_sum=OnesT*(W*Ones);
  Xavg=(1./W_sum(0,0))*(OnesT*(W*X));
  
  A=DMatrix(DMatrix::kTransposed,X)*(W*X)
    -W_sum(0,0)*(DMatrix(DMatrix::kTransposed,Xavg)*Xavg);
  if(!A.IsValid())return UNRECOVERABLE_ERROR;

  // The characteristic equation is 
  //   lambda^3+B2*lambda^2+lambda*B1+B0=0 
  //
  B2=-(A(0,0)+A(1,1)+A(2,2));
  B1=A(0,0)*A(1,1)-A(1,0)*A(0,1)+A(0,0)*A(2,2)-A(2,0)*A(0,2)+A(1,1)*A(2,2)
    -A(2,1)*A(1,2);
  B0=-A.Determinant();
  if(B0==0 || !finite(B0))return UNRECOVERABLE_ERROR;

  // The roots of the cubic equation are given by 
  //        lambda1= -B2/3 + S+T
  //        lambda2= -B2/3 - (S+T)/2 + i sqrt(3)/2. (S-T)
  //        lambda3= -B2/3 - (S+T)/2 - i sqrt(3)/2. (S-T)
  // where we define some temporary variables:
  //        S= (R+sqrt(Q^3+R^2))^(1/3)
  //        T= (R-sqrt(Q^3+R^2))^(1/3)
  //        Q=(3*B1-B2^2)/9
  //        R=(9*B2*B1-27*B0-2*B2^3)/54
  //        sum=S+T;
  //        diff=i*(S-T)
  // We divide Q and R by a safety factor to prevent multiplying together 
  // enormous numbers that cause unreliable results.

  Q=(3.*B1-B2*B2)/9.e4; 
  R=(9.*B2*B1-27.*B0-2.*B2*B2*B2)/54.e6;
  Q1=Q*Q*Q+R*R;
  if (Q1<0) Q1=sqrt(-Q1);
  else{
    return VALUE_OUT_OF_RANGE;
  }

  // DeMoivre's theorem for fractional powers of complex numbers:  
  //      (r*(cos(theta)+i sin(theta)))^(p/q)
  //                  = r^(p/q)*(cos(p*theta/q)+i sin(p*theta/q))
  //
  double temp=100.*pow(R*R+Q1*Q1,0.16666666666666666667);
  theta=atan2(Q1,R)/3.;
  sum=2.*temp*cos(theta);
  diff=-2.*temp*sin(theta);
  // Third root
  lambda_min=-B2/3.-sum/2.+sqrt(3.)/2.*diff;
 
  // Normal vector to plane
  CalcNormal(A,lambda_min,N1);
  // Store in private array
  N[0]=N1(0,0);
  N[1]=N1(1,0);
  N[2]=N1(2,0);
 
  // Distance to origin
  dist_to_origin=-(N1(0,0)*Xavg(0,0)+N1(1,0)*Xavg(0,1)+N1(2,0)*Xavg(0,2));

  // Center and radius of the circle
  xc=-N1(0,0)/2./N1(2,0);
  yc=-N1(1,0)/2./N1(2,0);
  rc=sqrt(1.-N1(2,0)*N1(2,0)-4.*dist_to_origin*N1(2,0))/2./fabs(N1(2,0));

  return NOERROR;
}

// Charge-finding routine with corrected CRPhi (see above)
double DRiemannFit::GetCharge(double rc_input){
 // Covariance matrices
  DMatrix CRPhi(hits.size(),hits.size());
  DMatrix CR(hits.size(),hits.size());
  // Auxiliary matrices for correcting for non-normal track incidence to FDC
  // The correction is 
  //    CRPhi'= C*CRPhi*C+S*CR*S, where S(i,i)=R_i*kappa/2
  //                                and C(i,i)=sqrt(1-S(i,i)^2)  
  DMatrix C(hits.size(),hits.size());
  DMatrix S(hits.size(),hits.size());
  for (unsigned int i=0;i<hits.size();i++){
    double rtemp=sqrt(hits[i]->x*hits[i]->x+hits[i]->y*hits[i]->y); 
    double stemp=rtemp/4./rc_input;
    double ctemp=1.-stemp*stemp;
    if (ctemp>0){
      S(i,i)=stemp;
      C(i,i)=sqrt(ctemp);
    }
    else{
      S(i,i)=0.;
      C(i,i)=1;      
    }
  }
  
  // Covariance matrices
  if (CovRPhi_==NULL){
    CovRPhi_ = new DMatrix(hits.size(),hits.size());
    for (unsigned int i=0;i<hits.size();i++){
      double Phi=atan2(hits[i]->y,hits[i]->x);
      CovRPhi_->operator()(i,i)
      =(Phi*cos(Phi)-sin(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covx
      +(Phi*sin(Phi)+cos(Phi))*(Phi*sin(Phi)+cos(Phi))*hits[i]->covy
      +2.*(Phi*sin(Phi)+cos(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covxy;
    }
  }
  if (CovR_==NULL){
    CovR_= new DMatrix(hits.size(),hits.size());
    for (unsigned int m=0;m<hits.size();m++){
      double Phi=atan2(hits[m]->y,hits[m]->x);
      CovR_->operator()(m,m)=cos(Phi)*cos(Phi)*hits[m]->covx
	+sin(Phi)*sin(Phi)*hits[m]->covy
	+2.*sin(Phi)*cos(Phi)*hits[m]->covxy;
    } 
  }
  for (unsigned int i=0;i<hits.size();i++){
    for (unsigned int j=0;j<hits.size();j++){
      CR(i,j)=CovR_->operator()(i, j);
      CRPhi(i,j)=CovRPhi_->operator()(i, j);
    }
  }
  // Correction for non-normal incidence of track on FDC 
  CRPhi=C*CRPhi*C+S*CR*S; 
  for (unsigned int i=0;i<hits.size();i++)
    for (unsigned int j=0;j<hits.size();j++)
      CovRPhi_->operator()(i,j)=CRPhi(i,j);
  return GetCharge();
}

 
// Linear regression to find charge
double DRiemannFit::GetCharge(){
  // Covariance matrices
  DMatrix CRPhi(hits.size(),hits.size());
  DMatrix CR(hits.size(),hits.size());
  if (CovRPhi_==NULL){
    CovRPhi_ = new DMatrix(hits.size(),hits.size());
    for (unsigned int i=0;i<hits.size();i++){
      double Phi=atan2(hits[i]->y,hits[i]->x);
      CovRPhi_->operator()(i,i)
	=(Phi*cos(Phi)-sin(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covx
      +(Phi*sin(Phi)+cos(Phi))*(Phi*sin(Phi)+cos(Phi))*hits[i]->covy
      +2.*(Phi*sin(Phi)+cos(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covxy;
    }
  }
  if (CovR_==NULL){
    CovR_= new DMatrix(hits.size(),hits.size());
    for (unsigned int m=0;m<hits.size();m++){
      double Phi=atan2(hits[m]->y,hits[m]->x);
      CovR_->operator()(m,m)=cos(Phi)*cos(Phi)*hits[m]->covx
	+sin(Phi)*sin(Phi)*hits[m]->covy
	+2.*sin(Phi)*cos(Phi)*hits[m]->covxy;
    } 
  }
  for (unsigned int i=0;i<hits.size();i++){
    for (unsigned int j=0;j<hits.size();j++){
      CR(i,j)=CovR_->operator()(i, j);
      CRPhi(i,j)=CovRPhi_->operator()(i, j);
    }
  }

  double phi_old=atan2(hits[0]->y,hits[0]->x);
  double sumv=0,sumy=0,sumx=0,sumxx=0,sumxy=0;
  for (unsigned int k=0;k<hits.size();k++){
    DRiemannHit_t *hit=hits[k];   
    double phi_z=atan2(hit->y,hit->x);

    // Check for problem regions near +pi and -pi
    if (fabs(phi_z-phi_old)>M_PI){  
      if (phi_old<0) phi_z-=2.*M_PI;
      else phi_z+=2.*M_PI;
    }
    double inv_var=(hit->x*hit->x+hit->y*hit->y)
      /(CRPhi(k,k)+phi_z*phi_z*CR(k,k));
    sumv+=inv_var;
    sumy+=phi_z*inv_var;
    sumx+=hit->z*inv_var;
    sumxx+=hit->z*hit->z*inv_var;
    sumxy+=phi_z*hit->z*inv_var;
    phi_old=phi_z;
  }
  double slope=(sumv*sumxy-sumy*sumx)/(sumv*sumxx-sumx*sumx); 
 
  // Guess particle charge (+/-1);
  if (slope<0.) return -1.;

  return 1.;
}


// Riemann Line fit: linear regression of s on z to determine the tangent of 
// the dip angle and the z position of the closest approach to the beam line.
// Computes intersection points along the helical path.
//
jerror_t DRiemannFit::FitLine(){
  // Get covariance matrix 
  DMatrix CR(hits.size(),hits.size());
  if (CovR_==NULL){ 
    CovR_= new DMatrix(hits.size(),hits.size());
    for (unsigned int m=0;m<hits.size();m++){
      double Phi=atan2(hits[m]->y,hits[m]->x);
      CovR_->operator()(m,m)=cos(Phi)*cos(Phi)*hits[m]->covx
      +sin(Phi)*sin(Phi)*hits[m]->covy
      +2.*sin(Phi)*cos(Phi)*hits[m]->covxy;
    } 
  }
  for (unsigned int i=0;i<hits.size();i++)
    for (unsigned int j=0;j<hits.size();j++)
      CR(i,j)=CovR_->operator()(i, j);

  // Fill vector of intersection points 
  double x_int0,temp,y_int0;
  double denom= N[0]*N[0]+N[1]*N[1];
  double numer;
  vector<int>bad(hits.size());
  int numbad=0;
  // Clear old projection vector
  projections.clear();
  for (unsigned int m=0;m<hits.size();m++){
    double r2=hits[m]->x*hits[m]->x+hits[m]->y*hits[m]->y;
    numer=dist_to_origin+r2*N[2];
    DRiemannHit_t *temphit = new DRiemannHit_t;
    temphit->z=hits[m]->z;

    if (r2==0){
      temphit->x=0.;
      temphit->y=0.;
      //      bad[m]=1;
    }
    else{
      x_int0=-N[0]*numer/denom;
      y_int0=-N[1]*numer/denom;
      temp=denom*r2-numer*numer;
      if (temp<0){  // Skip point if the intersection gives nonsense
	bad[m]=1;
	numbad++;
	temphit->x=x_int0;
	temphit->y=y_int0;
	projections.push_back(temphit);
	continue;
      }
      temp=sqrt(temp)/denom;
      
      // Choose sign of square root based on proximity to actual measurements
      double diffx1=x_int0+N[1]*temp-hits[m]->x;
      double diffy1=y_int0-N[0]*temp-hits[m]->y;
      double diffx2=x_int0-N[1]*temp-hits[m]->x;
      double diffy2=y_int0+N[0]*temp-hits[m]->y;
      if (diffx1*diffx1+diffy1*diffy1 > diffx2*diffx2+diffy2*diffy2){
	temphit->x=x_int0-N[1]*temp;
	temphit->y=y_int0+N[0]*temp;
      }
      else{
	temphit->x=x_int0+N[1]*temp;
	temphit->y=y_int0-N[0]*temp;
      }
    }
    projections.push_back(temphit);
  }  
  
  // All arc lengths are measured relative to some reference plane with a hit.
  // Don't use a "bad" hit for the reference...
  unsigned int start=0;
  for (unsigned int i=0;i<bad.size();i++){
    if (!bad[i]){
      start=i;
      break;
    }
  }

  // Linear regression to find z0, tanl   
  unsigned int n=projections.size();
  double sumv=0.,sumx=0.,sumy=0.,sumxx=0.,sumxy=0.;
  double sperp=0., chord=0,ratio=0, Delta;
  for (unsigned int k=start;k<n;k++){
    if (!bad[k])
      {
      double diffx=projections[k]->x-projections[start]->x;
      double diffy=projections[k]->y-projections[start]->y;
      chord=sqrt(diffx*diffx+diffy*diffy);
      ratio=chord/2./rc; 
      // Make sure the argument for the arcsin does not go out of range...
      if (ratio>1.) 
	sperp=2.*rc*(M_PI/2.);
      else
	sperp=2.*rc*asin(ratio);
      // Assume errors in s dominated by errors in R 
      sumv+=1./CR(k,k);
      sumy+=sperp/CR(k,k);
      sumx+=projections[k]->z/CR(k,k);
      sumxx+=projections[k]->z*projections[k]->z/CR(k,k);
      sumxy+=sperp*projections[k]->z/CR(k,k);
    }
  }
  chord=sqrt(projections[start]->x*projections[start]->x
	     +projections[start]->y*projections[start]->y);
  ratio=chord/2./rc; 
  // Make sure the argument for the arcsin does not go out of range...
  if (ratio>1.) 
    sperp=2.*rc*(M_PI/2.);
  else
    sperp=2.*rc*asin(ratio);
  Delta=sumv*sumxx-sumx*sumx;
  // Track parameter tan(lambda)
  tanl=-Delta/(sumv*sumxy-sumy*sumx); 

  // Vertex position
  zvertex=projections[start]->z-sperp*tanl;
  double zvertex_temp=projections[start]->z-(2.*rc*M_PI-sperp)*tanl;
  // Choose vertex z based on proximity of projected z-vertex to the center 
  // of the target
  if (fabs(zvertex-Z_TARGET)>fabs(zvertex_temp-Z_TARGET)) zvertex=zvertex_temp;

  return NOERROR;
}

