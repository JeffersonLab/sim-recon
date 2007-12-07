#include "DRiemannFit.h"
#include <TDecompLU.h>
#include <math.h>

#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define Z_TARGET 65.0
#define EPS 1.0e-8

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

// Riemann Circle fit:  points on a circle in x,y project onto a plane cutting
// the circular paraboloid surface described by (x,y,x^2+y^2).  Therefore the
// task of fitting points in (x,y) to a circle is transormed to the task of
// fitting planes in (x,y, w=x^2+y^2) space
//
jerror_t DRiemannFit::FitCircle(double BeamRMS,DMatrix *Cov){  
  if (hits.size()==0) return RESOURCE_UNAVAILABLE;

  DMatrix X(hits.size()+1,3);
  DMatrix Xavg(1,3);
  DMatrix A(3,3);
  double B0,B1,B2,Q,Q1,R,sum,diff;
  double theta,lambda_min=1.e8,lambda[3];
  int smallest_eigenvalue=0;
  // Column and row vectors of ones
  DMatrix Ones(hits.size()+1,1),OnesT(1,hits.size()+1);
  DMatrix W_sum(1,1);
  DMatrix W(hits.size()+1,hits.size()+1);
  // Eigenvectors
  DMatrix N1(3,1);
  DMatrix N2(3,1);
  DMatrix N3(3,1);
  DMatrix VN(3,3);
 
  // Covariance matrix
  DMatrix CRPhi(hits.size()+1,hits.size()+1);
  if (Cov==NULL){
    Cov=new DMatrix(hits.size()+1,hits.size()+1);
    for (unsigned int i=0;i<hits.size();i++){
      double Phi=atan2(hits[i]->y,hits[i]->x);
      CRPhi(i,i)
      =(Phi*cos(Phi)-sin(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covx
      +(Phi*sin(Phi)+cos(Phi))*(Phi*sin(Phi)+cos(Phi))*hits[i]->covy
      +2.*(Phi*sin(Phi)+cos(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covxy;
    }
    CRPhi(hits.size(),hits.size())=BeamRMS*BeamRMS;
    Cov=&CRPhi;
  }
  else{
    for (unsigned int i=0;i<hits.size()+1;i++)
      for (unsigned int j=0;j<hits.size()+1;j++)
	CRPhi(i,j)=Cov->operator()(i, j);

  }
 
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
  Ones(hits.size(),0)=OnesT(0,hits.size())=1.;

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

  // The characteristic equation is 
  //   lambda^3+B2*lambda^2+lambda*B1+B0=0 
  //
  B2=-(A(0,0)+A(1,1)+A(2,2));
  B1=A(0,0)*A(1,1)-A(1,0)*A(0,1)+A(0,0)*A(2,2)-A(2,0)*A(0,2)+A(1,1)*A(2,2)
    -A(2,1)*A(1,2);
  B0=-A.Determinant();

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
 
  // First root
  lambda[0]=-B2/3.+sum;
  if (lambda[0]<lambda_min&&lambda[0]>EPS){
    lambda_min=lambda[0];
    smallest_eigenvalue=0;
  }

  // Second root
  lambda[1]=-B2/3.-sum/2.-sqrt(3.)/2.*diff;
  if (lambda[1]<lambda_min&&lambda[1]>EPS){
    lambda_min=lambda[1];
    smallest_eigenvalue=1;
  }

  // Third root
  lambda[2]=-B2/3.-sum/2.+sqrt(3.)/2.*diff;
  if (lambda[2]<lambda_min&&lambda[2]>EPS){
    lambda_min=lambda[2];
    smallest_eigenvalue=2;
  }

  // Normal vector to plane
  CalcNormal(A,lambda_min,N1);
  N[0]=N1(0,0);
  N[1]=N1(1,0);
  N[2]=N1(2,0);

  // Distance to origin
  double dist_to_origin=-(N1(0,0)*Xavg(0,0)+N1(1,0)*Xavg(0,1)+N1(2,0)*Xavg(0,2));

  // Center and radius of the circle
  xc=-N1(0,0)/2./N1(2,0);
  yc=-N1(1,0)/2./N1(2,0);
  rc=sqrt(1.-N1(2,0)*N1(2,0)-4.*dist_to_origin*N1(2,0))/2./fabs(N1(2,0));
  
  // p_trans and phi
	double Bz_avg=-2.0; 
	double q = +1.0;
	p_trans = q*Bz_avg*rc*qBr2p; // qBr2p converts to GeV/c
	phi = atan2(yc,xc) - M_PI_2;
	if(p_trans<0.0){
		p_trans = -p_trans;
	}
	if(phi<0)phi+=2.0*M_PI;
	if(phi>=2.0*M_PI)phi-=2.0*M_PI;

  return NOERROR;
}
 
// Linear regression to find charge
double DRiemannFit::GetCharge(double BeamRMS,DMatrix *CovR, DMatrix *CovRPhi){
  // Covariance matrices
  DMatrix CRPhi(hits.size()+1,hits.size()+1);
  DMatrix CR(hits.size()+1,hits.size()+1);
  if (CovRPhi==NULL){
    CovRPhi=new DMatrix(hits.size()+1,hits.size()+1);
    for (unsigned int i=0;i<hits.size();i++){
      double Phi=atan2(hits[i]->y,hits[i]->x);
      CRPhi(i,i)
      =(Phi*cos(Phi)-sin(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covx
      +(Phi*sin(Phi)+cos(Phi))*(Phi*sin(Phi)+cos(Phi))*hits[i]->covy
      +2.*(Phi*sin(Phi)+cos(Phi))*(Phi*cos(Phi)-sin(Phi))*hits[i]->covxy;
    }
    CRPhi(hits.size(),hits.size())=BeamRMS*BeamRMS;
    CovRPhi=&CRPhi;
  }
  else{
    for (unsigned int i=0;i<hits.size()+1;i++)
      for (unsigned int j=0;j<hits.size()+1;j++)
	CRPhi(i,j)=CovRPhi->operator()(i, j);
  }
  if (CovR==NULL){
    CovR=new DMatrix(hits.size()+1,hits.size()+1);
    for (unsigned int m=0;m<hits.size();m++){
      double Phi=atan2(hits[m]->y,hits[m]->x);
      CR(m,m)=cos(Phi)*cos(Phi)*hits[m]->covx
      +sin(Phi)*sin(Phi)*hits[m]->covy
      +2.*sin(Phi)*cos(Phi)*hits[m]->covxy;
    } 
    CR(hits.size(),hits.size())=BeamRMS*BeamRMS;
    CovR=&CR;
  }
  else{
    for (unsigned int i=0;i<hits.size()+1;i++)
      for (unsigned int j=0;j<hits.size()+1;j++)
	CR(i,j)=CovR->operator()(i, j);
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
    double r2=hit->x*hit->x+hit->y*hit->y;
    double var=(CRPhi(k,k)+phi_z*phi_z*CR(k,k))/r2;
    sumv+=1./var;
    sumy+=phi_z/var;
    sumx+=hit->z/var;
    sumxx+=hit->z*hit->z/var;
    sumxy+=phi_z*hit->z/var;
    phi_old=phi_z;
  }
  double slope=(sumv*sumxy-sumy*sumx)/(sumv*sumxx-sumx*sumx); 
 
  // Guess particle charge (+/-1);
  if (slope<0.) return -1.;

  return 1.;
}
