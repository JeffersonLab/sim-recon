// Set of routines for fitting tracks in both the cdc and fdc assuming a
// uniform magnetic field.  The track therefore is assumed to follow a helical
// path.

#include <iostream>
#include <algorithm>
using namespace std;

#include <TDecompLU.h>
#include <math.h>

#include "DHelicalFit.h"
#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define Z_VERTEX 65.0
#define Z_MIN 45.0
#define Z_MAX 85.0

// The following is for sorting hits by z
class DHFHitLessThanZ{
	public:
		bool operator()(DHFHit_t* const &a, DHFHit_t* const &b) const {
			return a->z < b->z;
		}
};

bool DHFHitLessThanZ_C(DHFHit_t* const &a, DHFHit_t* const &b) {
	return a->z < b->z;
}
bool RiemannFit_hit_cmp(DHFHit_t *a,DHFHit_t *b){
  return (a->z>b->z);
}



//-----------------
// DHelicalFit (Constructor)
//-----------------
DHelicalFit::DHelicalFit(void)
{
	x0 = y0 = r0 = 0;
	chisq = 0;
	chisq_source = NOFIT;
	bfield = NULL;
	
	// For Riemann Fit
	CovR_=NULL;
	CovRPhi_=NULL;
	c_origin=0.;
	normal.SetXYZ(0.,0.,0.);
}

//-----------------
// DHelicalFit (Constructor)
//-----------------
DHelicalFit::DHelicalFit(const DHelicalFit &fit)
{
	Copy(fit);
}
//-----------------
// Copy
//-----------------
void DHelicalFit::Copy(const DHelicalFit &fit)
{
	x0 = fit.x0;
	y0 = fit.y0;
	q = fit.q;
	p = fit.p;
	p_trans = fit.p_trans;
	phi = fit.phi;
	theta = fit.theta;
	z_vertex = fit.z_vertex;
	chisq = fit.chisq;
	dzdphi = fit.dzdphi;
	chisq_source = fit.chisq_source;
	bfield = fit.GetMagneticFieldMap();
	Bz_avg = fit.GetBzAvg();
	z_mean = fit.GetZMean();
	phi_mean = fit.GetPhiMean();
	
	const vector<DHFHit_t*> myhits = fit.GetHits();
	for(unsigned int i=0; i<myhits.size(); i++){
		DHFHit_t *a = new DHFHit_t;
		*a = *myhits[i];
		hits.push_back(a);
	}
	if (fit.CovR_!=NULL){
	  CovR_ = new DMatrix(hits.size(),hits.size());
	  for (unsigned int i=0;i<hits.size();i++)
	    for (unsigned int j=0;j<hits.size();j++)
	      CovR_->operator()(i,j)=fit.CovR_->operator()(i,j);
	}
	else CovR_=NULL; 
	if (fit.CovRPhi_!=NULL){
	  CovRPhi_ = new DMatrix(hits.size(),hits.size());
	  for (unsigned int i=0;i<hits.size();i++)
	    for (unsigned int j=0;j<hits.size();j++)
	      CovRPhi_->operator()(i,j)=fit.CovRPhi_->operator()(i,j);
	}
	else CovRPhi_=NULL; 
}

//-----------------
// operator= (Assignment operator)
//-----------------
DHelicalFit& DHelicalFit::operator=(const DHelicalFit& fit)
{
	if(this == &fit)return *this;
	Copy(fit);

	return *this;
}

//-----------------
// DHelicalFit (Destructor)
//-----------------
DHelicalFit::~DHelicalFit()
{
	Clear();
}
//-----------------
// AddHit --  add an FDC hit to the list
//-----------------
jerror_t DHelicalFit::AddHit(const DFDCPseudo *fdchit){
  DHFHit_t *hit = new DHFHit_t;
  hit->x = fdchit->x;
  hit->y = fdchit->y;
  hit->z = fdchit->wire->origin.z();
  hit->covx=fdchit->covxx;
  hit->covy=fdchit->covyy;
  hit->covxy=fdchit->covxy;
  hit->chisq = 0.0;
  hits.push_back(hit);
  
  return NOERROR;

}


//-----------------
// AddHit
//-----------------
jerror_t DHelicalFit::AddHit(float r, float phi, float z)
{
	/// Add a hit to the list of hits using cylindrical coordinates

	return AddHitXYZ(r*cos(phi), r*sin(phi), z);
}

//-----------------
// AddHitXYZ
//-----------------
jerror_t DHelicalFit::AddHitXYZ(float x, float y, float z)
{
	/// Add a hit to the list of hits using Cartesian coordinates
	DHFHit_t *hit = new DHFHit_t;
	hit->x = x;
	hit->y = y;
	hit->z = z;
	hit->covx=1.;
	hit->covy=1.;
	hit->covxy=0.;
	hit->chisq = 0.0;
	hits.push_back(hit);

	return NOERROR;
}

// Add a hit to the list of hits using Cartesian coordinates
jerror_t DHelicalFit::AddHitXYZ(float x,float y, float z,float covx,
				float covy, float covxy){
  DHFHit_t *hit = new DHFHit_t;
  hit->x = x;
  hit->y = y;
  hit->z = z;
  hit->covx=covx;
  hit->covy=covy;
  hit->covxy=covxy;
  hits.push_back(hit);

  return NOERROR;
}


//-----------------
// PruneHit
//-----------------
jerror_t DHelicalFit::PruneHit(int idx)
{
	/// Remove the hit specified by idx from the list
	/// of hits. The value of idx can be anywhere from
	/// 0 to GetNhits()-1.
	if(idx<0 || idx>=(int)hits.size())return VALUE_OUT_OF_RANGE;

	delete hits[idx];
	hits.erase(hits.begin() + idx);

	return NOERROR;
}

//-----------------
// Clear
//-----------------
jerror_t DHelicalFit::Clear(void)
{
  /// Remove covariance matrices
  if (CovR_!=NULL) delete CovR_;
  if (CovRPhi_!=NULL) delete CovRPhi_; 
 
  /// Remove all hits
  for(unsigned int i=0; i<hits.size(); i++)delete hits[i];
  hits.clear();
 
  // Remove projections 
  for (unsigned int i=0;i<projections.size();i++)
    delete projections[i];
  projections.clear();

  return NOERROR;
}

//-----------------
// PrintChiSqVector
//-----------------
jerror_t DHelicalFit::PrintChiSqVector(void) const
{
	/// Dump the latest chi-squared vector to the screen.
	/// This prints the individual hits' chi-squared
	/// contributions in the order in which the hits were
	/// added. See PruneHits() for more detail.

	cout<<"Chisq vector from DHelicalFit: (source=";
	switch(chisq_source){
		case NOFIT:		cout<<"NOFIT";		break;
		case CIRCLE:	cout<<"CIRCLE";	break;
		case TRACK:		cout<<"TRACK";		break;
	};
	cout<<")"<<endl;
	cout<<"----------------------------"<<endl;

	for(unsigned int i=0;i<hits.size();i++){
		cout<<i<<"  "<<hits[i]->chisq<<endl;
	}
	cout<<"Total: "<<chisq<<endl<<endl;

	return NOERROR;
}

//-----------------
// FitCircle
//-----------------
jerror_t DHelicalFit::FitCircle(void)
{
	/// Fit the current set of hits to a circle
	///
	/// This takes the hits which have been added thus far via one
	/// of the AddHit() methods and fits them to a circle.
	/// The alogorithm used here calculates the parameters directly using
	/// a technique very much like linear regression. The key assumptions
	/// are:
	/// 1. The magnetic field is uniform and along z so that the projection
	///    of the track onto the X/Y plane will fall on a circle
	///    (this also implies no multiple-scattering or energy loss)
	/// 2. The vertex is at the target center (i.e. 0,0 in the coordinate
	///    system of the points passed to us.
	///
	/// IMPORTANT: The value of phi which results from this assumes
	/// the particle was positively charged. If the particle was
	/// negatively charged, then phi will be 180 degrees off. To
	/// correct this, one needs z-coordinate information to determine
	/// the sign of the charge.
	///
	/// ALSO IMPORTANT: This assumes a charge of +1 for the particle. If
	/// the particle actually has a charge of +2, then the resulting
	/// value of p_trans will be half of what it should be.

	float alpha=0.0, beta=0.0, gamma=0.0, deltax=0.0, deltay=0.0;
	chisq_source = NOFIT; // in case we return early
	
	// Loop over hits to calculate alpha, beta, gamma, and delta
	// if a magnetic field map was given, use it to find average Z B-field
	DHFHit_t *a = NULL;
	for(unsigned int i=0;i<hits.size();i++){
		a = hits[i];
		float x=a->x;
		float y=a->y;
		alpha += x*x;
		beta += y*y;
		gamma += x*y;
		deltax += x*(x*x+y*y)/2.0;
		deltay += y*(x*x+y*y)/2.0;		
	}
	
	// Calculate x0,y0 - the center of the circle
	double denom = alpha*beta-gamma*gamma;
	if(fabs(denom)<1.0E-20)return UNRECOVERABLE_ERROR;
	x0 = (deltax*beta-deltay*gamma)/denom;
	//y0 = (deltay-gamma*x0)/beta;
	y0 = (deltay*alpha-deltax*gamma)/denom;
	
	// Calculate the phi angle traversed by the track from the
	// origin to the last hit. NOTE: this can be off by a multiple
	// of 2pi!
	double delta_phi=0.0;
	if(a){ // a should be pointer to last hit from above loop
		delta_phi = atan2(a->y-y0, a->x-x0);
		if(delta_phi<0.0)delta_phi += 2.0*M_PI;
	}

	// Momentum depends on magnetic field. If bfield has been
	// set, we should use it to determine an average value of Bz
	// for this track. Otherwise, assume -2T.
	// Also assume a singly charged track (i.e. q=+/-1)
	// The sign of the charge will be determined below.
	Bz_avg=-2.0; 
	q = +1.0;
	r0 = sqrt(x0*x0 + y0*y0);
	p_trans = q*Bz_avg*r0*qBr2p; // qBr2p converts to GeV/c
	phi = atan2(y0,x0) - M_PI_2;
	if(p_trans<0.0){
		p_trans = -p_trans;
	}
	if(phi<0)phi+=2.0*M_PI;
	if(phi>=2.0*M_PI)phi-=2.0*M_PI;
	
	// Calculate the chisq
	ChisqCircle();
	chisq_source = CIRCLE;

	return NOERROR;
}

//-----------------
// ChisqCircle
//-----------------
double DHelicalFit::ChisqCircle(void)
{
	/// Calculate the chisq for the hits and current circle
	/// parameters.
	/// NOTE: This does not return the chi2/dof, just the
	/// raw chi2 with errs set to 1.0
	chisq = 0.0;
	for(unsigned int i=0;i<hits.size();i++){
		DHFHit_t *a = hits[i];
		float x = a->x - x0;
		float y = a->y - y0;
		float c = sqrt(x*x + y*y) - r0;
		c *= c;
		a->chisq = c;
		chisq+=c;
	}
	
	// Do NOT divide by DOF
	return chisq;
}

//-----------------
// FitCircleRiemann
//-----------------
jerror_t DHelicalFit::FitCircleRiemann(float BeamRMS)
{
	chisq_source = NOFIT; // in case we return early
	
	// Fake point at origin
	float beam_var=BeamRMS*BeamRMS;
	AddHitXYZ(0.,0.,Z_VERTEX,beam_var,beam_var,0.);

	jerror_t err = FitCircleRiemann();
	if(err!=NOERROR)return err;
	
	// Momentum depends on magnetic field. If bfield has been
	// set, we should use it to determine an average value of Bz
	// for this track. Otherwise, assume -2T.
	// Also assume a singly charged track (i.e. q=+/-1)
	// The sign of the charge will be determined below.
	Bz_avg=-2.0; 
	q = +1.0;
	p_trans = q*Bz_avg*r0*qBr2p; // qBr2p converts to GeV/c
	if(p_trans<0.0){
		p_trans = -p_trans;
	}
	phi=atan2(-x0,y0);  
	if(phi<0)phi+=2.0*M_PI;
        if(phi>=2.0*M_PI)phi-=2.0*M_PI;
	
	// Normal vector for plane intersecting Riemann surface 
	normal.SetXYZ(N[0],N[1],N[2]);

	return NOERROR;
}

//-----------------
// FitCircleRiemann
//-----------------
jerror_t DHelicalFit::FitCircleRiemann(void){  
/// Riemann Circle fit:  points on a circle in x,y project onto a plane cutting
/// the circular paraboloid surface described by (x,y,x^2+y^2).  Therefore the
/// task of fitting points in (x,y) to a circle is transormed to the task of
/// fitting planes in (x,y, w=x^2+y^2) space
///
  DMatrix X(hits.size(),3);
  DMatrix Xavg(1,3);
  DMatrix A(3,3);
  double B0,B1,B2,Q,Q1,R,sum,diff;
  double theta,lambda_min=0.;
  // Column and row vectors of ones
  DMatrix Ones(hits.size(),1),OnesT(1,hits.size());
  DMatrix W_sum(1,1);
  DMatrix W(hits.size(),hits.size());

  // Make sure hit list is ordered in z
  std::sort(hits.begin(),hits.end(),RiemannFit_hit_cmp);
 
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
 
  // Calculate the (normal) eigenvector corresponding to the eigenvalue lambda
  N[0]=1.;
  N[1]=(A(1,0)*A(0,2)-(A(0,0)-lambda_min)*A(1,2))
    /(A(0,1)*A(2,1)-(A(1,1)-lambda_min)*A(0,2));
  N[2]=(A(2,0)*(A(1,1)-lambda_min)-A(1,0)*A(2,1))
    /(A(1,2)*A(2,1)-(A(2,2)-lambda_min)*(A(1,1)-lambda_min));
  
  // Normalize: n1^2+n2^2+n3^2=1
  sum=0.;
  for (int i=0;i<3;i++){
    sum+=N[i]*N[i];
  }
  for (int i=0;i<3;i++){
    N[i]/=sqrt(sum);
  }
  /*
  DVector3 testn(1.,(A(1,0)*A(0,2)-(A(0,0)-lambda_min)*A(1,2))
		 /(A(0,1)*A(2,1)-(A(1,1)-lambda_min)*A(0,2)),
		 (A(2,0)*(A(1,1)-lambda_min)-A(1,0)*A(2,1))
		 /(A(1,2)*A(2,1)-(A(2,2)-lambda_min)*(A(1,1)-lambda_min)));
  testn.SetMag(1.0);
  testn.Print();
  printf("N %f %f %f\n",N[0],N[1],N[2]);
  */

  // Distance to origin
  c_origin=-(N[0]*Xavg(0,0)+N[1]*Xavg(0,1)+N[2]*Xavg(0,2));

  // Center and radius of the circle
  double one_over_2Nz=1./(2.*N[2]);
  //x0=-N[0]/2./N[2];
  //y0=-N[1]/2./N[2];
  //r0=sqrt(1.-N[2]*N[2]-4.*c_origin*N[2])/2./fabs(N[2]);
  x0=-N[0]*one_over_2Nz;
  y0=-N[1]*one_over_2Nz;
  r0=sqrt(1.-N[2]*N[2]-4.*c_origin*N[2])*fabs(one_over_2Nz);
 
  // Phi value at "vertex"
  phi=atan2(-x0,y0);  
  if(phi<0)phi+=2.0*M_PI;
  if(phi>=2.0*M_PI)phi-=2.0*M_PI;

  // Calculate the chisq
  ChisqCircle();
  chisq_source = CIRCLE;

  return NOERROR;
}


//-----------------
// FitCircleRiemannCorrected
//-----------------
// Riemann Circle fit with correction for non-normal track incidence
//
jerror_t DHelicalFit::FitCircleRiemannCorrected(float rc){
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
    S(i,i)=0.;
    C(i,i)=1.;
    if (rc>0){
      double rtemp=sqrt(hits[i]->x*hits[i]->x+hits[i]->y*hits[i]->y); 
      double stemp=rtemp/4./rc;
      double ctemp=1.-stemp*stemp;
      if (ctemp>0){
	S(i,i)=stemp;
	C(i,i)=sqrt(ctemp);
      }
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
  return FitCircleRiemann();
}


//-----------------
// GetChargeRiemann
//-----------------
// Charge-finding routine with corrected CRPhi (see above)
//
jerror_t DHelicalFit::GetChargeRiemann(float rc_input){
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
  return GetChargeRiemann();
}

//-----------------
// GetChargeRiemann
//-----------------
// Linear regression to find charge
//
jerror_t DHelicalFit::GetChargeRiemann(){
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
    DHFHit_t *hit=hits[k];   
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
  q=+1.;
  if (slope<0.) q= -1.;

  return NOERROR;
}

//-----------------
// FitLineRiemann
//-----------------
jerror_t DHelicalFit::FitLineRiemann(){
/// Riemann Line fit: linear regression of s on z to determine the tangent of 
/// the dip angle and the z position of the closest approach to the beam line.
/// Computes intersection points along the helical path.
///
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
    numer=c_origin+r2*N[2];
    DHFHit_t *temphit = new DHFHit_t;
    temphit->z=hits[m]->z;

    if (r2==0){
      temphit->x=0.;
      temphit->y=0.;
      //      bad[m]=1;
    }
    else{
      double ratio=numer/denom;
      x_int0=-N[0]*ratio;
      y_int0=-N[1]*ratio;
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
      double deltax=N[1]*temp;
      double deltay=-N[0]*temp;
      double x1=x_int0+deltax;
      double y1=y_int0+deltay;
      double x2=x_int0-deltax;
      double y2=y_int0-deltay;
      double diffx1=x1-hits[m]->x;
      double diffy1=y1-hits[m]->y;
      double diffx2=x2-hits[m]->x;
      double diffy2=y2-hits[m]->y;
      if (diffx1*diffx1+diffy1*diffy1 > diffx2*diffx2+diffy2*diffy2){
	temphit->x=x2;
	temphit->y=y2;
      }
      else{
	temphit->x=x1;
	temphit->y=y1;
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
  double sperp=0.,sperp_old=0., ratio=0, Delta;
  double z_last=0.,z=0.;
  for (unsigned int k=start;k<n;k++){
    if (!bad[k]){
      sperp_old=sperp;
      z_last=z;

      double diffx=projections[k]->x-projections[start]->x;
      double diffy=projections[k]->y-projections[start]->y;   
      ratio=sqrt(diffx*diffx+diffy*diffy)/(2.*r0); 
      // Make sure the argument for the arcsin does not go out of range...
      if (ratio>1.) 
	sperp=2.*r0*(M_PI/2.);
      else
	sperp=2.*r0*asin(ratio);
      if (sperp-sperp_old<0.){
	if (k==n-1) sperp=2.*r0*M_PI-sperp;
      }
      z=projections[k]->z;

      // Assume errors in s dominated by errors in R 
      double weight=1./CR(k,k);
      sumv+=weight;
      sumy+=sperp*weight;
      sumx+=projections[k]->z*weight;
      sumxx+=projections[k]->z*projections[k]->z*weight;
      sumxy+=sperp*projections[k]->z*weight;
    }
  }
  Delta=sumv*sumxx-sumx*sumx;
  // Track parameters tan(lambda) and z-vertex
  tanl=-Delta/(sumv*sumxy-sumy*sumx); 
  z_vertex=(sumxx*sumy-sumx*sumxy)/Delta;

  // Use the last z and s values to estimate the vertex position if the first
  // method gave a result beyond the extent of the target
  if (z_vertex<Z_MIN || z_vertex>Z_MAX){
    sperp-=sperp_old;
    double myz_vertex=z_last-sperp*tanl;
    if (fabs(myz_vertex-Z_VERTEX)<fabs(z_vertex-Z_VERTEX)) z_vertex=myz_vertex;
  }
  theta=M_PI_2-atan(tanl);

  return NOERROR;
}

//-----------------
// FitCircleStraightTrack
//-----------------
jerror_t DHelicalFit::FitCircleStraightTrack(void)
{
	/// This is a last resort!
	/// The circle fit can fail for high momentum tracks that are nearly
	/// straight tracks. In these cases (when pt~2GeV or more) there
	/// is not enough position resolution with wire positions only
	/// to measure the curvature of the track. 
	/// We can though, fit the X/Y points with a straight line in order
	/// to get phi and project out the stereo angles. 
	///
	/// For the radius of the circle (i.e. p_trans) we loop over momenta
	/// from 1 to 9GeV and just use the one with the smallest chisq.

	// Average the x's and y's of the individual hits. We should be 
	// able to do this via linear regression, but I can't get it to
	// work right now and I'm under the gun to get ready for a review
	// so I take the easy way. Note that we don't average phi's since
	// that will cause problems at the 0/2pi boundary.
	double X=0.0, Y=0.0;
	DHFHit_t *a = NULL;
	for(unsigned int i=0;i<hits.size();i++){
		a = hits[i];
		double r = sqrt(pow((double)a->x,2.0) + pow((double)a->y, 2.0));
		// weight by r to give outer points more influence. Note that
		// we really are really weighting by r^2 since x and y already
		// have a magnitude component.
		X += a->x*r; 
		Y += a->y*r;
	}
	phi = atan2(Y,X);
	if(phi<0)phi+=2.0*M_PI;
	if(phi>=2.0*M_PI)phi-=2.0*M_PI;

	// Search the chi2 space for values for p_trans, x0, ...
	SearchPtrans(9.0, 0.5);

#if 0
	// We do a simple linear regression here to find phi. This is
	// simplified by the intercept being zero (i.e. the track
	// passes through the beamline).
	float Sxx=0.0, Syy=0.0, Sxy=0.0;
	chisq_source = NOFIT; // in case we return early
	
	// Loop over hits to calculate Sxx, Syy, and Sxy
	DHFHit_t *a = NULL;
	for(unsigned int i=0;i<hits.size();i++){
		a = hits[i];
		float x=a->x;
		float y=a->y;
		Sxx += x*x;
		Syy += y*y;
		Sxy += x*y;
	}
	double A = 2.0*Sxy;
	double B = Sxx - Syy;
	phi = B>A ? atan2(A,B)/2.0 : 1.0/atan2(B,A)/2.0;
	if(phi<0)phi+=2.0*M_PI;
	if(phi>=2.0*M_PI)phi-=2.0*M_PI;
#endif

	return NOERROR;
}

//-----------------
// SearchPtrans
//-----------------
void DHelicalFit::SearchPtrans(double ptrans_max, double ptrans_step)
{
	/// Search the chi2 space as a function of the transverse momentum
	/// for the minimal chi2. The values corresponding to the minmal
	/// chi2 are left in chisq, x0, y0, r0, q, and p_trans.
	///
	/// This does NOT optimize the value of p_trans. It simply
	/// does a straight forward chisq calculation on a grid
	/// and keeps the best one. It is intended for use in finding
	/// a reasonable value for p_trans for straight tracks that
	/// are contained to less than p_trans_max which is presumably
	/// chosen based on a known &theta;.

	// Loop over several values for p_trans and calculate the
	// chisq for each.
	double min_chisq=1.0E6;
	double min_x0=0.0, min_y0=0.0, min_r0=0.0;
	for(double my_p_trans=ptrans_step; my_p_trans<=ptrans_max; my_p_trans+=ptrans_step){

		r0 = my_p_trans/(0.003*2.0);
		double alpha, my_chisq;

		// q = +1
		alpha = phi + M_PI_2;
		x0 = r0*cos(alpha);
		y0 = r0*sin(alpha);
		my_chisq = ChisqCircle();
		if(my_chisq<min_chisq){
			min_chisq=my_chisq;
			min_x0 = x0;
			min_y0 = y0;
			min_r0 = r0;
			q = +1.0;
			p_trans = my_p_trans;
		}

		// q = -1
		alpha = phi - M_PI_2;
		x0 = r0*cos(alpha);
		y0 = r0*sin(alpha);
		my_chisq = ChisqCircle();
		if(my_chisq<min_chisq){
			min_chisq=my_chisq;
			min_x0 = x0;
			min_y0 = y0;
			min_r0 = r0;
			q = -1.0;
			p_trans = my_p_trans;
		}
	}
	
	// Copy params from minimum chisq
	x0 = min_x0;
	y0 = min_y0;
	r0 = min_r0;
}

//-----------------
// QuickPtrans
//-----------------
void DHelicalFit::QuickPtrans(void)
{
	/// Quickly calculate a value of p_trans by looking
	/// for the hit furthest out and the hit closest
	/// to half that distance. Those 2 hits along with the
	/// origin are used to define a circle from which
	/// p_trans is calculated.
	
	// Find hit with largest R
	double R2max = 0.0;
	DHFHit_t *hit_max = NULL;
	for(unsigned int i=0; i<hits.size(); i++){
		// use cross product to decide if hits is to left or right
		double x = hits[i]->x;
		double y = hits[i]->y;
		double R2 = x*x + y*y;
		if(R2>R2max){
			R2max = R2;
			hit_max = hits[i];
		}
	}
	
	// Bullet proof
	if(!hit_max)return;
	
	// Find hit closest to half-way out
	double Rmid = 0.0;
	double Rmax = sqrt(R2max);
	DHFHit_t *hit_mid = NULL;
	for(unsigned int i=0; i<hits.size(); i++){
		// use cross product to decide if hits is to left or right
		double x = hits[i]->x;
		double y = hits[i]->y;
		double R = sqrt(x*x + y*y);
		if(fabs(R-Rmax/2.0) < fabs(Rmid-Rmax/2.0)){
			Rmid = R;
			hit_mid = hits[i];
		}
	}
	
	// Bullet proof
	if(!hit_mid)return;
	
	// Calculate p_trans from 3 points
	double x1 = hit_mid->x;
	double y1 = hit_mid->y;
	double x2 = hit_max->x;
	double y2 = hit_max->y;
	double r2 = sqrt(x2*x2 + y2*y2);
	double cos_phi = x2/r2;
	double sin_phi = y2/r2;
	double u1 =  x1*cos_phi + y1*sin_phi;
	double v1 = -x1*sin_phi + y1*cos_phi;
	double u2 =  x2*cos_phi + y2*sin_phi;
	double u0 = u2/2.0;
	double v0 = (u1*u1 + v1*v1 - u1*u2)/(2.0*v1);
	
	x0 = u0*cos_phi - v0*sin_phi;
	y0 = u0*sin_phi + v0*cos_phi;
	r0 = sqrt(x0*x0 + y0*y0);
	
	double B=-2.0;
	p_trans = fabs(0.003*B*r0);
	q = v1>0.0 ? -1.0:+1.0;
}

//-----------------
// GuessChargeFromCircleFit
//-----------------
jerror_t DHelicalFit::GuessChargeFromCircleFit(void)
{
	/// Adjust the sign of the charge (magnitude will stay 1) based on
	/// whether the hits tend to be to the right or to the left of the
	/// line going from the origin through the center of the circle.
	/// If the sign is flipped, the phi angle will also be shifted by
	/// +/- pi since the majority of hits are assumed to come from
	/// outgoing tracks.
	///
	/// This is just a guess since tracks can bend all the way
	/// around and have hits that look exactly like tracks from an
	/// outgoing particle of opposite charge. The final charge should
	/// come from the sign of the dphi/dz slope.
	
	// Simply count the number of hits whose phi angle relative to the
	// phi of the center of the circle are greater than pi.
	unsigned int N=0;
	for(unsigned int i=0; i<hits.size(); i++){
		// use cross product to decide if hits is to left or right
		double x = hits[i]->x;
		double y = hits[i]->y;
		if((x*y0 - y*x0) < 0.0)N++;
	}
	
	// Check if more hits are negative and make sign negative if so.
	if(N>hits.size()/2.0){
		q = -1.0;
		phi += M_PI;
		if(phi>2.0*M_PI)phi-=2.0*M_PI;
	}
	
	return NOERROR;
}

//-----------------
// FitTrack
//-----------------
jerror_t DHelicalFit::FitTrack(void)
{
	/// Find theta, sign of electric charge, total momentum and
	/// vertex z position.

	// Points must be in order of increasing Z
	sort(hits.begin(), hits.end(), DHFHitLessThanZ_C);

	// Fit to circle to get circle's center
	FitCircle();
	
	// The thing that is really needed is dphi/dz (where phi is the angle
	// of the point as measured from the center of the circle, not the beam
	// line). The relation between phi and z is linear so we use linear
	// regression to find the slope (dphi/dz). The one complication is
	// that phi is periodic so the value obtained via the x and y of a
	// specific point can be off by an integral number of 2pis. Handle this
	// by assuming the first point is measured before a full rotation
	// has occurred. Subsequent points should always be within 2pi of
	// the previous point so we just need to calculate the relative phi
	// between succesive points and keep a running sum. We do this in
	// the first loop were we find the mean z and phi values. The regression
	// formulae are calculated in the second loop.

	// Calculate phi about circle center for each hit
	Fill_phi_circle(hits, x0, y0);

	// Linear regression for z/phi relation
	float Sxx=0.0, Syy=0.0, Sxy=0.0;
	for(unsigned int i=0;i<hits.size();i++){
		DHFHit_t *a = hits[i];
		float deltaZ = a->z - z_mean;
		float deltaPhi = a->phi_circle - phi_mean;
		Syy += deltaZ*deltaZ;
		Sxx += deltaPhi*deltaPhi;
		Sxy += deltaZ*deltaPhi;
		//cout<<__FILE__<<":"<<__LINE__<<" deltaZ="<<deltaZ<<" deltaPhi="<<deltaPhi<<" Sxy(i)="<<deltaZ*deltaPhi<<endl;
	}
	float dzdphi = Syy/Sxy;
	dzdphi = Syy/Sxy;
	z_vertex = z_mean - phi_mean*dzdphi;
//cout<<__FILE__<<":"<<__LINE__<<" z_mean="<<z_mean<<" phi_mean="<<phi_mean<<" dphidz="<<dphidz<<" Sxy="<<Sxy<<" Syy="<<Syy<<" z_vertex="<<z_vertex<<endl;
	
	// Fill in the rest of the parameters
	return FillTrackParams();
}

//-----------------
// FitTrackRiemann
//-----------------
jerror_t DHelicalFit::FitTrackRiemann(float rc_input){
  jerror_t error=NOERROR;

  if (CovR_!=NULL) delete CovR_;
  if (CovRPhi_!=NULL) delete CovRPhi_; 
  CovR_=NULL;
  CovRPhi_=NULL;

  error=FitCircleRiemannCorrected(rc_input);
  error=FitLineRiemann();
  GetChargeRiemann();
  
  // Shift phi by pi if the charge is negative
  if (q<0) phi+=M_PI; 

  return error;
}



//-----------------
// FitTrack_FixedZvertex
//-----------------
jerror_t DHelicalFit::FitTrack_FixedZvertex(float z_vertex)
{
	/// Fit the points, but hold the z_vertex fixed at the specified value.
	///
	/// This just calls FitCircle and FitLine_FixedZvertex to find all parameters
	/// of the track

	// Fit to circle to get circle's center
	FitCircle();
	
	// Fit to line in phi-z plane and return error
	return FitLine_FixedZvertex(z_vertex);
}

//-----------------
// FitLine_FixedZvertex
//-----------------
jerror_t DHelicalFit::FitLine_FixedZvertex(float z_vertex)
{
	/// Fit the points, but hold the z_vertex fixed at the specified value.
	///
	/// This assumes FitCircle has already been called and the values
	/// in x0 and y0 are valid.
	///
	/// This just fits the phi-z angle by minimizing the chi-squared
	/// using a linear regression technique. As it turns out, the
	/// chi-squared weights points by their distances squared which
	/// leads to a quadratic equation for cot(theta) (where theta is 
	/// the angle in the phi-z plane). To pick the right solution of
	/// the quadratic equation, we pick the one closest to the linearly
	/// weighted angle. One could argue that we should just use the
	/// linear weighting here rather than the square weighting. The
	/// choice depends though on how likely the "out-lier" points
	/// are to really be from this track. If they are likely, to
	/// be, then we would want to leverage their "longer" lever arms
	/// with the squared weighting. If they are more likely to be bad
	/// points, then we would want to minimize their influence with
	/// a linear (or maybe even root) weighting. It is expected than
	/// for our use, the points are all likely valid so we use the
	/// square weighting.
	
	 // Points must be in order of increasing Z
	sort(hits.begin(), hits.end(), DHFHitLessThanZ_C);
	
	// Fit is being done for a fixed Z-vertex
	this->z_vertex = z_vertex;

	// Calculate phi about circle center for each hit
	Fill_phi_circle(hits, x0, y0);

	// Do linear regression on phi-Z
	float Sx=0, Sy=0;
	float Sxx=0, Syy=0, Sxy = 0;
	float r0 = sqrt(x0*x0 + y0*y0);
	for(unsigned int i=0; i<hits.size(); i++){
		DHFHit_t *a = hits[i];

		float dz = a->z - z_vertex;
		float dphi = a->phi_circle;
		Sx  += dz;
		Sy  += r0*dphi;
		Syy += r0*dphi*r0*dphi;
		Sxx += dz*dz;
		Sxy += r0*dphi*dz;
	}
	
	float k = (Syy-Sxx)/(2.0*Sxy);
	float s = sqrt(k*k + 1.0);
	float cot_theta1 = -k+s;
	float cot_theta2 = -k-s;
	float cot_theta_lin = Sx/Sy;
	float cot_theta;
	if(fabs(cot_theta_lin-cot_theta1) < fabs(cot_theta_lin-cot_theta2)){
		cot_theta = cot_theta1;
	}else{
		cot_theta = cot_theta2;
	}
	
	dzdphi = r0*cot_theta;
	
	// Fill in the rest of the paramters
	return FillTrackParams();
}

//------------------------------------------------------------------
// Fill_phi_circle
//------------------------------------------------------------------
jerror_t DHelicalFit::Fill_phi_circle(vector<DHFHit_t*> hits, float x0, float y0)
{
	float x_last = -x0;
	float y_last = -y0;
	float phi_last = 0.0;
	z_mean = phi_mean = 0.0;
	for(unsigned int i=0; i<hits.size(); i++){
		DHFHit_t *a = hits[i];

		float dx = a->x - x0;
		float dy = a->y - y0;
		float dphi = atan2f(dx*y_last - dy*x_last, dx*x_last + dy*y_last);
		float my_phi = phi_last +dphi;
		a->phi_circle = my_phi;

		z_mean += a->z;
		phi_mean += my_phi;
		
		x_last = dx;
		y_last = dy;
		phi_last = my_phi;
	}	
	z_mean /= (float)hits.size();
	phi_mean /= (float)hits.size();

	return NOERROR;
}

//------------------------------------------------------------------
// FillTrackParams
//------------------------------------------------------------------
jerror_t DHelicalFit::FillTrackParams(void)
{
	/// Fill in and tweak some parameters like q, phi, theta using
	/// other values already set in the class. This is used by
	/// both FitTrack() and FitTrack_FixedZvertex(). 

	float r0 = sqrt(x0*x0 + y0*y0);
	theta = atan(r0/fabs(dzdphi));
	
	// The sign of the electric charge will be opposite that
	// of dphi/dz. Also, the value of phi will be PI out of phase
	if(dzdphi<0.0){
		phi += M_PI;
		if(phi<0)phi+=2.0*M_PI;
		if(phi>=2.0*M_PI)phi-=2.0*M_PI;
	}else{
		q = -q;
	}
	
	// Re-calculate chi-sq (needs to be done)
	chisq_source = TRACK;
	
	// Up to now, the fit has assumed a forward going particle. In
	// other words, if the particle is going backwards, the helix does
	// actually still go through the points, but only if extended
	// backward in time. We use the average z of the hits compared
	// to the reconstructed vertex to determine if it was back-scattered.
	if(z_mean<z_vertex){
		// back scattered particle
		theta = M_PI - theta;
		phi += M_PI;
		if(phi<0)phi+=2.0*M_PI;
		if(phi>=2.0*M_PI)phi-=2.0*M_PI;
		q = -q;
	}

	// There is a problem with theta for data generated with a uniform
	// field. The following factor is emprical until I can figure out
	// what the source of the descrepancy is.
	//theta += 0.1*pow((double)sin(theta),2.0);
	
	p = fabs(p_trans/sin(theta));
	
	return NOERROR;
}

//------------------------------------------------------------------
// Print
//------------------------------------------------------------------
jerror_t DHelicalFit::Print(void) const
{
	cout<<"-- DHelicalFit Params ---------------"<<endl;
	cout<<"          x0 = "<<x0<<endl;
	cout<<"          y0 = "<<y0<<endl;
	cout<<"           q = "<<q<<endl;
	cout<<"           p = "<<p<<endl;
	cout<<"     p_trans = "<<p_trans<<endl;
	cout<<"         phi = "<<phi<<endl;
	cout<<"       theta = "<<theta<<endl;
	cout<<"    z_vertex = "<<z_vertex<<endl;
	cout<<"       chisq = "<<chisq<<endl;
	cout<<"chisq_source = ";
	switch(chisq_source){
		case NOFIT:		cout<<"NOFIT";		break;
		case CIRCLE:	cout<<"CIRCLE";	break;
		case TRACK:		cout<<"TRACK";		break;
	}
	cout<<endl;
	cout<<"     Bz(avg) = "<<Bz_avg<<endl;

	return NOERROR;
}


//------------------------------------------------------------------
// Dump
//------------------------------------------------------------------
jerror_t DHelicalFit::Dump(void) const
{
	Print();

	for(unsigned int i=0;i<hits.size();i++){
		DHFHit_t *v = hits[i];
		cout<<" x="<<v->x<<" y="<<v->y<<" z="<<v->z;
		cout<<" phi_circle="<<v->phi_circle<<" chisq="<<v->chisq<<endl;
	}
	
	return NOERROR;
}
