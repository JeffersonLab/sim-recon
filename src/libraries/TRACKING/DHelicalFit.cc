// Set of routines for fitting tracks in both the cdc and fdc assuming a
// uniform magnetic field.  The track therefore is assumed to follow a helical
// path.

#include <iostream>
#include <algorithm>
#include <cmath>
using namespace std;

#include "DANA/DApplication.h"
#include <TDecompLU.h>
#include <math.h>

#include "DHelicalFit.h"
#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c

#define ONE_THIRD  0.33333333333333333
#define SQRT3      1.73205080756887719
#define EPS 1e-8
#ifndef M_TWO_PI
#define M_TWO_PI 6.28318530717958647692
#endif

// The following is for sorting hits by z
class DHFHitLessThanZ{
	public:
		bool operator()(DHFHit_t* const &a, DHFHit_t* const &b) const {
			return a->z < b->z;
		}
};

inline bool DHFHitLessThanZ_C(DHFHit_t* const &a, DHFHit_t* const &b) {
	return a->z < b->z;
}
inline bool RiemannFit_hit_cmp(DHFHit_t *a,DHFHit_t *b){
  return (a->z<b->z);
}

bool DHFProjection_cmp(const DHFProjection_t &a,
		       const DHFProjection_t &b){
  if(fabs(a.z - b.z) > 0.0)
    return (a.z<b.z);
  return (a.xy.Mod2() > b.xy.Mod2());
}


//-----------------
// DHelicalFit (Constructor)
//-----------------
DHelicalFit::DHelicalFit(void)
{
  x0 = 0., y0 = 0., r0=0., tanl=0., z_vertex = 0., phi=0., theta = 0.,dzdphi=0.;
  h=1.;

  chisq = 0;
  ndof=0;
  chisq_source = NOFIT;
  bfield = NULL;
  
  // For Riemann Fit
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

void DHelicalFit::Reset(void)
{
	x0 = 0.0;
	y0 = 0.0;
	r0 = 0.0;
	h = 1.0;
	phi = 0.0;
	theta = 0.0;
	tanl = 0.0;
	z_vertex = 0.0;
	chisq = 0.0;
	ndof = 0;
	dzdphi = 0.0;
	chisq_source = NOFIT;
 
	normal.SetXYZ(0.,0.,0.);
	c_origin = 0.0;
	
	Clear();
	//hits.clear();
	z_mean = 0.0;
	phi_mean = 0.0;
  
	N[0] = 0.0;
	N[1] = 0.0;
	N[2] = 0.0;

	xavg[0] = 0.0;
	xavg[1] = 0.0;
	xavg[2] = 0.0;

	var_avg = 0.0;
}

//-----------------
// Copy
//-----------------
void DHelicalFit::Copy(const DHelicalFit &fit)
{
	x0 = fit.x0;
	y0 = fit.y0;
	r0 = fit.r0;
	h = fit.h;
	tanl=fit.tanl;
	phi = fit.phi;
	theta = fit.theta;
	z_vertex = fit.z_vertex;
	chisq = fit.chisq;
	ndof=fit.ndof;
	dzdphi = fit.dzdphi;
	chisq_source = fit.chisq_source;
	z_mean = fit.GetZMean();
	phi_mean = fit.GetPhiMean();
	
	const vector<DHFHit_t*> myhits = fit.GetHits();
	for(unsigned int i=0; i<myhits.size(); i++){
		DHFHit_t *a = new DHFHit_t;
		*a = *myhits[i];
		hits.push_back(a);
	}
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
  hit->x = fdchit->xy.X();
  hit->y = fdchit->xy.Y();
  hit->z = fdchit->wire->origin.z();
  hit->covx=fdchit->covxx;
  hit->covy=fdchit->covyy;
  hit->covxy=fdchit->covxy;
  hit->chisq = 0.0;
  hit->is_axial=false;

  double Phi=atan2(hit->y,hit->x);

  // Covariance matrix elements
  double u=fdchit->w;
  double v=fdchit->s;
  double temp1=u*Phi-v;
  double temp2=v*Phi+u;
  double var_u=0.08333; // 1.0 cm^2/12  
  double var_v=0.0075;// 0.3*0.3 cm^2/12 
  double one_over_R2=1./fdchit->xy.Mod2();
  hit->covrphi=one_over_R2*(var_v*temp1*temp1+var_u*temp2*temp2);
  hit->covr=one_over_R2*(var_u*u*u+var_v*v*v);

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
	hit->is_axial=false;
	
  double Phi=atan2(hit->y,hit->x);
  double cosPhi=cos(Phi);
  double sinPhi=sin(Phi);
  double Phi_cosPhi_minus_sinPhi=Phi*cosPhi-sinPhi;
  double Phi_sinPhi_plus_cosPhi=Phi*sinPhi+cosPhi;
  hit->covrphi=Phi_cosPhi_minus_sinPhi*Phi_cosPhi_minus_sinPhi*hit->covx
    +Phi_sinPhi_plus_cosPhi*Phi_sinPhi_plus_cosPhi*hit->covy
    +2.*Phi_sinPhi_plus_cosPhi*Phi_cosPhi_minus_sinPhi*hit->covxy;
  hit->covr=cosPhi*cosPhi*hit->covx+sinPhi*sinPhi*hit->covy
    +2.*sinPhi*cosPhi*hit->covxy;

	hits.push_back(hit);

	return NOERROR;
}

// Add a hit to the list of hits using Cartesian coordinates
jerror_t DHelicalFit::AddHitXYZ(float x,float y, float z,float covx,
				float covy, float covxy, bool is_axial){
  DHFHit_t *hit = new DHFHit_t;
  hit->x = x;
  hit->y = y;
  hit->z = z;
  hit->covx=covx;
  hit->covy=covy;
  hit->covxy=covxy;
  hit->is_axial=is_axial;

  double Phi=atan2(hit->y,hit->x);
  double cosPhi=cos(Phi);
  double sinPhi=sin(Phi);
  double Phi_cosPhi_minus_sinPhi=Phi*cosPhi-sinPhi;
  double Phi_sinPhi_plus_cosPhi=Phi*sinPhi+cosPhi;
  hit->covrphi=Phi_cosPhi_minus_sinPhi*Phi_cosPhi_minus_sinPhi*hit->covx
    +Phi_sinPhi_plus_cosPhi*Phi_sinPhi_plus_cosPhi*hit->covy
    +2.*Phi_sinPhi_plus_cosPhi*Phi_cosPhi_minus_sinPhi*hit->covxy;
  hit->covr=cosPhi*cosPhi*hit->covx+sinPhi*sinPhi*hit->covy
    +2.*sinPhi*cosPhi*hit->covxy;

  hits.push_back(hit);

  return NOERROR;
}

// Routine to add stereo hits once an initial circle fit has been done
jerror_t DHelicalFit::AddStereoHit(const DCDCWire *wire){
  DVector3 origin = wire->origin;
  DVector3 dir = wire->udir;
  double dx=origin.x()-x0;  
  double dy=origin.y()-y0;
  double ux=dir.x();
  double uy=dir.y();
  double temp1=ux*ux+uy*uy;
  double temp2=ux*dy-uy*dx;
  double b=-ux*dx-uy*dy;
  double r0_sq=r0*r0;
  double A=r0_sq*temp1-temp2*temp2;

  // Check that this wire intersects this circle
  if(A<0.0) return VALUE_OUT_OF_RANGE; // line along wire does not intersect circle, ever.

  // Calculate intersection points for the two roots 
  double B = sqrt(A);
  double dz1 = (b-B)/temp1;
  double dz2 = (b+B)/temp1;
  
  // At this point we must decide which value of alpha to use. 
  // For now, we just use the value closest to zero (i.e. closest to
  // the center of the wire).
  double dz=dz1;
  if (fabs(dz2)<fabs(dz1)){
    dz=dz2;
  }
  // distance along wire relative to origin
  double s=dz/cos(wire->stereo);
  
  //  if(DEBUG_LEVEL>15)
    _DBG_<<"s="<<s<<" ring="<<wire->ring<<" straw="<<wire->straw<<" stereo="<<wire->stereo<<endl;
  if(fabs(s) > 0.5*wire->L) return VALUE_OUT_OF_RANGE; // if wire doesn't cross circle, skip hit
		
  // Compute the position for this hit
  DVector3 pos = origin + s*dir;
 
  // Assume error on the radius of the circle is 10% of the radius and include
  // the cell size contribution 1.6*1.6/12
  double var_r=0.01*r0_sq+0.213; 

  DHFHit_t *hit = new DHFHit_t;
  hit->x=pos.x();
  hit->y=pos.y();
  hit->z=pos.z();
  hit->covx=0.213; // place holder
  hit->covy=0.213;
  hit->covxy=0.;
  hit->chisq=0.;
  hit->is_axial=false;

  double Phi=atan2(hit->y,hit->x);
  double cosPhi=cos(Phi);
  double sinPhi=sin(Phi);
  double Phi_cosPhi_minus_sinPhi=Phi*cosPhi-sinPhi;
  double Phi_sinPhi_plus_cosPhi=Phi*sinPhi+cosPhi;
  hit->covrphi=Phi_cosPhi_minus_sinPhi*Phi_cosPhi_minus_sinPhi*hit->covx
    +Phi_sinPhi_plus_cosPhi*Phi_sinPhi_plus_cosPhi*hit->covy;
  hit->covr=var_r;

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
  // Remove all hits
  for(unsigned int i=0; i<hits.size(); i++)delete hits[i];
  hits.clear();
 
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

	float alpha=0.0, beta=0.0, gamma=0.0, deltax=0.0, deltay=0.0;
	chisq_source = NOFIT; // in case we return early
	size_t num_hits=hits.size();
	ndof=num_hits-2;
	
	// Loop over hits to calculate alpha, beta, gamma, and delta
	// if a magnetic field map was given, use it to find average Z B-field
	DHFHit_t *a = NULL;
	for(unsigned int i=0;i<num_hits;i++){
	  a = hits[i];
	  float x=a->x;
	  float y=a->y;
	  float x_sq=x*x;
	  float y_sq=y*y;
	  float x_sq_plus_y_sq=x_sq+y_sq;
	  alpha += x_sq;
	  beta += y_sq;
	  gamma += x*y;
	  deltax += 0.5*x*x_sq_plus_y_sq;
	  deltay += 0.5*y*x_sq_plus_y_sq;
	}
	
	// Calculate x0,y0 - the center of the circle
	double denom = alpha*beta-gamma*gamma;
	if(fabs(denom)<1.0E-20)return UNRECOVERABLE_ERROR;
	x0 = (deltax*beta-deltay*gamma)/denom;
	y0 = (deltay*alpha-deltax*gamma)/denom;
	r0 = sqrt(x0*x0 + y0*y0);

	phi = atan2(y0,x0) - M_PI_2;
	if(phi<0)phi+=M_TWO_PI;
	if(phi>=M_TWO_PI)phi-=M_TWO_PI;
	
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
	/// raw chi2 
	chisq = 0.0;
	for(unsigned int i=0;i<hits.size();i++){
		DHFHit_t *a = hits[i];
		float x = a->x - x0;
		float y = a->y - y0;
		float x2=x*x,y2=y*y;
		float r2=x2+y2;
		float c = sqrt(r2) - r0;
		c *= c;
		double cov=(x2*a->covx+y2*a->covy+2.*x*y*a->covxy)/r2;
		c/=cov;
		a->chisq = c;
		chisq+=c;
	}
	
	// Do NOT divide by DOF
	return chisq;
}

//-----------------
// FitCircleRiemann
//-----------------
// Arguments are for inserting a "vertex" point
jerror_t DHelicalFit::FitCircleRiemann(float z_vertex,float BeamRMS)
{
	chisq_source = NOFIT; // in case we return early
	
	// Fake point at origin
	float beam_var=BeamRMS*BeamRMS;
	AddHitXYZ(0.,0.,z_vertex,beam_var,beam_var,0.);

	jerror_t err = FitCircleRiemann();
	if(err!=NOERROR)return err;

	// Number of degrees of freedom 
	ndof=hits.size()-3;

	phi=atan2(-x0,y0);  
	if(phi<0)phi+=M_TWO_PI;
        if(phi>=M_TWO_PI)phi-=M_TWO_PI;
	
	// Normal vector for plane intersecting Riemann surface 
	normal.SetXYZ(N[0],N[1],N[2]);

	return NOERROR;
}

//-----------------
// FitCircleRiemann
//-----------------
jerror_t DHelicalFit::FitCircleRiemann(float rc){  
/// Riemann Circle fit:  points on a circle in x,y project onto a plane cutting
/// the circular paraboloid surface described by (x,y,x^2+y^2).  Therefore the
/// task of fitting points in (x,y) to a circle is transormed to the task of
/// fitting planes in (x,y, w=x^2+y^2) space
///
  size_t num_hits=hits.size();
  DMatrix X(num_hits,3);
  DMatrix Xavg(1,3);
  DMatrix A(3,3);
  // vector of ones
  DMatrix OnesT(1,num_hits);
  double W_sum=0.;
  DMatrix W(num_hits,num_hits);

  // Make sure hit list is ordered in z
  std::sort(hits.begin(),hits.end(),RiemannFit_hit_cmp);
 
  // Covariance matrix
  DMatrix CRPhi(num_hits,num_hits);
  for (unsigned int i=0;i<num_hits;i++){
    CRPhi(i,i)=hits[i]->covrphi;
  }
  
  // Apply a correction for non-normal track incidence if we already have a 
  // guess for rc.
  if (rc>0.){
    DMatrix C(num_hits,num_hits);
    DMatrix S(num_hits,num_hits);
    DMatrix CR(num_hits,num_hits);
    for (unsigned int i=0;i<num_hits;i++){
      S(i,i)=0.;
      C(i,i)=1.;
      CR(i,i)=hits[i]->covr;

      double rtemp=sqrt(hits[i]->x*hits[i]->x+hits[i]->y*hits[i]->y); 
      double stemp=rtemp/(4.*rc);
      double ctemp=1.-stemp*stemp;
      if (ctemp>0){
	S(i,i)=stemp;
	C(i,i)=sqrt(ctemp);
      }
    }
    // Correction for non-normal incidence of track on FDC 
    CRPhi=C*CRPhi*C+S*CR*S;
  }

  // The goal is to find the eigenvector corresponding to the smallest 
  // eigenvalue of the equation
  //            lambda=n^T (X^T W X - W_sum Xavg^T Xavg)n
  // where n is the normal vector to the plane slicing the cylindrical 
  // paraboloid described by the parameterization (x,y,w=x^2+y^2),
  // and W is the weight matrix, assumed for now to be diagonal.
  // In the absence of multiple scattering, W_sum is the sum of all the 
  // diagonal elements in W.

  for (unsigned int i=0;i<num_hits;i++){
    X(i,0)=hits[i]->x;
    X(i,1)=hits[i]->y;
    X(i,2)=hits[i]->x*hits[i]->x+hits[i]->y*hits[i]->y;
    OnesT(0,i)=1.;
    W(i,i)=1./CRPhi(i,i);
    W_sum+=W(i,i);
  }
  Xavg=(1./W_sum)*(OnesT*(W*X));

  A=DMatrix(DMatrix::kTransposed,X)*(W*X)
    -W_sum*(DMatrix(DMatrix::kTransposed,Xavg)*Xavg);
  if(!A.IsValid())return UNRECOVERABLE_ERROR;

  // The characteristic equation is 
  //   lambda^3+B2*lambda^2+lambda*B1+B0=0 
  //
  double B2=-(A(0,0)+A(1,1)+A(2,2));
  double B1=A(0,0)*A(1,1)-A(1,0)*A(0,1)+A(0,0)*A(2,2)-A(2,0)*A(0,2)
    +A(1,1)*A(2,2)-A(2,1)*A(1,2);
  double B0=-A.Determinant();
  if(B0==0 || !isfinite(B0))return UNRECOVERABLE_ERROR;

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

  long double Q=(3.*B1-B2*B2)/9.e6;//4; 
  long double R=(9.*B2*B1-27.*B0-2.*B2*B2*B2)/54.e9;//6;
  long double Q1=Q*Q*Q+R*R;
  if (Q1<0) Q1=sqrtl(-Q1);
  else{
    return VALUE_OUT_OF_RANGE;
  }

  // DeMoivre's theorem for fractional powers of complex numbers:  
  //      (r*(cos(theta)+i sin(theta)))^(p/q)
  //                  = r^(p/q)*(cos(p*theta/q)+i sin(p*theta/q))
  //
  //double temp=100.*pow(R*R+Q1*Q1,0.16666666666666666667);
  long double temp=1000.*sqrtl(cbrtl(R*R+Q1*Q1));
  long double theta1=ONE_THIRD*atan2l(Q1,R);
  long double sum_over_2=temp*cosl(theta1);
  long double diff_over_2=-temp*sinl(theta1);
  // Third root
  long double lambda_min=-ONE_THIRD*B2-sum_over_2+SQRT3*diff_over_2;
  
  // Calculate the (normal) eigenvector corresponding to the eigenvalue lambda
  double A11_minus_lambda_min=A(1,1)-lambda_min;
  N[0]=1.;
  N[1]=(A(1,0)*A(0,2)-(A(0,0)-lambda_min)*A(1,2))
    /(A(0,1)*A(2,1)-A11_minus_lambda_min*A(0,2));
  N[2]=(A(2,0)*A11_minus_lambda_min-A(1,0)*A(2,1))
    /(A(1,2)*A(2,1)-(A(2,2)-lambda_min)*A11_minus_lambda_min);
  
  // Normalize: n1^2+n2^2+n3^2=1
  double denom=sqrt(N[0]*N[0]+N[1]*N[1]+N[2]*N[2]);
  for (int i=0;i<3;i++){
    N[i]/=denom;
  }

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
  if(phi<0)phi+=M_TWO_PI;
  if(phi>=M_TWO_PI)phi-=M_TWO_PI;

  // Calculate the chisq
  ChisqCircle();
  chisq_source = CIRCLE;
  ndof=hits.size()-3;

  return NOERROR;
}

//-----------------
// FitLineRiemann
//-----------------
jerror_t DHelicalFit::FitLineRiemann(){
/// Riemann Line fit: linear regression of s on z to determine the tangent of 
/// the dip angle and the z position of the closest approach to the beam line.
/// Computes intersection points along the helical path.
///   Note: this implementation assumes that the error in the z-position is 
/// negligible; practically speaking, this means it should only be used for FDC
/// hits...
  size_t num_hits=hits.size();
 
  // Fill vector of intersection points 
  double denom= N[0]*N[0]+N[1]*N[1];
  // projection vector
  vector<DHFProjection_t>projections;
  for (unsigned int m=0;m<num_hits;m++){
    if (hits[m]->is_axial) continue;
  
    double r2=hits[m]->x*hits[m]->x+hits[m]->y*hits[m]->y;
    double numer=c_origin+r2*N[2];
    double ratio=numer/denom;
    double x_int0=-N[0]*ratio;
    double y_int0=-N[1]*ratio;
    double temp=denom*r2-numer*numer;
    if (temp<0){  // Skip point if the intersection gives nonsense
      //_DBG_ << "Bad?" <<endl;
      continue;
    }
    temp=sqrt(temp)/denom;
  
    // Store projection data in a temporary structure
    DHFProjection_t temp_proj;
    temp_proj.z=hits[m]->z;
    temp_proj.covR=hits[m]->covr;
    
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
      temp_proj.xy=DVector2(x2,y2);
    }
    else{
      temp_proj.xy=DVector2(x1,y1);
    }
    projections.push_back(temp_proj);
  }  
  if (projections.size()==0) return RESOURCE_UNAVAILABLE;
  sort(projections.begin(),projections.end(),DHFProjection_cmp);

  // Linear regression to find z0, tanl   
  unsigned int n=projections.size();
  double sumv=0.,sumx=0.,sumy=0.,sumxx=0.,sumxy=0.;
  double sperp=0.,sperp_old=0., ratio=0, Delta;
  double z_last=0.,z=0.;
  DVector2 old_proj=projections[0].xy;
  double two_r0=2.*r0;
  for (unsigned int k=0;k<n;k++){
    sperp_old=sperp;
    z_last=z;
    double chord=(projections[k].xy-old_proj).Mod();
    ratio=chord/two_r0;

    // Make sure the argument for the arcsin does not go out of range...
    double ds=(ratio>1)? two_r0*M_PI_2 : two_r0*asin(ratio);
    sperp=sperp_old+ds;
    z=projections[k].z;
    
    double weight=1./projections[k].covR;
    sumv+=weight;
    sumy+=sperp*weight;
    sumx+=z*weight;
    sumxx+=z*z*weight;
    sumxy+=sperp*z*weight;
    
    // Store the current x and y projection values
    old_proj=projections[k].xy;
  }
  Delta=sumv*sumxx-sumx*sumx;
  double tanl_denom=sumv*sumxy-sumy*sumx;
  if (fabs(Delta)<EPS || fabs(tanl_denom)<EPS) return VALUE_OUT_OF_RANGE;

  // Track parameters tan(lambda) and z-vertex
  tanl=Delta/tanl_denom; 
 
  double chord=projections[0].xy.Mod();
  ratio=chord/two_r0;
  sperp=(ratio>1? two_r0*M_PI_2 : two_r0*asin(ratio));
  z_vertex=projections[0].z-sperp*tanl;

  /*
  if (z_vertex<Z_MIN){
    z_vertex=Z_MIN;
    tanl=(z_last-Z_MIN)/sperp;
  }
  else if (z_vertex>Z_MAX){
    z_vertex=Z_MAX;
    tanl=(z_last-Z_MAX)/sperp;
  }
  */
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
	size_t num_hits=hits.size();
	for(unsigned int i=0;i<num_hits;i++){
		a = hits[i];
		double r = sqrt(pow((double)a->x,2.0) + pow((double)a->y, 2.0));
		// weight by r to give outer points more influence. Note that
		// we really are really weighting by r^2 since x and y already
		// have a magnitude component.
		X += a->x*r; 
		Y += a->y*r;
	}
	phi = atan2(Y,X);
	if(phi<0)phi+=M_TWO_PI;
	if(phi>=M_TWO_PI)phi-=M_TWO_PI;

	// Search the chi2 space for values for x0, y0, ...
	SearchPtrans(9.0, 0.5);

#if 0
	// We do a simple linear regression here to find phi. This is
	// simplified by the intercept being zero (i.e. the track
	// passes through the beamline).
	float Sxx=0.0, Syy=0.0, Sxy=0.0;
	chisq_source = NOFIT; // in case we return early
	
	// Loop over hits to calculate Sxx, Syy, and Sxy
	DHFHit_t *a = NULL;
	for(unsigned int i=0;i<num_hits;i++){
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
	if(phi<0)phi+=M_TWO_PI;
	if(phi>=M_TWO_PI)phi-=M_TWO_PI;
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

		// h = +1
		alpha = phi + M_PI_2;
		x0 = r0*cos(alpha);
		y0 = r0*sin(alpha);
		my_chisq = ChisqCircle();
		if(my_chisq<min_chisq){
			min_chisq=my_chisq;
			min_x0 = x0;
			min_y0 = y0;
			min_r0 = r0;
			h = +1.0;
		}

		// h = -1
		alpha = phi - M_PI_2;
		x0 = r0*cos(alpha);
		y0 = r0*sin(alpha);
		my_chisq = ChisqCircle();
		if(my_chisq<min_chisq){
			min_chisq=my_chisq;
			min_x0 = x0;
			min_y0 = y0;
			min_r0 = r0;
			h = -1.0;
		}
	}
	
	// Copy params from minimum chisq
	x0 = min_x0;
	y0 = min_y0;
	r0 = min_r0;
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
	if(N>hits.size()/2){
		h = -1.0;
		phi += M_PI;
		if(phi>M_TWO_PI)phi-=M_TWO_PI;
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
	z_vertex = z_mean - phi_mean*dzdphi;
//cout<<__FILE__<<":"<<__LINE__<<" z_mean="<<z_mean<<" phi_mean="<<phi_mean<<" dphidz="<<dphidz<<" Sxy="<<Sxy<<" Syy="<<Syy<<" z_vertex="<<z_vertex<<endl;
	
	// Fill in the rest of the parameters
	return FillTrackParams();
}

//-----------------
// FitTrackRiemann
//-----------------
jerror_t DHelicalFit::FitTrackRiemann(float rc_input){
  jerror_t error=FitCircleRiemann(rc_input); 
  if (error!=NOERROR) return error;
  error=FitLineRiemann();
  FindSenseOfRotation();

  return error;
}


jerror_t DHelicalFit::FitCircleAndLineRiemann(float rc_input){
  jerror_t error=FitCircleRiemann(rc_input);
  if (error!=NOERROR) return error;
  error=FitLineRiemann();

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
	
	// the sense of rotation about the z axis will be opposite that
	// of dphi/dz. Also, the value of phi will be PI out of phase
	if(dzdphi<0.0){
		phi += M_PI;
		if(phi<0)phi+=M_TWO_PI;
		if(phi>=M_TWO_PI)phi-=M_TWO_PI;
	}else{
		h = -h;
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
		if(phi<0)phi+=M_TWO_PI;
		if(phi>=M_TWO_PI)phi-=M_TWO_PI;
		h = -h;
	}
	tanl=tan(M_PI_2-theta);

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
	cout<<"           h = "<<h<<endl;
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

// Find the sense of the rotation about the circle based on the proximity of 
// the hits to the helical trajectory for the two hypotheses 
void DHelicalFit::FindSenseOfRotation(void){
  // Try to use a plane with a single hit for the reference position
  unsigned int i=1;
  for (;i<hits.size()-1;i++){
    if (fabs(hits[i-1]->z-hits[i]->z)>0.01
	&& fabs(hits[i]->z-hits[i+1]->z)>0.01) break;
  }

  double Phi1=atan2(hits[i]->y-y0,hits[i]->x-x0);
  double z0=hits[i]->z;
  double plus_sum=0.,minus_sum=0.;
  for (unsigned int j=0;j<hits.size();j++){
    DHFHit_t *hit = hits[j];
    double dphi=(hit->z-z0)/(r0*tanl);

    double phiplus=Phi1+dphi;
    double dxplus=x0+r0*cos(phiplus)-hit->x;
    double dyplus=y0+r0*sin(phiplus)-hit->y;
    double dxplus2=dxplus*dxplus;
    double dyplus2=dyplus*dyplus;
    double d2plus=dxplus2+dyplus2;
    double varplus=(dxplus2*hit->covx+dyplus2*hit->covy
		    +2.*dyplus*dxplus*hit->covxy)/d2plus;
    plus_sum+=d2plus/varplus;
    
    double phiminus=Phi1-dphi;
    double dxminus=x0+r0*cos(phiminus)-hit->x;
    double dyminus=y0+r0*sin(phiminus)-hit->y;
    double dxminus2=dxminus*dxminus;
    double dyminus2=dyminus*dyminus;
    double d2minus=dxminus2+dyminus2;
    double varminus=(dxminus2*hit->covx+dyminus2*hit->covy
		     +2.*dyminus*dxminus*hit->covxy)/d2minus;
    minus_sum+=d2minus/varminus;

    
    //printf("z %f d+ %f %f d- %f %f\n",hit->z,d2plus,plus_sum,d2minus,minus_sum);

  }
  
  h=1.0;
  // Look for smallest sum to determine q
  if (minus_sum<plus_sum) h=-1.;
}
