// $Id$
// Author: David Lawrence  June 25, 2004
//
//
//

#include <iostream>
#include <algorithm>
using namespace std;

#include <math.h>

#include "DQuickFit.h"
#include "DRiemannFit.h"
#define qBr2p 0.003  // conversion for converting q*B*r to GeV/c
#define Z_VERTEX 65.0

// The following is for sorting hits by z
class DQFHitLessThanZ{
	public:
		bool operator()(DQFHit_t* const &a, DQFHit_t* const &b) const {
			return a->z < b->z;
		}
};

bool DQFHitLessThanZ_C(DQFHit_t* const &a, DQFHit_t* const &b) {
	return a->z < b->z;
}

//-----------------
// DQuickFit (Constructor)
//-----------------
DQuickFit::DQuickFit(void)
{
	x0 = y0 = 0;
	chisq = 0;
	chisq_source = NOFIT;
	bfield = NULL;

	c_origin=0.;
	normal.SetXYZ(0.,0.,0.);
}

//-----------------
// DQuickFit (Constructor)
//-----------------
DQuickFit::DQuickFit(const DQuickFit &fit)
{
	Copy(fit);
}
//-----------------
// Copy
//-----------------
void DQuickFit::Copy(const DQuickFit &fit)
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
	
	const vector<DQFHit_t*> myhits = fit.GetHits();
	for(unsigned int i=0; i<myhits.size(); i++){
		DQFHit_t *a = new DQFHit_t;
		*a = *myhits[i];
		hits.push_back(a);
	}
}

//-----------------
// operator= (Assignment operator)
//-----------------
DQuickFit& DQuickFit::operator=(const DQuickFit& fit)
{
	if(this == &fit)return *this;
	Copy(fit);

	return *this;
}

//-----------------
// DQuickFit (Destructor)
//-----------------
DQuickFit::~DQuickFit()
{
	Clear();
}

//-----------------
// AddHit
//-----------------
jerror_t DQuickFit::AddHit(float r, float phi, float z)
{
	/// Add a hit to the list of hits using cylindrical coordinates

	return AddHitXYZ(r*cos(phi), r*sin(phi), z);
}

//-----------------
// AddHitXYZ
//-----------------
jerror_t DQuickFit::AddHitXYZ(float x, float y, float z)
{
	/// Add a hit to the list of hits useing Cartesian coordinates
	DQFHit_t *hit = new DQFHit_t;
	hit->x = x;
	hit->y = y;
	hit->z = z;
	hit->chisq = 0.0;
	hits.push_back(hit);

	return NOERROR;
}

//-----------------
// PruneHit
//-----------------
jerror_t DQuickFit::PruneHit(int idx)
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
jerror_t DQuickFit::Clear(void)
{
	/// Remove all hits
	for(unsigned int i=0; i<hits.size(); i++)delete hits[i];
	hits.clear();

	return NOERROR;
}

//-----------------
// PrintChiSqVector
//-----------------
jerror_t DQuickFit::PrintChiSqVector(void) const
{
	/// Dump the latest chi-squared vector to the screen.
	/// This prints the individual hits' chi-squared
	/// contributions in the order in which the hits were
	/// added. See PruneHits() for more detail.

	cout<<"Chisq vector from DQuickFit: (source=";
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
jerror_t DQuickFit::FitCircle(void)
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
	chisq_source = NOFIT; // in case we reutrn early
	
	// Loop over hits to calculate alpha, beta, gamma, and delta
	// if a magnetic field map was given, use it to find average Z B-field
	DQFHit_t *a = NULL;
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
double DQuickFit::ChisqCircle(void)
{
	/// Calculate the chisq for the hits and current circle
	/// parameters.
	/// NOTE: This does not return the chi2/dof, just the
	/// raw chi2 with errs set to 1.0
	chisq = 0.0;
	for(unsigned int i=0;i<hits.size();i++){
		DQFHit_t *a = hits[i];
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
jerror_t DQuickFit::FitCircleRiemann(double BeamRMS)
{
	/// This is a temporary solution for doing a Riemann circle fit.
	/// It uses Simon's DRiemannFit class. This is not very efficient
	/// since it creates a new DRiemann fit object, copies the hits
	/// into it, and then copies the reults back to this DQuickFit
	/// object every time this is called. At some point, the DRiemannFit
	/// class will be merged into DQuickFit.

	chisq_source = NOFIT; // in case we reutrn early

	DRiemannFit rfit;
	for(unsigned int i=0; i<hits.size(); i++){
		DQFHit_t *hit = hits[i];
		rfit.AddHitXYZ(hit->x, hit->y, hit->z);
	}
	
	// Fake point at origin
	double beam_var=BeamRMS*BeamRMS;
	rfit.AddHit(0.,0.,Z_VERTEX,beam_var,beam_var,0.);

	jerror_t err = rfit.FitCircle(BeamRMS);
	if(err!=NOERROR)return err;
	
	x0 = rfit.xc;
	y0 = rfit.yc;
	phi = rfit.phi;

	//r0 = sqrt(x0*x0+y0*y0);
	r0=rfit.rc;

	rfit.GetPlaneParameters(c_origin,normal);
	
	// Momentum depends on magnetic field. If bfield has been
	// set, we should use it to determine an average value of Bz
	// for this track. Otherwise, assume -2T.
	// Also assume a singly charged track (i.e. q=+/-1)
	// The sign of the charge will be determined below.
	Bz_avg=-2.0; 
	q = +1.0;
	float r0 = sqrt(x0*x0 + y0*y0);
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
// FitCircleStraightTrack
//-----------------
jerror_t DQuickFit::FitCircleStraightTrack(void)
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
	DQFHit_t *a = NULL;
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
	chisq_source = NOFIT; // in case we reutrn early
	
	// Loop over hits to calculate Sxx, Syy, and Sxy
	DQFHit_t *a = NULL;
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
void DQuickFit::SearchPtrans(double ptrans_max, double ptrans_step)
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
void DQuickFit::QuickPtrans(void)
{
	/// Quickly calculate a value of p_trans by looking
	/// for the hit furthest out and the hit closest
	/// to half that distance. Those 2 hits along with the
	/// origin are used to define a circle from which
	/// p_trans is calculated.
	
	// Find hit with largest R
	double R2max = 0.0;
	DQFHit_t *hit_max = NULL;
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
	DQFHit_t *hit_mid = NULL;
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
jerror_t DQuickFit::GuessChargeFromCircleFit(void)
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
jerror_t DQuickFit::FitTrack(void)
{
	/// Find theta, sign of electric charge, total momentum and
	/// vertex z position.

	// Points must be in order of increasing Z
	sort(hits.begin(), hits.end(), DQFHitLessThanZ_C);

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
		DQFHit_t *a = hits[i];
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
// FitTrack_FixedZvertex
//-----------------
jerror_t DQuickFit::FitTrack_FixedZvertex(float z_vertex)
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
jerror_t DQuickFit::FitLine_FixedZvertex(float z_vertex)
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
	sort(hits.begin(), hits.end(), DQFHitLessThanZ_C);
	
	// Fit is being done for a fixed Z-vertex
	this->z_vertex = z_vertex;

	// Calculate phi about circle center for each hit
	Fill_phi_circle(hits, x0, y0);

	// Do linear regression on phi-Z
	float Sx=0, Sy=0;
	float Sxx=0, Syy=0, Sxy = 0;
	float r0 = sqrt(x0*x0 + y0*y0);
	for(unsigned int i=0; i<hits.size(); i++){
		DQFHit_t *a = hits[i];

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
jerror_t DQuickFit::Fill_phi_circle(vector<DQFHit_t*> hits, float x0, float y0)
{
	float x_last = -x0;
	float y_last = -y0;
	float phi_last = 0.0;
	z_mean = phi_mean = 0.0;
	for(unsigned int i=0; i<hits.size(); i++){
		DQFHit_t *a = hits[i];

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
jerror_t DQuickFit::FillTrackParams(void)
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
jerror_t DQuickFit::Print(void) const
{
	cout<<"-- DQuickFit Params ---------------"<<endl;
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
jerror_t DQuickFit::Dump(void) const
{
	Print();

	for(unsigned int i=0;i<hits.size();i++){
		DQFHit_t *v = hits[i];
		cout<<" x="<<v->x<<" y="<<v->y<<" z="<<v->z;
		cout<<" phi_circle="<<v->phi_circle<<" chisq="<<v->chisq<<endl;
	}
	
	return NOERROR;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------
//--------------------------- UNUSED --------------------------------
//-------------------------------------------------------------------
//-------------------------------------------------------------------


#if 0

//-----------------
// AddHits
//-----------------
jerror_t DQuickFit::AddHits(int N, DVector3 *v)
{
	/// Append a list of hits to the current list of hits using
	/// DVector3 objects. The DVector3 objects are copied internally
	/// so it is safe to delete the objects after calling AddHits()
	/// For 2D hits, the value of z will be ignored.

	for(int i=0; i<N; i++, v++){
		DVector3 *vec = new DVector3(*v);
		if(vec)
			hits.push_back(vec);
		else
			cerr<<__FILE__<<":"<<__LINE__<<" NULL vector in DQuickFit::AddHits. Hit dropped."<<endl;
	}

	return NOERROR;
}

//-----------------
// PruneHits
//-----------------
jerror_t DQuickFit::PruneHits(float chisq_limit)
{
	/// Remove hits whose individual contribution to the chi-squared
	/// value exceeds <i>chisq_limit</i>. The value of the individual
	/// contribution is calculated like:
	///
	///	\f$ r_0 = \sqrt{x_0^2 + y_0^2} \f$ <br>
	///	\f$ r_i = \sqrt{(x_i-x_0)^2 + (y_i-y_0)^2} \f$ <br>
	///	\f$ \chi_i^2 = (r_i-r_0)^2 \f$ <br>

	if(hits.size() != chisqv.size()){
		cerr<<__FILE__<<":"<<__LINE__<<" hits and chisqv do not have the same number of rows!"<<endl;
		cerr<<"Call FitCircle() or FitTrack() method first!"<<endl;

		return NOERROR;
	}
	
	// Loop over these backwards to make it easier to
	// handle the deletes
	for(int i=chisqv.size()-1; i>=0; i--){
		if(chisqv[i] > chisq_limit)PruneHit(i);
	}
	
	return NOERROR;
}

//-----------------
// PruneOutliers
//-----------------
jerror_t DQuickFit::PruneOutlier(void)
{
	/// Remove the point which is furthest from the geometric
	/// center of all the points in the X/Y plane.
	///
	/// This just calculates the average x and y values of
	/// registered hits. It finds the distance of every hit
	/// with respect to the geometric mean and removes the
	/// hit whose furthest from the mean.
	
	float X=0, Y=0;
	for(unsigned int i=0;i<hits.size();i++){
		DVector3 *v = hits[i];
		X += v->x();
		Y += v->y();
	}
	X /= (float)hits.size();
	Y /= (float)hits.size();
	
	float max =0.0;
	int idx = -1;
	for(unsigned int i=0;i<hits.size();i++){
		DVector3 *v = hits[i];
		float x = v->x()-X;
		float y = v->y()-Y;
		float dist_sq = x*x + y*y; // we don't need to take sqrt just to find max
		if(dist_sq>max){
			max = dist_sq;
			idx = i;
		}
	}
	if(idx>=0)PruneHit(idx);
	
	return NOERROR;
}

//-----------------
// PruneOutliers
//-----------------
jerror_t DQuickFit::PruneOutliers(int n)
{
	/// Remove the n points which are furthest from the geometric
	/// center of all the points in the X/Y plane.
	///
	/// This just calls PruneOutlier() n times. Since the mean
	/// can change each time, it should be recalculated after
	/// every hit is removed.
	
	for(int i=0;i<n;i++)PruneOutlier();
	
	return NOERROR;
}

//-----------------
// firstguess_curtis
//-----------------
jerror_t DCDC::firstguess_curtis(s_Cdc_trackhit_t *hits, int Npoints
	, float &theta ,float &phi, float &p, float &q)
{
	if(Npoints<3)return NOERROR;
	// pick out 3 points to calculate the circle with
	s_Cdc_trackhit_t *hit1, *hit2, *hit3;
	hit1 = hits;
	hit2 = &hits[Npoints/2];
	hit3 = &hits[Npoints-1];

	float x1 = hit1->x, y1=hit1->y;
	float x2 = hit2->x, y2=hit2->y;
	float x3 = hit3->x, y3=hit3->y;

	float b = (x2*x2+y2*y2-x1*x1-y1*y1)/(2.0*(x2-x1));
	b -= (x3*x3+y3*y3-x1*x1-y1*y1)/(2.0*(x3-x1));
	b /= ((y1-y2)/(x1-x2)) - ((y1-y3)/(x1-x3));
	float a = (x2*x2-y2*y2-x1*x1-y1*y1)/(2.0*(x2-x1));
	a -= b*(y2-y1)/(x2-x1);

	// Above is the method from Curtis's crystal barrel note 93, pg 72
	// Below here is just a copy of the code from David's firstguess
	// routine above (after the x0,y0 calculation)
	float x0=a, y0=b;

	// Momentum depends on magnetic field. Assume 2T for now.
	// Also assume a singly charged track (i.e. q=+/-1)
	// The sign of the charge will be deterined below.
	float B=-2.0*0.593; // The 0.61 is empirical
	q = +1.0;
	float r0 = sqrt(x0*x0 + y0*y0);
	float hbarc = 197.326;
	float p_trans = q*B*r0/hbarc; // are these the right units?
	
	// Assuming points are ordered in increasing z, the sign of the
	// cross-product between successive points will be the opposite
	// sign of the charge. Since it's possible to have a few "bad"
	// points, we don't want to rely on any one to determine this.
	// The method we use is to sum cross-products of the first and
	// middle points, the next-to-first and next-to-middle, etc.
	float xprod_sum = 0.0;
	int n_2 = Npoints/2; 
	s_Cdc_trackhit_t *v = hits;
	s_Cdc_trackhit_t *v2=&hits[n_2];
	for(int i=0;i<n_2;i++, v++, v2++){
		xprod_sum += v->x*v2->y - v2->x*v->y;
	}
	if(xprod_sum>0.0)q = -q;

	// Phi is pi/2 out of phase with x0,y0. The sign of the phase difference
	// depends on the charge
	phi = atan2(y0,x0);
	phi += q>0.0 ? -M_PI_2:M_PI_2;
	
	// Theta is determined by extrapolating the helix back to the target.
	// To do this, we need dphi/dz and a phi,z point. The easiest way to
	// get these is by a simple average (I guess).
	v = hits;
	v2=&v[1];
	float dzdphi =0.0;
	int Ndzdphipoints = 0;
	for(int i=0;i<Npoints-1;i++, v++, v2++){
		float myphi1 = atan2(v->y-y0,  v->x-x0);
		float myphi2 = atan2(v2->y-y0, v2->x-x0);
		float mydzdphi = (v2->z-v->z)/(myphi2-myphi1);
		if(finite(mydzdphi)){
			dzdphi+=mydzdphi;
			Ndzdphipoints++;
		}
	}
	if(Ndzdphipoints){
		dzdphi/=(float)Ndzdphipoints;
	}
	
	theta = atan(r0/fabs(dzdphi));
	p = -p_trans/sin(theta);

	return NOERROR;
}

#endif

