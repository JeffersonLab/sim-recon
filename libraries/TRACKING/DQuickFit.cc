// $Id$
// Author: David Lawrence  June 25, 2004
//
//
//

#include <iostream>
using namespace std;

#include <math.h>

#include "DQuickFit.h"
#include "DMagneticFieldMap.h"


// The following is for sorting hits by z
class DQFHitLessThanZ{
	public:
		bool operator()(DQFHit_t* const &a, DQFHit_t* const &b) const {
			return a->z < b->z;
		}
};

//-----------------
// DQuickFit
//-----------------
DQuickFit::DQuickFit()
{	
	x0 = y0 = 0;
	chisq = 0;
	chisq_source = NOFIT;
	bfield = NULL;

}

//-----------------
// ~DQuickFit
//-----------------
DQuickFit::~DQuickFit()
{
	for(unsigned int i=0; i<hits.size(); i++)delete hits[i];
	hits.clear();
}

//-----------------
// AddHit
//-----------------
derror_t DQuickFit::AddHit(float r, float phi, float z)
{
	/// Add a hit to the list of hits using cylindrical coordinates

	return AddHitXYZ(r*cos(phi), r*sin(phi), z);
}

//-----------------
// AddHitXYZ
//-----------------
derror_t DQuickFit::AddHitXYZ(float x, float y, float z)
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
derror_t DQuickFit::PruneHit(int idx)
{
	/// Remove the hit specified by idx from the list
	/// of hits. The value of idx can be anywhere from
	/// 0 to GetNhits()-1.
	if(idx<0 || idx>=(int)hits.size())return VALUE_OUT_OF_RANGE;

	hits.erase(hits.begin() + idx);

	return NOERROR;
}


//-----------------
// PrintChiSqVector
//-----------------
derror_t DQuickFit::PrintChiSqVector(void)
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
derror_t DQuickFit::FitCircle(void)
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
	///    (this also implies no multiple-scattering)
	/// 2. The vertex is at the target center (i.e. 0,0,0 in the coordinate
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
	
	// Loop over hits to calculate alpha, beta, gamma, and delta
	// if a magnetic field map was given, use it to find average Z B-field
	float Bz = 0.0;
	for(unsigned int i=0;i<hits.size();i++){
		DQFHit_t *a = hits[i];
		float x=a->x;
		float y=a->y;
		alpha += x*x;
		beta += y*y;
		gamma += x*y;
		deltax += x*(x*x+y*y)/2.0;
		deltay += y*(x*x+y*y)/2.0;
		
		if(bfield){
			D3Vector_t tmp;
			tmp = bfield->getQuick(a->x/2.54, a->y/2.54, 26.0+a->z/2.54);
			Bz += tmp.z;
		}
	}
	
	// Calculate x0,y0 - the center of the circle
	x0 = (deltax*beta-deltay*gamma)/(alpha*beta-gamma*gamma);
	y0 = (deltay-gamma*x0)/beta;

	// Momentum depends on magnetic field. Assume 2T for now.
	// Also assume a singly charged track (i.e. q=+/-1)
	// The sign of the charge will be determined below.
	Bz_avg=-2.0; 
	if(bfield)Bz_avg = Bz/(float)hits.size();
	q = +1.0;
	float r0 = sqrt(x0*x0 + y0*y0);
	p_trans = q*Bz_avg*r0*qBr2p; // qBr2p converts to GeV/c
	phi = atan2(y0,x0) + M_PI_2;
	if(p_trans<0.0){
		p_trans = -p_trans;
	}
	if(phi<0)phi+=2.0*M_PI;
	if(phi>=2.0*M_PI)phi-=2.0*M_PI;
	
	// Calculate the chisq
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
	chisq_source = CIRCLE;

	return NOERROR;
}

//-----------------
// FitTrack
//-----------------
derror_t DQuickFit::FitTrack(void)
{
	/// Find theta, sign of electric charge, total momentum and
	/// vertex z position.

	// Points must be in order of increasing Z
	sort(hits.begin(), hits.end(), DQFHitLessThanZ());

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
	dphidz = Sxy/Syy;
	z_vertex = z_mean - phi_mean*dzdphi;
//cout<<__FILE__<<":"<<__LINE__<<" z_mean="<<z_mean<<" phi_mean="<<phi_mean<<" dphidz="<<dphidz<<" Sxy="<<Sxy<<" Syy="<<Syy<<" z_vertex="<<z_vertex<<endl;
	
	// Fill in the rest of the paramters
	return FillTrackParams();
}

//-----------------
// FitTrack_FixedZvertex
//-----------------
derror_t DQuickFit::FitTrack_FixedZvertex(float z_vertex)
{
	/// Fit the points, but hold the z_vertex fixed at the specified value.
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
	sort(hits.begin(), hits.end(), DQFHitLessThanZ());
	
	// Fit is being done for a fixed Z-vertex
	this->z_vertex = z_vertex;

	// Fit to circle to get circle's center
	FitCircle();
	
	// Calculate phi about circle center for each hit
	Fill_phi_circle(hits, x0, y0);

	// Do linear regression on phi-Z
	float Sx=0, Sy=0;
	float Sxx=0, Syy=0, Sxy = 0;
	for(unsigned int i=0; i<hits.size(); i++){
		DQFHit_t *a = hits[i];

		float dz = a->z - z_vertex;
		float dphi = a->phi_circle;
		Sx  += dz;
		Sy  += dphi;
		Syy += dphi*dphi;
		Sxx += dz*dz;
		Sxy += dphi*dz;
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
	
	//dphidz = 1.0/cot_theta;
	dphidz = Sy/Sx;
	
	// Fill in the rest of the paramters
	return FillTrackParams();
}

//------------------------------------------------------------------
// Fill_phi_circle
//------------------------------------------------------------------
derror_t DQuickFit::Fill_phi_circle(vector<DQFHit_t*> hits, float x0, float y0)
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
//cout<<__FILE__<<":"<<__LINE__<<" a.x="<<a.x<<" a.y="<<a.y<<" dphi="<<dphi<<" my_phi="<<my_phi<<endl;
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
derror_t DQuickFit::FillTrackParams(void)
{
	/// Fill in and tweak some parameters like q, phi, theta using
	/// other values already set in the class. This is used by
	/// both FitTrack() and FitTrack_FixedZvertex(). 

	float r0 = sqrt(x0*x0 + y0*y0);
	theta = atan(r0*fabs(dphidz));
	p = fabs(p_trans/sin(theta));
//cout<<__FILE__<<":"<<__LINE__<<" theta="<<theta<<" dphidz="<<dphidz<<" p_trans="<<p_trans<<" p="<<p<<" r0="<<r0<<endl;
	// The sign of the electric charge will be opposite that
	// of dphi/dz. Also, the value of phi will be PI out of phase
	if(dphidz<0.0){
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
	// actually go still go through the points, but only if extended
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
	
	return NOERROR;
}

//------------------------------------------------------------------
// Print
//------------------------------------------------------------------
derror_t DQuickFit::Print(void)
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
derror_t DQuickFit::Dump(void)
{
	Print();

	for(unsigned int i=0;i<hits.size();i++){
		DQFHit_t *v = hits[i];
		cout<<" x="<<v->x<<" y="<<v->y<<" z="<<v->z;
		cout<<" phi_circle="<<v->phi_circle<<" chisq="<<chisq<<endl;
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
derror_t DQuickFit::AddHits(int N, TVector3 *v)
{
	/// Append a list of hits to the current list of hits using
	/// TVector3 objects. The TVector3 objects are copied internally
	/// so it is safe to delete the objects after calling AddHits()
	/// For 2D hits, the value of z will be ignored.

	for(int i=0; i<N; i++, v++){
		TVector3 *vec = new TVector3(*v);
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
derror_t DQuickFit::PruneHits(float chisq_limit)
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
derror_t DQuickFit::PruneOutlier(void)
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
		TVector3 *v = hits[i];
		X += v->x();
		Y += v->y();
	}
	X /= (float)hits.size();
	Y /= (float)hits.size();
	
	float max =0.0;
	int idx = -1;
	for(unsigned int i=0;i<hits.size();i++){
		TVector3 *v = hits[i];
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
derror_t DQuickFit::PruneOutliers(int n)
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
derror_t DCDC::firstguess_curtis(s_Cdc_trackhit_t *hits, int Npoints
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
	float dphidz =0.0;
	int Ndphidzpoints = 0;
	for(int i=0;i<Npoints-1;i++, v++, v2++){
		float myphi1 = atan2(v->y-y0,  v->x-x0);
		float myphi2 = atan2(v2->y-y0, v2->x-x0);
		float mydphidz = (myphi2-myphi1)/(v2->z-v->z);
		if(finite(mydphidz)){
			dphidz+=mydphidz;
			Ndphidzpoints++;
		}
	}
	if(Ndphidzpoints){
		dphidz/=(float)Ndphidzpoints;
	}
	
	theta = atan(r0*fabs(dphidz));
	p = -p_trans/sin(theta);

	return NOERROR;
}

#endif

