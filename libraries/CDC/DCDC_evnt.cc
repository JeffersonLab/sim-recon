// Author: David Lawrence  June 25, 2004
//
//
// DCDC event processor methods
//

#include <iostream>
using namespace std;

#include "DCDC.h"

static int qsort_cdchit(const void* arg1, const void* arg2);


//------------------------------------------------------------------
// evnt   -Reconstruction of single event
//------------------------------------------------------------------
derror_t DCDC::evnt(int eventnumber)
{

	// Sort CDC hits by track, then z
	qsort(hddm->CDChits->CDChit, hddm->CDChits->nrows, sizeof(CDChit_t), qsort_cdchit);
		
	// Fit CDC tracks
	CDChit_t *CDChit = hddm->CDChits->CDChit;
	int track_start = 0;
	for(int j=0;j<hddm->CDChits->nrows-1;j++, CDChit++){
		
		if(CDChit->track!=CDChit[1].track || j==hddm->CDChits->nrows-2){
			// Do a first guess calculation of track parameters
			float theta,phi,p,q,x0,y0;
			firstguess(&hddm->CDChits->CDChit[track_start] ,j+1-track_start ,theta, phi, p, q, x0, y0);
			track_start = j+1;
			
			// Fill in CDCtrack
			CDCtracks_t *CDCtracks = hddm->CDCtracks;
			CDCtracks->nrows++;
			CDCtracks->Grow();
			CDCtrack_t *CDCtrack = &CDCtracks->CDCtrack[CDCtracks->nrows-1];

			CDCtrack->q = q;
			CDCtrack->dir.SetXYZ(sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta));
			CDCtrack->p.SetVect(p*CDCtrack->dir);
			CDCtrack->track = CDChit->track;
			CDCtrack->x0 = x0;
			CDCtrack->y0 = y0;
		}
	}

	return NOERROR;
}

//------------------------------------------------------------------
// qsort_cdc_trackhits
//------------------------------------------------------------------
static int qsort_cdchit(const void* arg1, const void* arg2)
{
	CDChit_t *a = (CDChit_t*)arg1;
	CDChit_t *b = (CDChit_t*)arg2;
	
	// Sort by track first ...
	if(a->track != b->track)return a->track - b->track;
	
	// ... then z
	if(a->pos.z() == b->pos.z())return 0;
	return b->pos.z() > a->pos.z() ? 1:-1;
}

//-----------------
// firstguess
//-----------------
derror_t DCDC::firstguess(CDChit_t *hits, int Npoints
	, float &theta ,float &phi, float &p, float &q, float &x0, float &y0)
{
	/// This will determine starting parameters for the fit of a CDC track.
	/// The alogorithm used here calculates the parameters directly using
	/// a technique very much like linear regression. The key assumptions
	/// are:
	/// 1. The magnetic field is uniform and along z so that the projection
	///    of the track onto the X/Y plane will fall on a circle
	/// 2. The vertex is at the target center (i.e. 0,0,0 in the coordinate
	///    system of the TVector3 points passed to us.

	float alpha=0.0, beta=0.0, gamma=0.0, deltax=0.0, deltay=0.0;
	
	// Loop over hits to calculate alpha, beta, gamma, and delta
	CDChit_t *v=hits;
	for(int i=0;i<Npoints;i++, v++){
		float x=v->pos.x();
		float y=v->pos.y();
		alpha += x*x;
		beta += y*y;
		gamma += x*y;
		deltax += x*(x*x+y*y)/2.0;
		deltay += y*(x*x+y*y)/2.0;
	}
	
	// Calculate x0,y0 - the center of the circle
	x0 = (deltax*beta-deltay*gamma)/(alpha*beta-gamma*gamma);
	y0 = (deltay-gamma*x0)/beta;

	// Momentum depends on magnetic field. Assume 2T for now.
	// Also assume a singly charged track (i.e. q=+/-1)
	// The sign of the charge will be deterined below.
	float B=-2.0*0.593; // The 0.5931 is empirical
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
	v = hits;
	CDChit_t *v2=&hits[n_2];
	for(int i=0;i<n_2;i++, v++, v2++){
		xprod_sum += v->pos.x()*v2->pos.y() - v2->pos.x()*v->pos.y();
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
		float myphi1 = atan2(v->pos.y()-y0,  v->pos.x()-x0);
		float myphi2 = atan2(v2->pos.y()-y0, v2->pos.x()-x0);
		float mydphidz = (myphi2-myphi1)/(v2->pos.z()-v->pos.z());
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


#if 0
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
