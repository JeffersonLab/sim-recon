// Author: David Lawrence  June 25, 2004
//
//
// DCDC event processor methods
//

#include <iostream>
using namespace std;

#include "DCDC.h"

static int qsort_cdc(const void* arg1, const void* arg2);


//------------------------------------------------------------------
// evnt   -Reconstruction of single event
//------------------------------------------------------------------
derror_t DCDC::evnt(int eventnumber)
{
	// Convert CDC hits to simpler format
	copy_to_cdchits();
	
	// Allocate memory for cdc_tracks. (there can't be more tracks than hits)
	hddm->cdc_tracks = make_s_Cdc_tracks(hddm->cdc_trackhits->mult);
	
	// Fit CDC tracks
	s_Cdc_trackhits_t *cdc_trackhits = hddm->cdc_trackhits;
	int track_start = 0;
	for(int j=1;j<cdc_trackhits->mult;j++){
		if(cdc_trackhits->in[j].track!=cdc_trackhits->in[j-1].track || j==cdc_trackhits->mult-1){
			// For now, we just use the "first guess" parameters until
			// a real fitting routine is written
			float theta,phi,p,q;
			firstguess(&cdc_trackhits->in[track_start], j-track_start
				,theta, phi, p, q);
			track_start = j;
			
			// Fill in cdctrack
			s_Cdc_track_t *cdctrack = &hddm->cdc_tracks->in[hddm->cdc_tracks->mult++];
			cdctrack->theta = theta;
			cdctrack->phi = phi;
			cdctrack->p = p;
			cdctrack->q = q;
			cdctrack->px = p*sin(theta)*cos(phi);
			cdctrack->py = p*sin(theta)*sin(phi);
			cdctrack->pz = p*cos(theta);
			cdctrack->track = cdc_trackhits->in[j-1].track;
		}
	}

	return NOERROR;
}

//------------------------------------------------------------------
// copy_to_cdchits 
//------------------------------------------------------------------
derror_t DCDC::copy_to_cdchits(void)
{
	/// Copy from from s_CdcPoints_t to s_Cdc_trackhits_t
	
	// We need to allocate memory for the s_Cdc_trackhits_t, but
	// we don't know how many hits there are and it just doesn't seem
	// worth going through "N" levels to get the exact number. Just
	// allocate 1000.
	hddm->cdc_trackhits = make_s_Cdc_trackhits(1000);
	int Ntrackhits = hddm->cdc_trackhits->mult = 0;
	s_Cdc_trackhit_t *trackhits = hddm->cdc_trackhits->in;

	// Loop over Physics Events
	s_PhysicsEvents_t* PE = hddm->physicsEvents;
	if(!PE)return NOERROR;
	for(int i=0; i<PE->mult; i++){
	
		s_Rings_t *rings=NULL;
		s_HitView_t *HV = PE->in[i].hitView;
		if(PE->in[i].hitView)
			if(PE->in[i].hitView->centralDC)
				rings = PE->in[i].hitView->centralDC->rings;
		if(rings){
			for(int j=0;j<rings->mult;j++){
				float radius = rings->in[j].radius;
				s_Straws_t *straws = rings->in[j].straws;
				if(straws){
					for(int k=0;k<straws->mult;k++){
						float phim = straws->in[k].phim;
						s_CdcPoints_t *cdcPoints = straws->in[k].cdcPoints;
						if(cdcPoints){						
							for(int m=0;m<cdcPoints->mult;m++){
								float	r = cdcPoints->in[m].r;
								float phi = cdcPoints->in[m].phi;
								float z = cdcPoints->in[m].z;
								if(Ntrackhits>=1000)continue;
								
								trackhits[Ntrackhits].x = r*cos(phi);
								trackhits[Ntrackhits].y = r*sin(phi);
								trackhits[Ntrackhits].z = z;
								trackhits[Ntrackhits].t = 0.0;
								trackhits[Ntrackhits].track = cdcPoints->in[m].track;
								Ntrackhits++;
							}
						}
					}
				}
			}
		}
	}
	
	hddm->cdc_trackhits->mult = Ntrackhits;
	
	// Sort by track number, then z-position (ascending for both)
	qsort(hddm->cdc_trackhits->in, Ntrackhits, sizeof(s_Cdc_trackhit_t), qsort_cdc);

	return NOERROR;
}

//------------------------------------------------------------------
// qsort_cdc_trackhits
//------------------------------------------------------------------
static int qsort_cdc(const void* arg1, const void* arg2)
{
	s_Cdc_trackhit_t *a = (s_Cdc_trackhit_t*)arg1;
	s_Cdc_trackhit_t *b = (s_Cdc_trackhit_t*)arg2;
	
	// Sort by track first ...
	if(a->track != b->track)return a->track - b->track;
	
	// ... then z
	if(a->z == b->z)return 0;
	return b->z > a->z ? 1:-1;
}

//-----------------
// firstguess
//-----------------
derror_t DCDC::firstguess(s_Cdc_trackhit_t *hits, int Npoints
	, float &theta ,float &phi, float &p, float &q)
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
	s_Cdc_trackhit_t *v=hits;
	for(int i=0;i<Npoints;i++, v++){
		float x=v->x;
		float y=v->y;
		alpha += x*x;
		beta += y*y;
		gamma += x*y;
		deltax += x*(x*x+y*y)/2.0;
		deltay += y*(x*x+y*y)/2.0;
	}
	
	// Calculate x0,y0 - the center of the circle
	float x0 = (deltax*beta-deltay*gamma)/(alpha*beta-gamma*gamma);
	float y0 = (deltay-gamma*x0)/beta;

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
	v = hits;
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
