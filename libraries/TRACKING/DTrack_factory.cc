// $Id$
//
//    File: DTrack_factory.cc
// Created: Sun Apr  3 12:28:45 EDT 2005
// Creator: davidl (on Darwin Harriet.local 7.8.0 powerpc)
//

#include <math.h>

#include <TVector3.h>
#include <JANA/JEventLoop.h>

#include "GlueX.h"
#include "DANA/DApplication.h"
#include "DMagneticFieldStepper.h"
#include "DTrackCandidate.h"
#include "DTrack_factory.h"
#include "DTrackHit.h"
#include "DReferenceTrajectory.h"

//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(JEventLoop *loop, int runnumber)
{
	//dgeom = loop->GetJApplication()->GetJGeometry(runnumber);
	//bfield = NULL;
	//bfield = dgeom->GetDMagneticFieldMap();
	
	bfield = new DMagneticFieldMap(); // temporary until new geometry scheme is worked out
	
	MAX_HIT_DIST = 10.0; // cm
	
	dparms.SetDefaultParameter("TRK:MAX_HIT_DIST",	MAX_HIT_DIST);
	
	dparms.GetParameter("TRK:TRACKHIT_SOURCE",	TRACKHIT_SOURCE);

	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get the track candidates and hits
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DTrackHit*> trackhits;
	loop->Get(trackcandidates);
	loop->Get(trackhits,TRACKHIT_SOURCE.c_str());

	// Loop over track candidates
	for(unsigned int i=0; i<trackcandidates.size(); i++){

		// Fit the track		
		DTrack *track = FitTrack(trackcandidates[i], trackhits);
		
		// If fit is successful, then store the track
		if(track)_data.push_back(track);
	}

	return NOERROR;
}

//------------------
// FitTrack
//------------------
DTrack* DTrack_factory::FitTrack(const DTrackCandidate *trackcandidate, vector<const DTrackHit*> trackhits)
{
	// Generate reference trajectory for this track
	DReferenceTrajectory *rt = new DReferenceTrajectory(bfield, trackcandidate);
		
	// Determine the distance of each hit to the reference trajectory
	// This is not an efficient way of doing this, but it will work for
	// now and should be optimized later.
	for(unsigned int j=0; j<trackhits.size(); j++){
		const DTrackHit *trackhit = trackhits[j];
		double min_delta2 = 1.0E6;
		int min_step_index = -1;

		// Loop over swim steps. Don't allow the first or last points 
		// to be considered as the closest so we're guaranteed to
		// have points on either side of the closest swim point.
		swim_step_t *swim_step = rt->swim_steps;
		swim_step++;
		for(int i=1; i<rt->Nswim_steps-1; i++, swim_step++){
			double deltaZ = swim_step->pos.z() - trackhit->z;
			if(fabs(deltaZ) > MAX_HIT_DIST)continue;

			double deltaX = swim_step->pos.x() - trackhit->x;
			if(fabs(deltaX) > MAX_HIT_DIST)continue;

			double deltaY = swim_step->pos.y() - trackhit->y;
			if(fabs(deltaY) > MAX_HIT_DIST)continue;

			double delta2 = deltaX*deltaX + deltaY*deltaY + deltaZ*deltaZ;
			
			if(delta2 < min_delta2){
				min_delta2 = delta2;
				min_step_index = i;
			}
		}
		
		// If we didn't find a qualifying step, skip to the next hit
		if(min_step_index == -1)continue;

		// Get the distance of the hit to the R.T.
		TVector3 hit(trackhit->x, trackhit->y, trackhit->z);
		TVector3 delta = GetDistToRT(hit, &rt->swim_steps[min_step_index]);
		
	}
		
	// Delete DReferenceTrajectory object
	delete rt;

	// Create new DTrack object
	DTrack *track = new DTrack;
		
	// Copy in starting values
	track->q = trackcandidate->q;
	track->p = trackcandidate->p;
	track->theta = trackcandidate->theta;
	track->phi = trackcandidate->phi;
	track->x = 0.0;
	track->y = 0.0;
	track->z = trackcandidate->z_vertex;
	track->candidateid = trackcandidate->id;
		
	return track;
}

//------------------
// DTrack_factory::
//------------------
void DTrack_factory::TransformToRTframe(TVector3 &v, swim_step_t *swim_step)
{
	/// Transforms the given vector (replacing its contents) with
	/// values pointing to the same point in the R.T. frame at the
	/// given swim step.

}

//------------------
// GetDistToRT
//------------------
TVector3 DTrack_factory::GetDistToRT(TVector3 &hit, swim_step_t *s2)
{
	/// Calculate the distance of the given hit vector(in the lab
	/// reference frame) to the Reference Trajectory which the
	/// given swim step belongs to. This uses the momentum directions
	/// and positions of the swim_step and it's nearest neighbors
	/// to define a curve and calculate the distance of the hit
	/// from it. The swim step should be the closest one to the hit.
	/// It should also be part of an array of swim steps so that there
	/// both one before and one after it in memory.
	
	// Get the distances of the swim step's nearest neighbors
	// on either side.
	swim_step_t *s1=s2, *s3=s2;
	s1--;	// step before
	s3++;	// step after
	TVector3 delta_s1 = s1->pos - hit;
	TVector3 delta_s3 = s3->pos - hit;
	double delta2_s1 = delta_s1.Mag2();
	double delta2_s3 = delta_s3.Mag2();

	// s_nn is pointer to nearest neighbor of s2 that is closest to hit
	swim_step_t *s_nn = delta2_s1<delta2_s3 ? s1:s3;

	// The point of closest approach on the R.T. for the given hit
	// should now be between the swim steps s1 and s_nn. Each swim
	// step has a momentum vector which is pointing in the
	// same direction as the tangent of the R.T. at that step.
	//
	// We split this up into two 2-D problems: Defining the 
	// parameters of a parabola in the x/y plane and separately
	// in the y/z plane. Note that the values of the hit vector
	// and of the steps are all given in the lab frame. However,
	// here we need to work in the R.T. frame of s2.
	//
	// We assume smooth functions of the form
	//
	// x = Q*y^2 + R*y + S
	//
	// z = T*y^2 + U*y +V
	//
	// Note that since the curve must go through the swim step at
	// s2, the values of S  and V are zero by definition.
	// Also, because the slope of the curve in the x/y plane
	// at the origin is zero, the value of R is also zero.
	//
	// Thus, the only parameter we need to find to define the shape
	// of the R.T. in the X/Y plane is "Q". This is given by:
	//
	// Q = (k1 - k2)/(2*(y1 - y2)) = k1/(2 * y1)
	//
	// where k1,k2 are the slopes of the tangents in the X/Y plane
	// at steps s_nn and s2 respectively. y1 and y2 are the "y"
	// positions of the  steps in the s2 R.T. frame. Note that
	// use of the facts that y2=0, k2=0 were used.
	//
	// Therefore we need to find y1 and k1 in the s2 R.T. frame.
	// k1 comes just from the momentum at s_nn. y1 is just the
	// position at s_nn. Both though are calculated in the s2 R.T.
	// frame.
	//
	// The value of U is given by the slope of the trajectory in
	// the y/z plane at s2. Thus comes just from the momentum
	// at s2.
	//
	// Start by finding the projection of the s_nn momentum
	// vector onto the x/y plane.
	double delta_x = s_nn->mom.Dot(s2->xdir);
	double delta_y = s_nn->mom.Dot(s2->ydir);
	double delta_z = s_nn->mom.Dot(s2->zdir);
	double k1_xy = delta_x/delta_y;
	double k1_zy = delta_z/delta_y;

	// Now find the projection of the s2 momentum vector
	// on the y/z plane to get "e"
	delta_y = s2->mom.Dot(s2->ydir);
	delta_z = s2->mom.Dot(s2->zdir);
	double U = delta_z/delta_y;

	// Find y1 in the R.T. frame at s2
	TVector3 tmpv = s_nn->pos - s2->pos; // shift so origin is that of R.T. at s2
	double y1 = tmpv.Dot(s2->ydir);
	
	double Q = k1_xy/2.0/y1;
	double T = (k1_zy - U)/2.0/y1;

	// Now we need to find the coordinates of the point on the
	// RT closest to the hit. We do this by taking the derivative
	// of the distance squared, setting it equal to zero, and
	// solving for y. This leads to a 3rd order polynomial:
	//
	// Ay^3 + By^2 + Cy + D = 0
	//
	// where:
	//     A = 2(Q^2 + T^2)
	//     B = 3TU
	//     C = U^2 - 2Qx - 2Tz +1
	//     D = -(y + Uz)
	//     x,y,z = coordinates of hit in s2 RT frame
	//
	
	// Transform "hit" from lab coordinates to R.T. coordinates at s2
	tmpv = hit - s2->pos;
	TVector3 hit_rt(tmpv.Dot(s2->xdir), tmpv.Dot(s2->ydir), tmpv.Dot(s2->zdir));
	double xh = hit_rt.x();
	double yh = hit_rt.y();
	double zh = hit_rt.z();
	double A = 2.0*(Q*Q + T*T);
	double B = 3.0*T*U;
	double C = U*U - 2.0*Q*xh - 2.0*T*zh +1.0;
	double D = -(yh + U*zh);
	
	// OK, at this point there is a potential problem. If the trajectory's
	// momentum is in essentially the same direction for both s2 and s_nn,
	// then the above problem is really linear and Q and T are both
	// essentially zero. This means both A and B are zero above and
	// the solution for y is just simply y=-D/C. If the momenta are NOT
	// in the same direction, then we must solve the 3rd order polynomial.
	// The problem is that the first step in solving a 3rd order poly
	// is to divide by the cubic term's coefficient so that you get a new
	// equation whos cubic coeeficient is one. In the case that the
	// cubic coefficient is at or near zero, this blows up the the
	// other coefficients leading to an incorrect result. Thus, we
	// need to decide here whether to approximate the trajectory between
	// s2 and s_nn as a line or a parabola. We do this by checking the
	// angle between the two momenta.
	double y;
	if(s2->mom.Angle(s_nn->mom) <= 1.0E-5){
		// momentum doesn't really change between these two points.
		// use linear approximation
		y = -D/C;
	}else{
	
		// Now we have to solve the 3rd order polynomial for the y value of
		// the point of closest approach on the RT. This is a well documented
		// procedure. Essentially, when you have an equation of the form:
		//
		//  x^3 + a2*x^2 + a1*x + a0 = 0;
		//
		// a change of variables is made such that w = x + a2/3 which leads
		// to a third order poly with no w^2 term:
		//
		//  w^3 + 3.0*b*w + 2*c = 0 
		//
		// where:
		//    b = a1/3 - (a2^2)/9
		//    c = a0/2 - a1*a2/6  + (a2^3)/27
		//
		// The one real root of this is:
		//
		//  w0 = q - p
		//
		// where:
		//    q^3 = d - c
		//    p^3 = d + c
		//    d^2 = b^3 + c^2
	
		double a0 = D/A;
		double a1 = C/A;
		double a2 = B/A;

		double b = a1/3.0 - a2*a2/9.0;
		double c = a0/2.0 - a1*a2/6.0 + a2*a2*a2/27.0;
	
		double d = sqrt(pow(b, 3.0) + pow(c, 2.0));
		double q = pow(d - c, 1.0/3.0);
		double p = pow(d + c, 1.0/3.0);
	
		double w0 = q - p;
		y = w0 - a2/3.0;
	}
		
	// Calculate the y and z coordinates of the RT DOCA point 
	// from the x coordinate
	double x = Q*y*y;
	double z = T*y*y + U*y;
	// Return vector pointing from DOCA point to hit in RT frame
	TVector3 doca_point(x,y,z);

	return hit_rt - doca_point;
}

		
//------------------
// toString
//------------------
const string DTrack_factory::toString(void)
{
	// Ensure our Get method has been called so _data is up to date
	GetNrows();
	if(_data.size()<=0)return string(); // don't print anything if we have no data!

	printheader("row: q:       p:   theta:   phi:    x:    y:    z:");
	
	for(unsigned int i=0; i<_data.size(); i++){

		DTrack *track = _data[i];

		printnewrow();
		
		printcol("%x", i);
		printcol("%+d", (int)track->q);
		printcol("%3.3f", track->p);
		printcol("%1.3f", track->theta);
		printcol("%1.3f", track->phi);
		printcol("%2.2f", track->x);
		printcol("%2.2f", track->y);
		printcol("%2.2f", track->z);

		printrow();
	}
	
	return _table;
}




//======================================================================
//======================================================================
//           UNUSED CODEBELOW HERE
//======================================================================
//======================================================================

#if 0

#define USE_MINUIT 0

#if USE_MINUIT
// This is for MINUIT which is really not suited for a multi-threaded
// event processing environment. We do a few backflips here ...
static void FCN(int &npar, double *derivatives, double &chisq, double *par, int iflag);

static pthread_mutex_t trk_mutex = PTHREAD_MUTEX_INITIALIZER;
static JEventLoop *g_loop = NULL;
static const DMagneticFieldMap *g_bfield = NULL;
static JFactory<DTrackHit> *g_fac_trackhit;
static const DTrackCandidate *g_trackcandidate;
static double min_chisq, min_par[6];
static TMinuit *minuit=NULL;

typedef struct{
	double x,y,z;
	double dxdz, dydz;
	double d2xdz2, d2ydz2;
}trk_step_t;

#endif // USE_MINUIT

//------------------
// init
//------------------
jerror_t DTrack_factory::init(void)
{
#if USE_MINUIT
	pthread_mutex_lock(&trk_mutex);
	if(minuit==NULL){
		int icondn;
		minuit = new TMinuit();
		minuit->mninit(0,0,0);
		minuit->SetFCN(FCN);
		minuit->mncomd("CLEAR", icondn);
		minuit->mncomd("SET PRI -1", icondn);
		minuit->mncomd("SET NOWARNINGS", icondn);
		minuit->mncomd("SET BATCH", icondn);
		minuit->mncomd("SET STRATEGY 2", icondn);
		
		minuit->mncomd("SET PRI -1", icondn);
		minuit->mncomd("SET NOWARNINGS", icondn);
		minuit->mncomd("SET BATCH", icondn);
		minuit->mncomd("SET STRATEGY 2", icondn);

		minuit->DefineParameter(0,"p",		0.0, 0.02, 0.0, 0.0);
		minuit->DefineParameter(1,"theta",	0.0, 0.17, 0.0, 0.0);
		minuit->DefineParameter(2,"phi",		0.0, 0.17, 0.0, 0.0);
		minuit->DefineParameter(3,"x",		0.0, 0.50, 0.0, 0.0);
		minuit->DefineParameter(4,"y",		0.0, 0.50, 0.0, 0.0);
		minuit->DefineParameter(5,"z",		0.0, 7.0, 0.0, 0.0);

		minuit->mncomd("FIX 4", icondn);
		minuit->mncomd("FIX 5", icondn);
		minuit->mncomd("FIX 6", icondn);
	}
	pthread_mutex_unlock(&trk_mutex);

#endif // USE_MINUIT

	return NOERROR;
}

//------------------
// brun
//------------------
jerror_t DTrack_factory::brun(JEventLoop *loop, int runnumber)
{
	//dgeom = loop->GetJApplication()->GetJGeometry(runnumber);
	//bfield = NULL;
	//bfield = dgeom->GetDMagneticFieldMap();
	
	bfield = new DMagneticFieldMap(); // temporary until new geometry scheme is worked out
	
	return NOERROR;
}

//------------------
// evnt
//------------------
jerror_t DTrack_factory::evnt(JEventLoop *loop, int eventnumber)
{
	// Get the track candidates and hits
	vector<const DTrackCandidate*> trackcandidates;
	vector<const DTrackHit*> trackhits;
	loop->Get(trackcandidates);
#if USE_MINUIT
	JFactory<DTrackHit> *fac_trackhit =
#endif //USE_MINUIT
	loop->Get(trackhits,"MC");
		
	for(unsigned int i=0; i<trackcandidates.size(); i++){
		const DTrackCandidate *trackcandidate = trackcandidates[i];
		DTrack *track = new DTrack;
		
		// Copy in starting values
		track->q = trackcandidate->q;
		track->p = trackcandidate->p;
		track->theta = trackcandidate->theta;
		track->phi = trackcandidate->phi;
		track->x = 0.0;
		track->y = 0.0;
		track->z = trackcandidate->z_vertex;
		track->candidateid = trackcandidate->id;
		
#if USE_MINUIT
		// Lock mutex so only one thread at a time accesses minuit
		pthread_mutex_lock(&trk_mutex);
		g_loop = loop;
		g_bfield = bfield;
		g_fac_trackhit = fac_trackhit;
		g_trackcandidate = trackcandidate;
		
		// Calculate chi-square value
		double par[6];
		par[0] = track->p;
		par[1] = track->theta;
		par[2] = track->phi;
		par[3] = track->x;
		par[4] = track->y;
		par[5] = track->z;
		int npar=6;
		double dfdpar[6];
		double chisq=0.0;;
		int iflag = 2;
		FCN(npar, dfdpar, chisq, par, iflag);
		
		// Set the initial parameter values
		for(int i=0;i<6;i++){
			int icondn;
			char cmd[256];
			sprintf(cmd, "SET PAR %d %f", i+1, par[i]);
			min_par[i] = par[i];
			minuit->mncomd(cmd, icondn);
		}		
		
		min_chisq = 1000.0;
		//minuit->Migrad();
		//minuit->mncomd("MINIMIZE 100", icondn);
		//int icondn;
		minuit->mncomd("SEEK 100 1", icondn);

//cout<<__FILE__<<":"<<__LINE__<<" initial:"<<chisq<<"  final:"<<min_chisq<<endl;
//for(int i=0; i<6; i++){
//	cout<<__FILE__<<":"<<__LINE__<<" -- p"<<i<<"  initial:"<<par[i]<<"  final:"<<min_par[i]<<endl;
//}

		// Unlock mutex so only one thread at a time accesses minuit
		pthread_mutex_unlock(&trk_mutex);
		
		if(min_chisq < chisq){
			track->p = min_par[0];
			track->theta = min_par[1];
			track->phi = min_par[2];
			track->x = min_par[3];
			track->y = min_par[4];
			track->z = min_par[5];
		}

#endif // USE_MINUIT

		_data.push_back(track);
	}

	return NOERROR;
}


#if USE_MINUIT
//------------------
// FCN   -- for MINUIT
//------------------
void FCN(int &npar, double *derivatives, double &chisq, double *par, int iflag)
{
	// Step a singley charged particle with the momentum and vertex
	// position given by par through the magnetic field, calculating
	// the chisq from the closest hits. We need to keep track of the
	// positions and first and second derivatives of each step along
	// the trajectory. We will first step through field from vertex 
	// until we either hit the approximate BCAL, UPV, or FCAL positions.

	// Minuit will sometimes try zero momentum tracks. Immediately return
	// a large ChiSq when that happens
	if(par[0] < 0.050){chisq=1000.0; return;}
	
	vector<double> chisqv;
	TVector3 p;
	p.SetMagThetaPhi(par[0],par[1],par[2]);
	TVector3 pos(par[3],par[4],par[5]);
	DMagneticFieldStepper *stepper = new DMagneticFieldStepper(g_bfield, g_trackcandidate->q, &pos, &p);
	stepper->SetStepSize(0.1);
	
	// Step until we hit a boundary
	vector<trk_step_t> trk_steps;
	for(int i=0; i<10000; i++){
		stepper->Step(&pos);
		trk_step_t trk_step;
		trk_step.x = pos.X();
		trk_step.y = pos.Y();
		trk_step.z = pos.Z();
		trk_steps.push_back(trk_step);
		
		//if(pos.Perp()>65.0){cout<<__FILE__<<":"<<__LINE__<<" hit BCAL"<<endl;break;} // ran into BCAL
		//if(pos.Z()>650.0){cout<<__FILE__<<":"<<__LINE__<<" hit FCAL"<<endl;break;} // ran into FCAL
		//if(pos.Z()<-50.0){cout<<__FILE__<<":"<<__LINE__<<" hit UPV"<<endl;break;} // ran into UPV
		if(pos.Perp()>65.0){break;} // ran into BCAL
		if(pos.Z()>650.0){break;} // ran into FCAL
		if(pos.Z()<-50.0){break;} // ran into UPV
	}
	delete stepper;

#if 0
	// Calculate first derivatives in trk_steps
	for(unsigned int i=1; i<trk_steps.size()-1; i++){
		// We do this by averaging the slopes of the lines connecting
		// this point to the previous and nexxt ones. We assume
		// that the stepper used small enough steps to vary the
		// positions smoothly.
		trk_step_t *prev = &trk_steps[i-1];
		trk_step_t *curr = &trk_steps[i];
		trk_step_t *next = &trk_steps[i+1];
		curr->dxdz = ((curr->x-prev->x)/(curr->z-prev->z) + (next->x-curr->x)/(next->z-curr->z))/2.0;
		curr->dydz = ((curr->y-prev->y)/(curr->z-prev->z) + (next->y-curr->y)/(next->z-curr->z))/2.0;
	}

	// Calculate second derivatives in trk_steps
	for(unsigned int i=2; i<trk_steps.size()-2; i++){
		// We do this by averaging the slopes of the lines connecting
		// this point to the previous and nexxt ones. We assume
		// that the stepper used small enough steps to vary the
		// positions smoothly.
		trk_step_t *prev = &trk_steps[i-1];
		trk_step_t *curr = &trk_steps[i];
		trk_step_t *next = &trk_steps[i+1];
		curr->d2xdz2 = ((curr->dxdz-prev->dxdz)/(curr->z-prev->z) + (next->dxdz-curr->dxdz)/(next->z-curr->z))/2.0;
		curr->d2ydz2 = ((curr->dydz-prev->dydz)/(curr->z-prev->z) + (next->dydz-curr->dydz)/(next->z-curr->z))/2.0;
	}
#endif //0 
	
	// Loop over the track hits for this candidate using
	// the closest trk_step to determine the distance
	// from the hit to the track.
	vector<const DTrackHit*> trackhits;
	vector<void*>& vptrs = g_fac_trackhit->Get();
	for(unsigned int i=0; i<vptrs.size(); i++){
		const DTrackHit *trackhit = (const DTrackHit*)vptrs[i];
		if(!trackhit)continue;
		
		// Ignore hits outside of the CDC and FDC for now
		if(trackhit->z <= 0.0)continue;
		if(trackhit->z > 405.0)continue;
		double r = sqrt(pow((double)trackhit->x,2.0) + pow((double)trackhit->y,2.0));
		if(r >= 65.0)continue;
		
		// We don't really want the hit closest physically in space.
		// Rather, we want the distance measured in standard deviations
		// of detector resolutions. This depends upon the detector type
		// from which the hit came. Note that this also depends on the
		// step size used in the stepper. The FDC in particular needs
		// much larger sigma values due to detector resolution alone.
		double sigmaxy;
		double sigmaz;
		switch(trackhit->system){
			case SYS_CDC:	sigmaxy = 0.80;	sigmaz = 7.50;	break;
			case SYS_FDC:	sigmaxy = 0.50;	sigmaz = 0.50;	break;
			default:			sigmaxy = 10.00;	sigmaz = 10.00;
		}

		double xh = trackhit->x;
		double yh = trackhit->y;
		double zh = trackhit->z;

		// Loop over all steps to find the closest one to this hit
		double r2_min = 1.0E6;
		trk_step_t *trk_step_min=NULL;
		for(unsigned int j=0; j<trk_steps.size(); j++){
			trk_step_t *trk_step = &trk_steps[j];
			double xs = trk_step->x;
			double ys = trk_step->y;
			double zs = trk_step->z;
			
			double r2 = pow((xh-xs)/sigmaxy,2.0) + pow((yh-ys)/sigmaxy,2.0) + pow((zh-zs)/sigmaz,2.0);
//cout<<__FILE__<<":"<<__LINE__<<" dx="<<xh-xs<<" dy="<<yh-ys<<" dz="<<zh-zs<<"   z="<<zs<<endl;
			if(r2<r2_min){
				r2_min = r2;
				trk_step_min = trk_step;
			}
		}
//cout<<__FILE__<<":"<<__LINE__<<" ##### r_min= "<<sqrt(r2_min)<<"  z="<<trk_step_min->z<<endl;
		
		// Now calculate the distance from the track to the hit
		double d = sqrt(r2_min); // use the distance to the step for now
		
		// If the hit and the step are more than 3 sigma apart, then
		// don't include this hit in the ChiSq.
		if(d>3.0)continue;

		// Add this hit to the chisq vector
		chisqv.push_back(d*d);
	}

	// Sum up total chisq/dof for this event
	chisq = 0.0;
	for(unsigned int i=0; i<chisqv.size(); i++)chisq += chisqv[i];
	chisq /= (double)chisqv.size() - 4.0;

//chisq = pow((par[0]-1.0)/0.100,2.0) + pow((par[5]-65.0)/7.0,2.0);

	// If NO hits were close to this track (how does that even happen?)
	// then the size of chisqv will be zero and chisq will now be
	// "nan". In this case, set the chisq to a very large value to
	// indicate a bad fit.
	if(chisqv.size() < 5)chisq=1000.0;

	if(chisq<min_chisq){
		min_chisq = chisq;
		for(int i=0;i<6;i++)min_par[i] = par[i];
	}

//cout<<__FILE__<<":"<<__LINE__<<" chisq= "<<chisq<<"  chisqv.size()="<<chisqv.size()<<endl;
}
#endif // USE_MINUIT

#endif // 0
