
#include <iostream>
using namespace std;

#include <math.h>

#include "cdc_fithelix.h"
#include <TMinuit.h>
#include <TEllipse.h>
#include <TCanvas.h>

static TVector3 *DATAV=NULL;
static int NDATAV = 0;
static TEllipse *ellipse=NULL;

//-----------------
// cdc_fithelix
//-----------------
derror_t cdc_fithelix(TVector3 *v, int Npoints)
{
	float p,theta, phi, q;
	
	cdc_firstguess(v, Npoints, theta, phi, p, q);
	cout<<"theta="<<theta<<"  phi="<<phi<<"  q="<<q<<"  Npoints="<<Npoints<<endl;


#if 0
	DATAV = v;
	NDATAV = Npoints;

	// Create a TMinuit object to do the fit
	TMinuit minuit(3);
	minuit.SetFCN(HelixChisq);
	
	Double_t arglist[10];
	Int_t ierflg=0;
	
	arglist[0] = 1;
	minuit.mnexcm("SET ERR", arglist, 1, ierflg);
	
	// Starting parameters and step sizes
	// par0=p  par1=phi  par2=theta
	Double_t vstart[3]	= { 3.0, 0.0, 0.2 };
	Double_t step[3]		= { 0.5, 0.2, 0.1 };
	minuit.mnparm(0, "p", vstart[0], step[0], 0,0,ierflg);
	
	// Mimimize
	arglist[0] = 500;
	arglist[1] = 1.0;
	minuit.mnexcm("MIGRAD", arglist, 2, ierflg);
	
	// Print results
	Double_t amin,edm,errdef;
	Int_t nvpar,nparx,icstat;
	minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);
	minuit.mnprin(3,amin);
#endif

	return NOERROR;
}

//-----------------
// cdc_firstguess
//-----------------
derror_t cdc_firstguess(TVector3 *points, int Npoints, float &theta
	,float &phi, float &p, float &q)
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
	TVector3 *v=points;
	for(int i=0;i<Npoints;i++, v++){
		float x=v->x();
		float y=v->y();
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
	float B=-2.0;
	q = +1.0;
	float r0 = sqrt(x0*x0 + y0*y0);
	float hbarc = 197.326;
	float p_trans = q*B*r0/hbarc*.61; // are these the right units?
	
	// Assuming points are ordered in increasing z, the sign of the
	// cross-product between successive points will be the opposite
	// sign of the charge. Since it's possible to have a few "bad"
	// points, we don't want to rely on any one to determine this.
	// The method we use is to sum cross-products of the first and
	// middle points, the next-to-first and next-to-middle, etc.
	float xprod_sum = 0.0;
	int n_2 = Npoints/2; 
	v = points;
	TVector3 *v2=&points[n_2];
	for(int i=0;i<n_2;i++, v++, v2++){
		xprod_sum += v->x()*v2->y() - v2->x()*v->y();
	}
	if(xprod_sum>0.0)q = -q;

	// Phi is pi/2 out of phase with x0,y0. The sign of the phase difference
	// depends on the charge
	phi = atan2(y0,x0);
	phi += q>0.0 ? M_PI_2:-M_PI_2;
	
	// Theta is determined by extrapolating the helix back to the target.
	// To do this, we need dphi/dz and a phi,z point. The easiest way to
	// get these is by a simple average (I guess).
	v = points;
	v2=&v[1];
	float dphidz =0.0;
	int Ndphidzpoints = 0;
	for(int i=0;i<Npoints-1;i++, v++, v2++){
		float myphi1 = atan2(v->y()-y0,  v->x()-x0);
		float myphi2 = atan2(v2->y()-y0, v2->x()-x0);
		float mydphidz = (myphi2-myphi1)/(v2->z()-v->z());
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
	
	if(ellipse)delete ellipse;
	ellipse = new TEllipse(x0, y0, r0, r0);
	ellipse->Draw();

	float px = p_trans*cos(phi);
	float py = p_trans*sin(phi);
	float pz = p*cos(theta);
cout<<endl;
cout<<__FILE__<<":"<<__LINE__<<" px="<<px<<"  py="<<py<<"  pz="<<pz<<"  E="<<sqrt(p*p+0.14*0.14)<<endl;
cout<<endl;

	return NOERROR;
}

//-----------------
// HelixChisq
//-----------------
void HelixChisq(Int_t& npar, Double_t *x, Double_t &chisq, Double_t *par, Int_t iflag)
{
	float x0;

	chisq = 0.0;
	for(int i=0;i<NDATAV;i++){
		
	}

	chisq = 1.0;
}

