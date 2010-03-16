// $Id$
//
//    File: DHoughFind.cc
// Created: Wed Oct 31 05:59:26 EDT 2007
// Creator: davidl (on Darwin Amelia.local 8.10.1 i386)
//

#include <algorithm>
#include <iostream>
#include <cmath>
using namespace std;

#include "DHoughFind.h"

//---------------------------------
// DHoughFind    (Constructor)
//---------------------------------
DHoughFind::DHoughFind()
{
	hist = NULL;
}

//---------------------------------
// DHoughFind    (Constructor)
//---------------------------------
DHoughFind::DHoughFind(double xmin, double xmax, double ymin, double ymax, unsigned int Nbinsx, unsigned int Nbinsy)
{
	hist = NULL;
	SetLimits(xmin, xmax, ymin, ymax, Nbinsx, Nbinsy);
}

//---------------------------------
// ~DHoughFind    (Destructor)
//---------------------------------
DHoughFind::~DHoughFind()
{
	if(hist)delete[] hist;
}

//---------------------------------
// SetLimits
//---------------------------------
void DHoughFind::SetLimits(double xmin, double xmax, double ymin, double ymax, unsigned int Nbinsx, unsigned int Nbinsy)
{
	// If hist already has memory allocated, check if we need to reallocate
	if(hist){
		if(this->Nbinsx*this->Nbinsy != Nbinsx*Nbinsy){
			delete[] hist;
			hist = NULL;
		}
	}
	
	// Allocate memory if necessary
	if(!hist)hist = new double[Nbinsx*Nbinsy];

	this->xmin = xmin;
	this->xmax = xmax;
	this->ymin = ymin;
	this->ymax = ymax;
	this->Nbinsx = Nbinsx;
	this->Nbinsy = Nbinsy;
	
	bin_widthx = (xmax-xmin)/(double)(Nbinsx-1);
	bin_widthy = (ymax-ymin)/(double)(Nbinsy-1);
	bin_size = bin_widthx<bin_widthy ? bin_widthx:bin_widthy; // used for characteristic "size" of a bin
	
	ResetHist();
}

//---------------------------------
// ResetHist
//---------------------------------
void DHoughFind::ResetHist(void)
{
	for(unsigned int i=0; i<Nbinsx*Nbinsy; i++)hist[i]=0.0;
	//memset(hist, 0, Nbinsx*Nbinsy*sizeof(double));

	imax_binx = Nbinsx/2;
	imax_biny = Nbinsy/2;
	max_bin_content = 0.0;
}

//---------------------------------
// GetMaxBinLocation
//---------------------------------
DVector2 DHoughFind::GetMaxBinLocation(void)
{
	double x = xmin + (0.5+(double)imax_binx)*bin_widthx;
	double y = ymin + (0.5+(double)imax_biny)*bin_widthy;

	DVector2 pos(x,y);
	
	return pos;
}

//---------------------------------
// GetMaxBinContent
//---------------------------------
double DHoughFind::GetMaxBinContent(void)
{
	return max_bin_content;
}

//---------------------------------
// GetSigmaX
//---------------------------------
double DHoughFind::GetSigmaX(void)
{
	return bin_widthx/sqrt(12.0);
}

//---------------------------------
// GetSigmaY
//---------------------------------
double DHoughFind::GetSigmaY(void)
{
	return bin_widthy/sqrt(12.0);
}

//---------------------------------
// Find
//---------------------------------
DVector2 DHoughFind::Find(void)
{
	return Find(points);
}

//---------------------------------
// Find
//---------------------------------
DVector2 DHoughFind::Find(const vector<DVector2> &points)
{
	/// Loop over "points" transforming them into lines and filling the 2-D
	/// histogram.

		
	// When determining the bins we just stepped into and out of,
	// we add a small step in the g direction to push us just over
	// the boundary we're currently on. We calculate the size of
	// that step here so we don't have to over and over inside the loop
	double small_step = bin_size*1.0E-3;

//_DBG_<<"points.size()="<<points.size()<<endl;
	for(unsigned int i=0; i<points.size(); i++){
		const DVector2 &point = points[i];
		DVector2 g(point.Y(), -point.X()); // perp. to point
		g /= g.Mod(); // Make unit vector
		DVector2 pos = point/2.0; // initialize to a point on the line
		
		// Find intersection with an edge of the histogram area that
		// we can use as a starting point.
//_DBG_<<"============= point "<<i<<"==============="<<endl;
		double beta = FindBeta(xmin, ymin, xmax-xmin, ymax-ymin, pos, g);
		pos += beta*g;
		
		// Find the indexes of the bin we're about to step across.
		// Note: we must also figure out if we need to step
		// in the +g or -g direction to go "into" the histogram. Since
		// we'll keep stepping in this direction, we adjust the sign
		// of g as needed to always step in the scan direction
		int ix, iy; // bin indexes
		FindIndexes(pos + small_step*g, ix, iy);
//_DBG_<<"beta="<<beta<<" ix="<<ix<<" iy="<<iy<<endl;
		if(ix<0 || ix>=(int)Nbinsx-1 || iy<0 || iy>=(int)Nbinsy-1){
			g *= -1.0;
			FindIndexes(pos + small_step*g, ix, iy);

			// If we're still out of range, then this line doesn't intersect our histo ever
			if(ix<0 || ix>=(int)Nbinsx-1 || iy<0 || iy>=(int)Nbinsy-1){
				continue;
			}
		}
		
		unsigned int Niterations=0;
		do{
			
			// Find distance to boundary of next bin
//_DBG_<<"   --- iteration "<<Niterations<<" ---"<<endl;
			beta = FindBeta(xmin+(double)ix*bin_widthx, ymin+(double)iy*bin_widthy, bin_widthx, bin_widthy, pos, g);
			
			// Beta too large indicates problem
			if(beta*beta > (bin_widthx*bin_widthx + bin_widthy*bin_widthy))break;
			
			// increment histo for bin we just stepped across
			if(ix<0 || ix>=(int)Nbinsx || iy<0 || iy>=(int)Nbinsy)break; // must have left the histo
			int index = ix + iy*Nbinsy;
//_DBG_<<"before: index="<<index<<" ix="<<ix<<" iy="<<iy<<endl;
			hist[index] += fabs(beta);
			if(hist[index]>max_bin_content){
				max_bin_content = hist[index];
				imax_binx = ix;
				imax_biny = iy;
			}
			
			// Step to next boundary and find indexes of next bin
			pos += beta*g;
			FindIndexes(pos + small_step*g, ix, iy);

		}while(++Niterations<2*Nbinsx);
//_DBG_<<"after:  ix="<<ix<<" iy="<<iy<<" Niterations="<<Niterations<<endl;
	}

	return GetMaxBinLocation();
}

//---------------------------------
// AddPoint
//---------------------------------
void DHoughFind::AddPoint(const DVector2 &point)
{
	points.push_back(point);
}

//---------------------------------
// AddPoint
//---------------------------------
void DHoughFind::AddPoint(const double &x, const double &y)
{
	DVector2 point(x,y);
	AddPoint(point);
}

//---------------------------------
// AddPoints
//---------------------------------
void DHoughFind::AddPoints(const vector<DVector2> &points)
{
	for(unsigned int i=0; i<points.size(); i++)this->points.push_back(points[i]);
}

//---------------------------------
// ClearPoints
//---------------------------------
void DHoughFind::ClearPoints(void)
{
	points.clear();
}

#if 0
//---------------------------------
// FindIndexes
//---------------------------------
void DHoughFind::FindIndexes(const DVector2 &pos, int &ix, int &iy)
{
	ix = (int)floor((pos.X()-xmin)/bin_widthx);
	iy = (int)floor((pos.Y()-ymin)/bin_widthy);
}

//---------------------------------
// FindBeta
//---------------------------------
double DHoughFind::FindBeta(double xlo, double ylo, double widthx, double widthy, DVector2 &pos, DVector2 &step)
{
	DVector2 a0(xlo, ylo);
	DVector2 xdir(1.0, 0.0);
	DVector2 ydir(0.0, 1.0);
	DVector2 stepdir = step/step.Mod();

	//vector<double> beta(4);
	double beta[4];
	beta[0] = FindBeta(a0, xdir, pos, stepdir);
	beta[1] = FindBeta(a0+widthx*xdir, ydir, pos, stepdir);
	beta[2] = FindBeta(a0+widthx*xdir+widthy*ydir, -1.0*xdir, pos, stepdir);
	beta[3] = FindBeta(a0+widthy*ydir, -1.0*ydir, pos, stepdir);
	//beta.push_back(FindBeta(a0, xdir, pos, stepdir));
	//beta.push_back(FindBeta(a0+widthx*xdir, ydir, pos, stepdir));
	//beta.push_back(FindBeta(a0+widthx*xdir+widthy*ydir, -1.0*xdir, pos, stepdir));
	//beta.push_back(FindBeta(a0+widthy*ydir, -1.0*ydir, pos, stepdir));

	// Here, we want to choose the closest boundary of the bin which is not
	// the boundary the point "pos" lies on. The values of beta give the
	// distances to the lines which define the boundaries, but not clipped
	// to the bin itself.

	// Loop over beta values
	DVector2 bin_center(xlo+widthx/2.0, ylo+widthy/2.0);
	double min_dist_to_center = 1.0E6;
	double min_beta = 1.0E6;
	for(unsigned int i=0; i<4; i++){
		if(fabs(beta[i])<1.0E-6*bin_size)continue; // This is most likely due to our starting pos being on a boundary already. Ignore it.
		
		DVector2 delta_center = pos + beta[i]*stepdir - bin_center;
		if(delta_center.Mod() < min_dist_to_center){
			min_dist_to_center = delta_center.Mod();
			min_beta = beta[i];
		}
	}
	
	return min_beta;
}

//---------------------------------
// FindBeta
//---------------------------------
double DHoughFind::FindBeta(const DVector2 &a, const DVector2 &b, const DVector2 &c, const DVector2 &d)
{
	/// Given the 2-D vectors a,b,c, and d find the scaler value "beta" such that
	///
	/// a + alpha*b = c + beta*d
	///
	/// It vectors b and d should be unit vectors.
	/// The value of alpha is not returned, but can be calculated via:
	///
	/// alpha = beta*(d.b) - b.(a-c)
	///
	/// If the vectors b and d are parallel, then the return value will be inf
	/// or whatever division by zero returns on the system.
#if 0
	double ex = a.X() - c.X();
	double ey = a.Y() - c.Y();
	double bx = b.X();
	double by = b.Y();
	double dx = d.X();
	double dy = d.Y();
	double b_dot_d = bx*dx + by*dy;
	double d_dot_e = dx*ex + dy*ey;
	double b_dot_e = bx*ex + by*ey;
	return (d_dot_e - b_dot_d*b_dot_e)/(1.0-b_dot_d*b_dot_d);
#endif

	return (d*(a-c) - (d*b)*(b*(a-c)))/(1.0-pow(b*d, 2.0));
}
#endif

//---------------------------------
// PrintHist
//---------------------------------
void DHoughFind::PrintHist(void)
{
	/// Dump an ASCII representation of the normalized histogram
	/// to the screen. This will dump a table of Nbinsx x Nbinsy
	/// characters to the screen as a sort of cheap visual of
	/// the histgram. It is only intended for debugging.
	
	// X-axis labels
	string space(Nbinsx-10, ' ');
	cout<<xmin<<space<<xmax<<endl;
	
	
	for(unsigned int i=0; i<Nbinsy; i++){
		string row(Nbinsx,' ');
		for(unsigned int j=0; j<Nbinsx; j++){

			int index = j + i*Nbinsx;
			unsigned int d = (unsigned int)(floor(15.0*hist[index]/max_bin_content));
			char c[16];
			sprintf(c, "%x", d);
			row[j] = c[0] == '0' ? '.':c[0];
		}
		cout<<row<<" "<<ymin+(double)i*bin_widthy<<endl;
	}
}

