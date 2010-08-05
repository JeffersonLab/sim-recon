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
	max_bin_valid = false;
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

	for(unsigned int i=0; i<points.size(); i++){
		const DVector2 &point = points[i];
		DVector2 g(point.Y(), -point.X()); // perp. to point
		g /= g.Mod(); // Make unit vector
		DVector2 pos = point/2.0; // initialize to a point on the line
		
		// Find intersection with an edge of the histogram area that
		// we can use as a starting point.
		double beta = FindBeta(xmin, ymin, xmax-xmin, ymax-ymin, pos, g);
		pos += beta*g;
		
		// Find the indexes of the bin we're about to step across.
		// Note: we must also figure out if we need to step
		// in the +g or -g direction to go "into" the histogram. Since
		// we'll keep stepping in this direction, we adjust the sign
		// of g as needed to always step in the scan direction
		int ix, iy; // bin indexes
		FindIndexes(pos + small_step*g, ix, iy);
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
			beta = FindBeta(xmin+(double)ix*bin_widthx, ymin+(double)iy*bin_widthy, bin_widthx, bin_widthy, pos, g);
			
			// Beta too large indicates problem
			if(beta*beta > (bin_widthx*bin_widthx + bin_widthy*bin_widthy))break;
			
			// increment histo for bin we just stepped across
			if(ix<0 || ix>=(int)Nbinsx || iy<0 || iy>=(int)Nbinsy)break; // must have left the histo
			int index = ix + iy*Nbinsy;
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
	}

	return GetMaxBinLocation();
}

//---------------------------------
// Fill
//---------------------------------
void DHoughFind::Fill(double x, double sigmax, double y, double sigmay)
{
	/// Increment the histogram bins according to the given location
	/// and sigmas. This will increment each bin by calculating the
	/// product of 2 gaussians using the coordinates of the center of
	/// the bin. Only bins within 4 sigma in each dimension are incremented.
	///
	/// Note that the histogram is NOT reset prior to filling. This is
	/// to allow accumulation over multiple calls.
	
	int ixmin = (int)(((x-xmin) - 4.0*sigmax)/bin_widthx -0.5);
	int ixmax = (int)(((x-xmin) + 4.0*sigmax)/bin_widthx +0.5);
	int iymin = (int)(((y-ymin) - 4.0*sigmay)/bin_widthy -0.5);
	int iymax = (int)(((y-ymin) + 4.0*sigmay)/bin_widthy +0.5);
	if(ixmin<0)ixmin=0;
	if(iymin<0)iymin=0;
	if(ixmax>(int)Nbinsx)ixmax=Nbinsx;
	if(iymax>(int)Nbinsy)iymax=Nbinsy;
	
	// Loop over bins
	double x_bin = xmin + (0.5+(double)ixmin)*bin_widthx;
	for(int i=ixmin; i<ixmax; i++, x_bin+=bin_widthx){
		double y_bin = ymin + (0.5+(double)iymin)*bin_widthy;
		for(int j=iymin; j<iymax; j++, y_bin+=bin_widthy){
			double k_x = (x - x_bin)/sigmax;
			double k_y = (y - y_bin)/sigmay;
			double gauss_x = exp(k_x*k_x)/sigmax; // divide by sigma to make constant area
			double gauss_y = exp(k_y*k_y)/sigmay; // divide by sigma to make constant area
			
			int index = i + j*Nbinsy;
			hist[index] += gauss_x*gauss_y;
		}
	}
}

//---------------------------------
// GetMaxBinLocation (static)
//---------------------------------
DVector2 DHoughFind::GetMaxBinLocation(vector<const DHoughFind*> &houghs)
{
	/// This routine is designed to be called statically via:
	///
	///     DHoughFind::GetMaxBinLocation(houghs);
	///
	/// It will find the location of the center of the bin with the maximum
	/// content based on the sum of the input DHough objects. It does this
	/// dynamically without maintaining a sum histogram.
	///
	/// WARNING: If you call this as a method of an existing object, (e.g.
	/// like this myhough->GetMacBinLocation(houghs); ) the object is NOT
	/// used unless it explicitly appears in the "houghs" list!
	///
	/// It is left to the caller to ensure the limits and number of bins for each
	/// DHough object are the same.

	if(houghs.size()<1)return DVector2(0.0, 0.0);
	
	unsigned int Nbinsx = houghs[0]->Nbinsx;
	unsigned int Nbinsy = houghs[0]->Nbinsy;
	
	unsigned int imax_binx=0, imax_biny=0;
	double max_bin_content = 0.0;
	for(unsigned int i=0; i<Nbinsx; i++){
		for(unsigned int j=0; j<Nbinsy; j++){
			double tot = 0.0;
			for(unsigned int k=0; k<houghs.size(); k++){
				unsigned int index = i + Nbinsy*j;
				tot += houghs[k]->hist[index];
			}
			if(tot>max_bin_content){
				max_bin_content = tot;
				imax_binx = i;
				imax_biny = j;
			}
		}
	}
	
	double x = houghs[0]->xmin + (0.5+(double)imax_binx)*houghs[0]->bin_widthx;
	double y = houghs[0]->ymin + (0.5+(double)imax_biny)*houghs[0]->bin_widthy;
	
	return DVector2(x, y);
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

